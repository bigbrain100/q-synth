// Q-SYNTH library: synthetic Q participants + factor recovery.
// Self-contained: correlation/PCA/varimax/factorScores are open implementations in
// engine.mjs, and centroid is in centroid.mjs. Both extraction paths are validated
// against a real KADE export (PCA bit-exact; centroid structure col-corr >=0.9995;
// see validation/).
import { correlationMatrix, pcaExtract, varimax, autoFlag } from './engine.mjs';
import { factorScores, generateCompositeSort } from './engine.mjs';
import { centroidExtract } from './centroid.mjs';

// ---------- 可重現亂數 ----------
export function makeRng(seed) {
  let s = seed >>> 0;
  const u = () => { s = (s + 0x6D2B79F5) | 0; let t = Math.imul(s ^ (s >>> 15), 1 | s); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
  const n = () => { let a = 0, b = 0; while (a === 0) a = u(); while (b === 0) b = u(); return Math.sqrt(-2 * Math.log(a)) * Math.cos(2 * Math.PI * b); };
  return { u, n };
}

// ---------- 基本工具 ----------
export function zscore(a) {
  const m = a.reduce((s, x) => s + x, 0) / a.length;
  const sd = Math.sqrt(a.reduce((s, x) => s + (x - m) * (x - m), 0) / (a.length - 1)) || 1;
  return a.map(x => (x - m) / sd);
}
export function tuckerPhi(a, b) {
  let dot = 0, na = 0, nb = 0;
  for (let i = 0; i < a.length; i++) { dot += a[i] * b[i]; na += a[i] * a[i]; nb += b[i] * b[i]; }
  return dot / (Math.sqrt(na * nb) || 1);
}
// 總和 = M、約鐘形、每欄 ≥1 的強制分佈(-4..+4)
export function buildDistribution(M) {
  const positions = [-4, -3, -2, -1, 0, 1, 2, 3, 4];
  const w = positions.map(p => Math.exp(-(p * p) / (2 * 2.2 * 2.2)));
  const wsum = w.reduce((a, b) => a + b, 0);
  const counts = w.map(x => Math.max(1, Math.round((x / wsum) * M)));
  let diff = M - counts.reduce((a, b) => a + b, 0);
  const order = [4, 3, 5, 2, 6, 1, 7, 0, 8]; let oi = 0;
  while (diff !== 0 && oi < 10000) { const idx = order[oi % order.length]; if (diff > 0) { counts[idx]++; diff--; } else if (counts[idx] > 1) { counts[idx]--; diff++; } oi++; }
  return counts;
}
// 複合對稱相關矩陣 (diag 1, off-diag rho) 的 Cholesky 下三角
export function cholCS(K, rho) {
  const R = Array.from({ length: K }, (_, i) => Array.from({ length: K }, (_, j) => (i === j ? 1 : rho)));
  const L = Array.from({ length: K }, () => new Array(K).fill(0));
  for (let i = 0; i < K; i++) for (let j = 0; j <= i; j++) {
    let s = 0; for (let k = 0; k < j; k++) s += L[i][k] * L[j][k];
    L[i][j] = i === j ? Math.sqrt(Math.max(1e-12, R[i][i] - s)) : (R[i][j] - s) / L[j][j];
  }
  return L;
}
const clamp = (x, lo, hi) => Math.max(lo, Math.min(hi, x));
// 依模式把連續向量渲染成排序:'forced'=壓格(ipsative);'free'=round/clamp(非 ipsative)
export function render(vec, mode, distribution) {
  return mode === 'forced' ? generateCompositeSort(vec, distribution) : vec.map(v => clamp(Math.round(v), -4, 4));
}

// ---------- 產生器 ----------
// K 個彼此相關 rho 的真觀點(壓格 + zscore)
export function makeIdeals(rng, K, M, rho, distribution) {
  const L = cholCS(K, rho);
  const Z = Array.from({ length: K }, () => Array.from({ length: M }, () => rng.n()));
  const idealZ = [];
  for (let k = 0; k < K; k++) {
    const y = new Array(M);
    for (let m = 0; m < M; m++) { let v = 0; for (let j = 0; j <= k; j++) v += L[k][j] * Z[j][m]; y[m] = v; }
    idealZ.push(zscore(generateCompositeSort(y, distribution)));
  }
  return idealZ;
}
// N 個受試者:每因素 perFactorN 人;共享變異 c;可摻非定義者(純噪音,不屬任何因素)
export function makeParticipants(rng, idealZ, perFactorN, c, mode, distribution, nonLoaderFrac = 0) {
  const K = idealZ.length, M = idealZ[0].length;
  const sorts = [], trueType = [];
  for (let k = 0; k < K; k++) for (let p = 0; p < perFactorN; p++) {
    const noisy = idealZ[k].map(v => Math.sqrt(c) * v + Math.sqrt(1 - c) * rng.n());
    sorts.push(render(noisy, mode, distribution)); trueType.push(k);
  }
  const nNon = Math.round(nonLoaderFrac * sorts.length / (1 - nonLoaderFrac || 1));
  for (let i = 0; i < nNon; i++) { sorts.push(render(Array.from({ length: M }, () => rng.n()), mode, distribution)); trueType.push(-1); }
  return { sorts, trueType };
}

// ---------- 因素保留規則 ----------
export function eigsOf(sorts) { return pcaExtract(correlationMatrix(sorts), 1).allEigenvalues; }
export function kaiserK(eigs) { return eigs.filter(e => e > 1).length; }
// 平行分析門檻(依 N,M,mode;pct=0.95 標準、0.99 保守),只依 null 分佈故每 (N,M,mode) 算一次
export function paThreshold(rng, N, M, mode, distribution, B = 60, pct = 0.95) {
  const cols = Array.from({ length: N }, () => []);
  for (let b = 0; b < B; b++) {
    const rs = Array.from({ length: N }, () => render(Array.from({ length: M }, () => rng.n()), mode, distribution));
    const e = eigsOf(rs); for (let i = 0; i < N; i++) cols[i].push(e[i]);
  }
  return cols.map(col => { col.sort((a, b) => a - b); return col[Math.floor(pct * (B - 1))]; });
}
export function paK(eigs, thr) { let k = 0; for (let i = 0; i < eigs.length; i++) if (eigs[i] > thr[i]) k++; return k; }

// ---------- 配對:精確 max-min 指派 ----------
// phi[t][f] = 植入 t 與還原 f 的 |φ|。回傳「讓最差配對 |φ| 最大」的一對一指派下,
// 每個植入觀點分到的 |φ|(還原因子數不足時,落單的植入觀點計 0)。
// 以「二分候選門檻 + 二部圖增廣配對」求解,是 README 所述判準
// 「每個植入觀點都由一個*不同的*還原因子在門檻以上配對」的精確實作
// (貪婪配對在衝突時可能漏掉存在的合法指派)。
export function maxMinAssignment(phi) {
  const K = phi.length, F0 = phi[0] ? phi[0].length : 0;
  const F = F0 + K; // 末尾補 K 個 φ=0 的虛擬因子,保證完美配對存在
  const val = (t, f) => (f < F0 ? phi[t][f] : 0);
  const matching = thr => {
    const matchOf = new Array(F).fill(-1);
    const tryAug = (u, seen) => {
      for (let v = 0; v < F; v++) {
        if (val(u, v) >= thr && !seen[v]) {
          seen[v] = true;
          if (matchOf[v] === -1 || tryAug(matchOf[v], seen)) { matchOf[v] = u; return true; }
        }
      }
      return false;
    };
    let size = 0;
    for (let u = 0; u < K; u++) if (tryAug(u, new Array(F).fill(false))) size++;
    return { size, matchOf };
  };
  const cands = [...new Set([0, ...phi.flat()])].sort((a, b) => a - b);
  let lo = 0, hi = cands.length - 1, best = matching(0);
  while (lo < hi) { // 最大的仍可完美配對的門檻
    const mid = (lo + hi + 1) >> 1;
    const m = matching(cands[mid]);
    if (m.size === K) { lo = mid; best = m; } else hi = mid - 1;
  }
  const phis = new Array(K).fill(0);
  for (let v = 0; v < F; v++) if (best.matchOf[v] >= 0) phis[best.matchOf[v]] = val(best.matchOf[v], v);
  return phis;
}

// ---------- 還原(PCA/centroid + varimax + flag + factorScores)+ 配對 ----------
// 因子陣列只由「自動 flag 的定義排序」計算(engine.mjs autoFlag:|l|>1.96/√M 且
// 該因子平方負荷 > 其他因子平方負荷總和),與 qmethod/KADE 的 auto-flag 一致;
// scoring='all' 保留舊行為(全體排序加權)供方法學比較。
// 回傳:{ estK, phis(最適指派下每個植入觀點的 |φ|) }
export function recover(sorts, idealZ, nFactors, extraction = 'pca', scoring = 'flag') {
  const K = idealZ.length;
  const nf = Math.max(1, Math.min(nFactors, sorts.length - 1));
  const cm = correlationMatrix(sorts);
  const { loadings } = extraction === 'centroid' ? centroidExtract(cm, nf) : pcaExtract(cm, nf);
  const rot = nf >= 2 ? varimax(loadings).loadings : loadings;
  const M = sorts[0].length;
  const flags = scoring === 'flag' ? autoFlag(rot, M) : null;
  const rec = [];
  for (let f = 0; f < nf; f++) {
    const col = rot.map((r, p) => (flags && !flags[p][f] ? 0 : r[f])); // 未 flag → 權重 0
    rec.push(factorScores(sorts, col, 0));
  }
  const phi = idealZ.map(iz => rec.map(r => Math.abs(tuckerPhi(iz, r))));
  return { estK: nf, phis: maxMinAssignment(phi) };
}

// ---------- Wilson 95% 信賴區間 ----------
export function wilson(x, n, z = 1.96) {
  if (n === 0) return [0, 0];
  const p = x / n, z2 = z * z;
  const denom = 1 + z2 / n;
  const center = (p + z2 / (2 * n)) / denom;
  const half = (z * Math.sqrt(p * (1 - p) / n + z2 / (4 * n * n))) / denom;
  return [Math.max(0, center - half), Math.min(1, center + half)];
}
