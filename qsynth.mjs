// Q-SYNTH library: synthetic Q participants + factor recovery.
// Self-contained: correlation/PCA/varimax/factorScores are open implementations in
// engine.mjs, and centroid is in centroid.mjs. Both extraction paths are validated
// against a real KADE export (PCA bit-exact; centroid structure col-corr >=0.9995;
// see validation/).
import { correlationMatrix, pcaExtract, varimax } from './engine.mjs';
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

// ---------- 還原(PCA + varimax + factorScores)+ 配對 ----------
// 回傳:{ estK, phis(每個植入觀點的最佳 |φ|) }
export function recover(sorts, idealZ, nFactors, extraction = 'pca') {
  const K = idealZ.length;
  const nf = Math.max(1, Math.min(nFactors, sorts.length - 1));
  const cm = correlationMatrix(sorts);
  const { loadings } = extraction === 'centroid' ? centroidExtract(cm, nf) : pcaExtract(cm, nf);
  const rot = nf >= 2 ? varimax(loadings).loadings : loadings;
  const rec = [];
  for (let f = 0; f < nf; f++) rec.push(factorScores(sorts, rot.map(r => r[f]), 0));
  const used = new Set(), phis = [];
  for (let t = 0; t < K; t++) {
    let bestAbs = 0, bf = -1;
    for (let f = 0; f < rec.length; f++) { if (used.has(f)) continue; const ph = Math.abs(tuckerPhi(idealZ[t], rec[f])); if (ph > bestAbs) { bestAbs = ph; bf = f; } }
    if (bf >= 0) used.add(bf); phis.push(bestAbs);
  }
  return { estK: nf, phis };
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
