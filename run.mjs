// Q-SYNTH 掃描器(平行,worker_threads)。以 PER-CELL 決定性種子跑整個設計空間:
// 每個設計格的亂數種子由該格參數雜湊而來,PA 門檻的種子由 (N,M,mode) 雜湊而來。
// 因此結果與「用幾個 worker、格子跑的順序」完全無關——既可平行、又完全可重現
//(同種子→同輸出),中斷後重跑天生冪等。輸出 results.csv(含 meanKaiserK)+ 圖。
//
// 用法:  node run.mjs [R]                 R=每格重抽次數,預設 1000
//         WORKERS=8 node run.mjs 1000      自訂 worker 數
import { Worker, isMainThread, parentPort, workerData } from 'node:worker_threads';
import os from 'node:os';
import { fileURLToPath } from 'node:url';
import { writeFileSync, mkdirSync } from 'node:fs';
import {
  makeRng, buildDistribution, makeIdeals, makeParticipants,
  eigsOf, paThreshold, paK, kaiserK, recover, wilson,
} from './qsynth.mjs';

// ---- 共用設定(main 與 worker 都載入)----
const R = parseInt(process.argv[2] || '1000', 10);
const B = 60;
const Ks = [2, 3, 4, 5, 6, 7];
const Ms = [20, 30, 40, 60, 80];
const perFactors = [4, 6, 8, 12];
const cs = [0.4, 0.6, 0.8];
const rhos = [0.0, 0.2, 0.4];
const modes = ['forced', 'free'];
const nonLoaders = [0.0, 0.2];
const retentions = ['oracle', 'pa'];
const totalNof = (pf, nlf, K) => { const base = pf * K; return base + Math.round(nlf * base / (1 - nlf || 1)); };

const SEED_BASE = 20260711;
function hseed(str) { let h = (2166136261 ^ SEED_BASE) >>> 0; for (let i = 0; i < str.length; i++) { h ^= str.charCodeAt(i); h = Math.imul(h, 16777619); } return h >>> 0; }

function enumerateCells() {
  const cells = []; let i = 0;
  for (const M of Ms) for (const K of Ks) for (const pf of perFactors) for (const nlf of nonLoaders) {
    const N = totalNof(pf, nlf, K);
    for (const mode of modes) for (const c of cs) for (const rho of rhos) for (const ret of retentions)
      cells.push({ i: i++, M, K, pf, nlf, N, mode, c, rho, ret });
  }
  return cells;
}

// 處理一批格子(worker 用)。distribution 依 M 快取;PA 依 (N,M,mode) 決定性種子快取。
// R/B 由參數傳入(worker 執行緒看不到 process.argv,不可依賴模組層級的 R)。
function processCells(cells, R, B) {
  const distCache = new Map(), paCache = new Map();
  const getDist = M => { if (!distCache.has(M)) distCache.set(M, buildDistribution(M)); return distCache.get(M); };
  const getPA = (N, M, mode, dist) => {
    const key = `${N}|${M}|${mode}`;
    if (!paCache.has(key)) paCache.set(key, paThreshold(makeRng(hseed(`pa|${key}`)), N, M, mode, dist, B, 0.95));
    return paCache.get(key);
  };
  const out = [];
  for (const { i, M, K, pf, nlf, N, mode, c, rho, ret } of cells) {
    const dist = getDist(M);
    const thr = ret === 'pa' ? getPA(N, M, mode, dist) : null;
    const rng = makeRng(hseed(`cell|${M}|${K}|${pf}|${nlf}|${mode}|${c}|${rho}|${ret}`));
    let r85 = 0, r90 = 0, r95 = 0, sumK = 0, sumKaiser = 0;
    for (let r = 0; r < R; r++) {
      const idealZ = makeIdeals(rng, K, M, rho, dist);
      const { sorts } = makeParticipants(rng, idealZ, pf, c, mode, dist, nlf);
      const eigs = eigsOf(sorts);
      const nf = ret === 'oracle' ? K : paK(eigs, thr);
      sumKaiser += kaiserK(eigs);
      const { estK, phis } = recover(sorts, idealZ, nf);
      sumK += estK;
      const worst = Math.min(...phis);
      if (worst >= 0.85) r85++;
      if (worst >= 0.90) r90++;
      if (worst >= 0.95) r95++;
    }
    const [lo, hi] = wilson(r90, R);
    out.push({ i, M, N, pf, K, c, rho, mode, nlf, ret, R, rec85: r85 / R, rec90: r90 / R, rec95: r95 / R, lo, hi, meanEstK: sumK / R, meanKaiserK: sumKaiser / R });
  }
  return out;
}

// ================= WORKER =================
if (!isMainThread) {
  parentPort.postMessage(processCells(workerData.cells, workerData.R, workerData.B));

// ================= MAIN =================
} else {
  const cells = enumerateCells();
  const W = Math.max(1, Math.min(cells.length, parseInt(process.env.WORKERS || '0', 10) || Math.max(1, os.cpus().length - 2)));
  const chunkSize = Math.ceil(cells.length / W);
  const chunks = [];
  for (let s = 0; s < cells.length; s += chunkSize) chunks.push(cells.slice(s, s + chunkSize));
  const __filename = fileURLToPath(import.meta.url);
  const t0 = Date.now();
  console.log(`Q-SYNTH 平行掃描:${cells.length} 格 × R=${R},${chunks.length} workers…`);

  const results = await Promise.all(chunks.map(chunk => new Promise((resolve, reject) => {
    const w = new Worker(__filename, { workerData: { cells: chunk, R, B } });
    let acc = null;
    w.on('message', m => { acc = m; });
    w.on('error', reject);
    w.on('exit', code => code === 0 ? resolve(acc) : reject(new Error('worker exit ' + code)));
  })));
  const rows = results.flat().sort((a, b) => a.i - b.i);
  const secs = ((Date.now() - t0) / 1000).toFixed(1);

  // ---- CSV ----
  const header = ['M', 'N', 'perFactorN', 'K', 'c', 'rho', 'mode', 'nonLoaderFrac', 'retention', 'R',
    'rec85', 'rec90', 'rec95', 'rec90_lo', 'rec90_hi', 'meanEstK', 'meanKaiserK'];
  const csv = [header.join(',')].concat(rows.map(r =>
    [r.M, r.N, r.pf, r.K, r.c, r.rho, r.mode, r.nlf, r.ret, r.R,
      r.rec85.toFixed(3), r.rec90.toFixed(3), r.rec95.toFixed(3), r.lo.toFixed(3), r.hi.toFixed(3), r.meanEstK.toFixed(2), r.meanKaiserK.toFixed(2)].join(',')));
  writeFileSync(new URL('./results.csv', import.meta.url), csv.join('\n') + '\n');

  // ---- 圖(K=3 切片,與序列版同)----
  const cell = (M, pf, c, rho, mode, nlf, ret, K = 3) => rows.find(r =>
    r.M === M && r.pf === pf && r.c === c && r.rho === rho && r.mode === mode && r.nlf === nlf && r.ret === ret && r.K === K);
  const colors = ['#1b7', '#28c', '#e63'];
  function lineChart(title, xs, series, yLabel = '還原率 @|φ|≥0.90') {
    const W2 = 520, H = 380, ml = 60, mr = 130, mt = 44, mb = 50, pw = W2 - ml - mr, ph = H - mt - mb;
    const xAt = i => ml + (xs.length === 1 ? pw / 2 : (i / (xs.length - 1)) * pw), yAt = v => mt + (1 - v) * ph;
    const e = [`<svg xmlns="http://www.w3.org/2000/svg" width="${W2}" height="${H}" font-family="sans-serif" font-size="12">`, `<rect width="${W2}" height="${H}" fill="white"/>`];
    e.push(`<text x="${W2 / 2}" y="24" text-anchor="middle" font-size="15" font-weight="bold">${title}</text>`);
    for (let g = 0; g <= 10; g += 2) { const y = yAt(g / 10); e.push(`<line x1="${ml}" y1="${y}" x2="${ml + pw}" y2="${y}" stroke="#eee"/><text x="${ml - 8}" y="${y + 4}" text-anchor="end" fill="#555">${g * 10}%</text>`); }
    e.push(`<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt + ph}" stroke="#333"/><line x1="${ml}" y1="${mt + ph}" x2="${ml + pw}" y2="${mt + ph}" stroke="#333"/>`);
    xs.forEach((x, i) => e.push(`<text x="${xAt(i)}" y="${mt + ph + 20}" text-anchor="middle" fill="#555">${x}</text>`));
    e.push(`<text x="${ml + pw / 2}" y="${H - 8}" text-anchor="middle" fill="#333">每因素人數(persons/factor)</text>`);
    e.push(`<text x="16" y="${mt + ph / 2}" text-anchor="middle" transform="rotate(-90 16 ${mt + ph / 2})" fill="#333">${yLabel}</text>`);
    series.forEach((s, si) => {
      e.push(`<polyline points="${s.ys.map((v, i) => `${xAt(i)},${yAt(v)}`).join(' ')}" fill="none" stroke="${s.color}" stroke-width="2.5"/>`);
      s.ys.forEach((v, i) => e.push(`<circle cx="${xAt(i)}" cy="${yAt(v)}" r="3.5" fill="${s.color}"/>`));
      const ly = mt + 6 + si * 20;
      e.push(`<line x1="${ml + pw + 14}" y1="${ly}" x2="${ml + pw + 34}" y2="${ly}" stroke="${s.color}" stroke-width="2.5"/><text x="${ml + pw + 38}" y="${ly + 4}" fill="#333">${s.label}</text>`);
    });
    e.push('</svg>'); return e.join('\n');
  }
  function figFor(rho, mode, nlf, fname, subtitle) {
    const series = cs.map((c, i) => ({ label: `c=${c}`, color: colors[i], ys: perFactors.map(pf => cell(40, pf, c, rho, mode, nlf, 'pa').rec90) }));
    writeFileSync(new URL(`./figures/${fname}`, import.meta.url), lineChart(`Q 因素還原率(M=40,ρ=${rho},${subtitle},平行分析)`, perFactors, series));
  }
  mkdirSync(new URL('./figures/', import.meta.url), { recursive: true });
  figFor(0.0, 'forced', 0.0, 'fig_rho0_forced.svg', '強制、無非定義者');
  figFor(0.4, 'forced', 0.0, 'fig_rho04_forced.svg', '強制、無非定義者');
  figFor(0.0, 'forced', 0.2, 'fig_rho0_nonload20.svg', '強制、20% 非定義者');

  // ---- 摘要 ----
  console.log(`完成:${rows.length} 格 × R=${R},${chunks.length} workers,耗時 ${secs}s → results.csv + figures/`);
  console.log('摘要(M=40,c=0.6,強制,無非定義者,平行分析,還原@|φ|≥0.90 [95% CI])');
  console.log('每因素人數 |   ρ=0.0            ρ=0.4');
  for (const pf of perFactors) {
    const f = rho => { const p = cell(40, pf, 0.6, rho, 'forced', 0.0, 'pa'); return `${Math.round(p.rec90 * 100)}% [${Math.round(p.lo * 100)}–${Math.round(p.hi * 100)}]`.padStart(15); };
    console.log(`   ${String(pf).padStart(2)}      |${f(0.0)}   ${f(0.4)}`);
  }
}
