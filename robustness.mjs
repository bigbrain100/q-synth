// Q-SYNTH robustness sweep (parallel) — is the recovery standard robust to the
// factor-extraction method? Compares PCA (default, bit-exact to KADE) against
// centroid (traditional Q, KADE-validated structure >=0.9995) on the SAME
// synthetic datasets (paired design), so any difference is the extraction method
// alone. Retention fixed to oracle (nf = K planted), mode=forced, no non-loaders.
//
// PER-CELL deterministic seeding (matches run.mjs): each (M,pf,c,rho) cell seeds
// its own generator, so results are independent of worker count and order, and
// the parallel run is bit-identical to a serial one.
//
// Usage:  node robustness.mjs [R]            R = replicates per cell, default 1000
//         WORKERS=8 node robustness.mjs 1000  custom worker count
import { Worker, isMainThread, parentPort, workerData } from 'node:worker_threads';
import os from 'node:os';
import { fileURLToPath } from 'node:url';
import { writeFileSync, mkdirSync } from 'node:fs';
import { makeRng, buildDistribution, makeIdeals, makeParticipants, recover, wilson } from './qsynth.mjs';

const R = parseInt(process.argv[2] || '1000', 10);
const K = 3;
const Ms = [20, 40];
const perFactors = [4, 6, 8, 12];
const cs = [0.4, 0.6, 0.8];
const rhos = [0.0, 0.2, 0.4];
const extractions = ['pca', 'centroid'];

const SEED_BASE = 20260712;
function hseed(str) { let h = (2166136261 ^ SEED_BASE) >>> 0; for (let i = 0; i < str.length; i++) { h ^= str.charCodeAt(i); h = Math.imul(h, 16777619); } return h >>> 0; }

function enumerateCells() {
  const cells = []; let i = 0;
  for (const M of Ms) for (const pf of perFactors) for (const c of cs) for (const rho of rhos)
    cells.push({ i: i++, M, pf, N: pf * K, c, rho });
  return cells;
}

function processCells(cells, R) {
  const distCache = new Map();
  const getDist = M => { if (!distCache.has(M)) distCache.set(M, buildDistribution(M)); return distCache.get(M); };
  const out = [];
  for (const { i, M, pf, N, c, rho } of cells) {
    const dist = getDist(M);
    const rng = makeRng(hseed(`rob|${M}|${pf}|${c}|${rho}`));
    const acc = { pca: { r90: 0, sum: 0 }, centroid: { r90: 0, sum: 0 } };
    for (let r = 0; r < R; r++) {
      const idealZ = makeIdeals(rng, K, M, rho, dist);
      const { sorts } = makeParticipants(rng, idealZ, pf, c, 'forced', dist, 0);
      for (const ext of extractions) {
        const { phis } = recover(sorts, idealZ, K, ext);
        const worst = Math.min(...phis);
        acc[ext].sum += worst;
        if (worst >= 0.90) acc[ext].r90++;
      }
    }
    for (const ext of extractions) {
      const [lo, hi] = wilson(acc[ext].r90, R);
      out.push({ i, M, N, pf, c, rho, ext, R, rec90: acc[ext].r90 / R, lo, hi, meanMinPhi: acc[ext].sum / R });
    }
  }
  return out;
}

// ================= WORKER =================
if (!isMainThread) {
  parentPort.postMessage(processCells(workerData.cells, workerData.R));

// ================= MAIN =================
} else {
  const cells = enumerateCells();
  const W = Math.max(1, Math.min(cells.length, parseInt(process.env.WORKERS || '0', 10) || Math.max(1, os.cpus().length - 2)));
  const chunkSize = Math.ceil(cells.length / W);
  const chunks = [];
  for (let s = 0; s < cells.length; s += chunkSize) chunks.push(cells.slice(s, s + chunkSize));
  const __filename = fileURLToPath(import.meta.url);
  const t0 = Date.now();
  console.log(`Q-SYNTH 穩健性平行掃描:${cells.length} 格 × R=${R} × 2 萃取,${chunks.length} workers…`);

  const results = await Promise.all(chunks.map(chunk => new Promise((resolve, reject) => {
    const w = new Worker(__filename, { workerData: { cells: chunk, R } });
    let acc = null;
    w.on('message', m => { acc = m; });
    w.on('error', reject);
    w.on('exit', code => code === 0 ? resolve(acc) : reject(new Error('worker exit ' + code)));
  })));
  const rows = results.flat().sort((a, b) => a.i - b.i || (a.ext < b.ext ? -1 : 1));
  const secs = ((Date.now() - t0) / 1000).toFixed(1);

  // ---- CSV ----
  const header = ['M', 'N', 'perFactorN', 'K', 'c', 'rho', 'extraction', 'R', 'rec90', 'rec90_lo', 'rec90_hi', 'meanMinPhi'];
  const csv = [header.join(',')].concat(rows.map(r =>
    [r.M, r.N, r.pf, K, r.c, r.rho, r.ext, r.R, r.rec90.toFixed(3), r.lo.toFixed(3), r.hi.toFixed(3), r.meanMinPhi.toFixed(3)].join(',')));
  writeFileSync(new URL('./robustness.csv', import.meta.url), csv.join('\n') + '\n');

  const cell = (M, pf, c, rho, ext) => rows.find(r => r.M === M && r.pf === pf && r.c === c && r.rho === rho && r.ext === ext);

  // ---- paired figure: PCA vs centroid recovery curves (M=40, rho=0.4) ----
  function pairFig(M, rho, fname, subtitle) {
    const Wd = 560, H = 400, ml = 62, mr = 150, mt = 46, mb = 54, pw = Wd - ml - mr, ph = H - mt - mb;
    const xAt = i => ml + (i / (perFactors.length - 1)) * pw, yAt = v => mt + (1 - v) * ph;
    const e = [`<svg xmlns="http://www.w3.org/2000/svg" width="${Wd}" height="${H}" font-family="sans-serif" font-size="12">`];
    e.push(`<rect width="${Wd}" height="${H}" fill="white"/>`);
    e.push(`<text x="${Wd / 2}" y="24" text-anchor="middle" font-size="15" font-weight="bold">PCA vs centroid 還原率(M=${M},${subtitle})</text>`);
    for (let g = 0; g <= 10; g += 2) { const y = yAt(g / 10); e.push(`<line x1="${ml}" y1="${y}" x2="${ml + pw}" y2="${y}" stroke="#eee"/>`); e.push(`<text x="${ml - 8}" y="${y + 4}" text-anchor="end" fill="#555">${g * 10}%</text>`); }
    e.push(`<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt + ph}" stroke="#333"/>`);
    e.push(`<line x1="${ml}" y1="${mt + ph}" x2="${ml + pw}" y2="${mt + ph}" stroke="#333"/>`);
    perFactors.forEach((x, i) => e.push(`<text x="${xAt(i)}" y="${mt + ph + 20}" text-anchor="middle" fill="#555">${x}</text>`));
    e.push(`<text x="${ml + pw / 2}" y="${H - 8}" text-anchor="middle" fill="#333">每因素人數(persons/factor)</text>`);
    e.push(`<text x="16" y="${mt + ph / 2}" text-anchor="middle" transform="rotate(-90 16 ${mt + ph / 2})" fill="#333">還原率 @|φ|≥0.90</text>`);
    const styles = { pca: { c: '#28c', dash: '' }, centroid: { c: '#e63', dash: '6 4' } };
    const cList = [0.4, 0.6, 0.8];
    let si = 0;
    for (const cc of cList) for (const ext of extractions) {
      const ys = perFactors.map(pf => cell(M, pf, cc, rho, ext).rec90);
      const st = styles[ext];
      const pts = ys.map((v, i) => `${xAt(i)},${yAt(v)}`).join(' ');
      e.push(`<polyline points="${pts}" fill="none" stroke="${st.c}" stroke-width="2.2" stroke-dasharray="${st.dash}" opacity="${0.4 + 0.3 * cList.indexOf(cc)}"/>`);
      ys.forEach((v, i) => e.push(`<circle cx="${xAt(i)}" cy="${yAt(v)}" r="3" fill="${st.c}" opacity="${0.4 + 0.3 * cList.indexOf(cc)}"/>`));
      const ly = mt + 6 + si * 18;
      e.push(`<line x1="${ml + pw + 12}" y1="${ly}" x2="${ml + pw + 30}" y2="${ly}" stroke="${st.c}" stroke-width="2.2" stroke-dasharray="${st.dash}" opacity="${0.4 + 0.3 * cList.indexOf(cc)}"/>`);
      e.push(`<text x="${ml + pw + 34}" y="${ly + 4}" fill="#333">c=${cc} ${ext}</text>`);
      si++;
    }
    e.push('</svg>');
    writeFileSync(new URL(`./figures/${fname}`, import.meta.url), e.join('\n'));
  }
  mkdirSync(new URL('./figures/', import.meta.url), { recursive: true });
  pairFig(40, 0.4, 'fig_robust_pca_vs_centroid.svg', 'ρ=0.4,強制,oracle');

  // ---- console summary ----
  console.log(`完成:${rows.length} 列 × R=${R},${chunks.length} workers,耗時 ${secs}s → robustness.csv`);
  let maxGap = 0, maxGapAt = '';
  for (const M of Ms) for (const rho of rhos) {
    console.log(`M=${M}, ρ=${rho}(強制,oracle,還原@|φ|≥0.90)`);
    console.log('  c \\ perF |' + perFactors.map(p => String(p).padStart(13)).join(''));
    for (const c of cs) {
      const cells2 = perFactors.map(pf => {
        const p = cell(M, pf, c, rho, 'pca').rec90, q = cell(M, pf, c, rho, 'centroid').rec90;
        const gap = q - p; if (Math.abs(gap) > Math.abs(maxGap)) { maxGap = gap; maxGapAt = `M=${M} ρ=${rho} c=${c} pf=${pf}`; }
        return `${Math.round(p * 100)}/${Math.round(q * 100)}`.padStart(13);
      });
      console.log(`   ${c}     |` + cells2.join(''));
    }
    console.log('');
  }
  console.log('每格 PCA%/centroid%。最大差距(centroid−PCA)= ' + (maxGap >= 0 ? '+' : '') + `${Math.round(maxGap * 100)}pp @ ${maxGapAt}`);
}
