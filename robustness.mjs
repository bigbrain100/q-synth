// Q-SYNTH robustness sweep — is the M×N recovery standard robust to the choice
// of factor-extraction method? Compares PCA (the default path,
// bit-exact to KADE) against centroid (traditional Q method, KADE-validated
// structure ≥0.9995 col-corr; see validation/) on the SAME synthetic datasets
// (paired design), so any difference is the extraction method alone.
//
// Retention is fixed to oracle (nf = K planted factors) to isolate the
// extraction effect from factor-retention rules. mode=forced, no non-loaders.
//
// Usage:  node robustness.mjs [R]     (R = replicates per cell, default 400)
import { writeFileSync, mkdirSync } from 'node:fs';
import {
  makeRng, buildDistribution, makeIdeals, makeParticipants, recover, wilson,
} from './qsynth.mjs';

const R = parseInt(process.argv[2] || '400', 10);
const K = 3;
const Ms = [20, 40];
const perFactors = [4, 6, 8, 12];
const cs = [0.4, 0.6, 0.8];
const rhos = [0.0, 0.4];
const extractions = ['pca', 'centroid'];

const rng = makeRng(20260712);
const header = ['M', 'N', 'perFactorN', 'K', 'c', 'rho', 'extraction', 'R', 'rec90', 'rec90_lo', 'rec90_hi', 'meanMinPhi'];
const rows = [];

const t0 = Date.now();
for (const M of Ms) {
  const distribution = buildDistribution(M);
  for (const pf of perFactors) {
    const N = pf * K;
    for (const c of cs) for (const rho of rhos) {
      // paired: draw R datasets once, score both extraction methods on each
      const acc = { pca: { r90: 0, sum: 0 }, centroid: { r90: 0, sum: 0 } };
      for (let r = 0; r < R; r++) {
        const idealZ = makeIdeals(rng, K, M, rho, distribution);
        const { sorts } = makeParticipants(rng, idealZ, pf, c, 'forced', distribution, 0);
        for (const ext of extractions) {
          const { phis } = recover(sorts, idealZ, K, ext);
          const worst = Math.min(...phis);
          acc[ext].sum += worst;
          if (worst >= 0.90) acc[ext].r90++;
        }
      }
      for (const ext of extractions) {
        const [lo, hi] = wilson(acc[ext].r90, R);
        rows.push({ M, N, pf, c, rho, ext, R, rec90: acc[ext].r90 / R, lo, hi, meanMinPhi: acc[ext].sum / R });
      }
    }
  }
}
const secs = ((Date.now() - t0) / 1000).toFixed(1);

// ---- CSV ----
const csv = [header.join(',')].concat(rows.map(r =>
  [r.M, r.N, r.pf, K, r.c, r.rho, r.ext, r.R,
    r.rec90.toFixed(3), r.lo.toFixed(3), r.hi.toFixed(3), r.meanMinPhi.toFixed(3)].join(',')));
writeFileSync(new URL('./robustness.csv', import.meta.url), csv.join('\n') + '\n');

const cell = (M, pf, c, rho, ext) => rows.find(r => r.M === M && r.pf === pf && r.c === c && r.rho === rho && r.ext === ext);

// ---- paired figure: PCA vs centroid recovery curves (M=40, ρ=0.4) ----
function pairFig(M, rho, fname, subtitle) {
  const W = 560, H = 400, ml = 62, mr = 150, mt = 46, mb = 54, pw = W - ml - mr, ph = H - mt - mb;
  const xAt = i => ml + (i / (perFactors.length - 1)) * pw, yAt = v => mt + (1 - v) * ph;
  const e = [`<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" font-family="sans-serif" font-size="12">`];
  e.push(`<rect width="${W}" height="${H}" fill="white"/>`);
  e.push(`<text x="${W / 2}" y="24" text-anchor="middle" font-size="15" font-weight="bold">PCA vs centroid 還原率(M=${M},${subtitle})</text>`);
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
console.log(`Q-SYNTH 穩健性掃描(PCA vs centroid,配對):${rows.length} 列 × R=${R},耗時 ${secs}s`);
console.log('→ robustness.csv;figures/fig_robust_pca_vs_centroid.svg\n');
let maxGap = 0, maxGapAt = '';
for (const M of Ms) for (const rho of rhos) {
  console.log(`M=${M}, ρ=${rho}(強制,oracle,還原@|φ|≥0.90)`);
  console.log('  c \\ perF |' + perFactors.map(p => String(p).padStart(13)).join(''));
  for (const c of cs) {
    const cells = perFactors.map(pf => {
      const p = cell(M, pf, c, rho, 'pca').rec90, q = cell(M, pf, c, rho, 'centroid').rec90;
      const gap = q - p; if (Math.abs(gap) > Math.abs(maxGap)) { maxGap = gap; maxGapAt = `M=${M} ρ=${rho} c=${c} pf=${pf}`; }
      return `${Math.round(p * 100)}/${Math.round(q * 100)}`.padStart(13);
    });
    console.log(`   ${c}     |` + cells.join(''));
  }
  console.log('');
}
console.log('每格顯示 PCA%/centroid%。最大差距(centroid−PCA)= ' + (maxGap >= 0 ? '+' : '') + `${Math.round(maxGap * 100)}pp @ ${maxGapAt}`);
