// Q-SYNTH 整合掃描器(投稿級)—— 掃 M×每因素人數×c×ρ×保留規則×分佈型態×非定義者比例,
// 附 Wilson CI,輸出 results.csv 與 SVG 折線圖。
// 用法:node run.mjs [R]     (R = 每格重抽次數,預設 500)
import { writeFileSync, mkdirSync } from 'node:fs';
import {
  makeRng, buildDistribution, makeIdeals, makeParticipants,
  eigsOf, paThreshold, paK, recover, wilson,
} from './qsynth.mjs';

const R = parseInt(process.argv[2] || '500', 10);
const K = 3, B = 60;
const Ms = [20, 40];
const perFactors = [4, 6, 8, 12];
const cs = [0.4, 0.6, 0.8];
const rhos = [0.0, 0.4];
const modes = ['forced', 'free'];
const nonLoaders = [0.0, 0.2];
const retentions = ['oracle', 'pa'];

const rng = makeRng(20260711);
const totalNof = (pf, nlf) => { const base = pf * K; return base + Math.round(nlf * base / (1 - nlf || 1)); };
const paCache = {};
function getPA(N, M, mode, distribution) {
  const key = `${N}|${M}|${mode}`;
  if (!paCache[key]) paCache[key] = paThreshold(rng, N, M, mode, distribution, B, 0.95);
  return paCache[key];
}

const header = ['M', 'N', 'perFactorN', 'K', 'c', 'rho', 'mode', 'nonLoaderFrac', 'retention', 'R',
  'rec85', 'rec90', 'rec95', 'rec90_lo', 'rec90_hi', 'meanEstK'];
const rows = [];

const t0 = Date.now();
for (const M of Ms) {
  const distribution = buildDistribution(M);
  for (const pf of perFactors) for (const nlf of nonLoaders) {
    const N = totalNof(pf, nlf);
    for (const mode of modes) {
      const thr = getPA(N, M, mode, distribution);
      for (const c of cs) for (const rho of rhos) for (const ret of retentions) {
        let r85 = 0, r90 = 0, r95 = 0, sumK = 0;
        for (let r = 0; r < R; r++) {
          const idealZ = makeIdeals(rng, K, M, rho, distribution);
          const { sorts } = makeParticipants(rng, idealZ, pf, c, mode, distribution, nlf);
          const nf = ret === 'oracle' ? K : paK(eigsOf(sorts), thr);
          const { estK, phis } = recover(sorts, idealZ, nf);
          sumK += estK;
          const worst = Math.min(...phis);
          if (worst >= 0.85) r85++;
          if (worst >= 0.90) r90++;
          if (worst >= 0.95) r95++;
        }
        const [lo, hi] = wilson(r90, R);
        rows.push({
          M, N, pf, K, c, rho, mode, nlf, ret, R,
          rec85: r85 / R, rec90: r90 / R, rec95: r95 / R, lo, hi, meanEstK: sumK / R,
        });
      }
    }
  }
}
const secs = ((Date.now() - t0) / 1000).toFixed(1);

// ---- CSV ----
const csv = [header.join(',')].concat(rows.map(r =>
  [r.M, r.N, r.pf, r.K, r.c, r.rho, r.mode, r.nlf, r.ret, r.R,
    r.rec85.toFixed(3), r.rec90.toFixed(3), r.rec95.toFixed(3), r.lo.toFixed(3), r.hi.toFixed(3), r.meanEstK.toFixed(2)].join(',')));
writeFileSync(new URL('./results.csv', import.meta.url), csv.join('\n') + '\n');

// ---- SVG 折線圖 ----
function lineChart(title, xs, series, yLabel = '還原率 @|φ|≥0.90') {
  const W = 520, H = 380, ml = 60, mr = 130, mt = 44, mb = 50;
  const pw = W - ml - mr, ph = H - mt - mb;
  const xAt = i => ml + (xs.length === 1 ? pw / 2 : (i / (xs.length - 1)) * pw);
  const yAt = v => mt + (1 - v) * ph;
  const e = [`<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" font-family="sans-serif" font-size="12">`];
  e.push(`<rect width="${W}" height="${H}" fill="white"/>`);
  e.push(`<text x="${W / 2}" y="24" text-anchor="middle" font-size="15" font-weight="bold">${title}</text>`);
  for (let g = 0; g <= 10; g += 2) { const y = yAt(g / 10); e.push(`<line x1="${ml}" y1="${y}" x2="${ml + pw}" y2="${y}" stroke="#eee"/>`); e.push(`<text x="${ml - 8}" y="${y + 4}" text-anchor="end" fill="#555">${g * 10}%</text>`); }
  e.push(`<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt + ph}" stroke="#333"/>`);
  e.push(`<line x1="${ml}" y1="${mt + ph}" x2="${ml + pw}" y2="${mt + ph}" stroke="#333"/>`);
  xs.forEach((x, i) => e.push(`<text x="${xAt(i)}" y="${mt + ph + 20}" text-anchor="middle" fill="#555">${x}</text>`));
  e.push(`<text x="${ml + pw / 2}" y="${H - 8}" text-anchor="middle" fill="#333">每因素人數(persons/factor)</text>`);
  e.push(`<text x="16" y="${mt + ph / 2}" text-anchor="middle" transform="rotate(-90 16 ${mt + ph / 2})" fill="#333">${yLabel}</text>`);
  series.forEach((s, si) => {
    const pts = s.ys.map((v, i) => `${xAt(i)},${yAt(v)}`).join(' ');
    e.push(`<polyline points="${pts}" fill="none" stroke="${s.color}" stroke-width="2.5"/>`);
    s.ys.forEach((v, i) => e.push(`<circle cx="${xAt(i)}" cy="${yAt(v)}" r="3.5" fill="${s.color}"/>`));
    const ly = mt + 6 + si * 20;
    e.push(`<line x1="${ml + pw + 14}" y1="${ly}" x2="${ml + pw + 34}" y2="${ly}" stroke="${s.color}" stroke-width="2.5"/>`);
    e.push(`<text x="${ml + pw + 38}" y="${ly + 4}" fill="#333">${s.label}</text>`);
  });
  // Brown 4–5 參考帶
  const bx1 = xAt(0), bx2 = xAt(1);
  e.push(`<rect x="${bx1}" y="${mt}" width="${bx2 - bx1}" height="${ph}" fill="#f3a" opacity="0.06"/>`);
  e.push(`<text x="${(bx1 + bx2) / 2}" y="${mt + ph - 6}" text-anchor="middle" fill="#c39" font-size="11">Brown 4–5</text>`);
  e.push('</svg>');
  return e.join('\n');
}

const colors = ['#1b7', '#28c', '#e63'];
const cell = (M, pf, c, rho, mode, nlf, ret) => rows.find(r =>
  r.M === M && r.pf === pf && r.c === c && r.rho === rho && r.mode === mode && r.nlf === nlf && r.ret === ret);
function figFor(rho, mode, nlf, fname, subtitle) {
  const series = cs.map((c, i) => ({
    label: `c=${c}`, color: colors[i],
    ys: perFactors.map(pf => cell(40, pf, c, rho, mode, nlf, 'pa').rec90),
  }));
  const svg = lineChart(`Q 因素還原率(M=40,ρ=${rho},${subtitle},平行分析)`, perFactors, series);
  writeFileSync(new URL(`./figures/${fname}`, import.meta.url), svg);
}
mkdirSync(new URL('./figures/', import.meta.url), { recursive: true });
figFor(0.0, 'forced', 0.0, 'fig_rho0_forced.svg', '強制、無非定義者');
figFor(0.4, 'forced', 0.0, 'fig_rho04_forced.svg', '強制、無非定義者');
figFor(0.0, 'forced', 0.2, 'fig_rho0_nonload20.svg', '強制、20% 非定義者');

// ---- 主控台摘要 ----
console.log(`Q-SYNTH 投稿級掃描完成:${rows.length} 格 × R=${R},耗時 ${secs}s`);
console.log('→ results.csv;figures/{fig_rho0_forced,fig_rho04_forced,fig_rho0_nonload20}.svg\n');
console.log('摘要(M=40,c=0.6,強制,無非定義者,平行分析,還原@|φ|≥0.90 [95% CI])');
console.log('每因素人數 |   ρ=0.0            ρ=0.4');
console.log('-----------|----------------------------------');
for (const pf of perFactors) {
  const f = rho => { const p = cell(40, pf, 0.6, rho, 'forced', 0.0, 'pa'); return `${Math.round(p.rec90 * 100)}% [${Math.round(p.lo * 100)}–${Math.round(p.hi * 100)}]`.padStart(15); };
  console.log(`   ${String(pf).padStart(2)}      |${f(0.0)}   ${f(0.4)}`);
}
console.log('\n非定義者(20%)對還原的影響、自由 vs 強制、oracle vs pa、c=0.4/0.8、φ 三門檻 → results.csv');