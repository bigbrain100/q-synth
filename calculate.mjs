// Q-SYNTH design calculator.
//
// Given the kind of viewpoints a planned Q study expects, recommend how many
// defining sorts per factor (and how many statements) are needed for the factors
// to be recovered reliably. It reads the simulation lookup table (results.csv)
// and returns the smallest design that reaches a target recovery rate.
//
// Usage:
//   node calculate.mjs --c 0.6 --rho 0.2 --M 40 --K 3 --target 0.90
//   node calculate.mjs --c 0.6 --rho 0.4 --threshold 0.95 --mode forced --retention pa
//
// Options (all optional; sensible defaults shown):
//   --c          expected communality of defining sorts (0-1), default 0.6
//   --rho        expected correlation between viewpoints (0-1), default 0.2.
//                Higher rho = viewpoints OVERLAP more = harder to separate = needs more sorts.
//                (Low rho means distinct viewpoints.)
//   --M          number of statements you plan to use, default 40
//   --K          number of expected viewpoints (factors), default 3
//   --target     target recovery rate (0-1), default 0.90
//   --threshold  congruence threshold to read (0.85|0.90|0.95), default 0.90
//   --mode       forced | free, default forced
//   --retention  pa (parallel analysis) | oracle, default pa
//   --nonloaders 0 | 0.2 expected share of non-defining participants, default 0
//
// The simulated grid is coarse (c in {0.4,0.6,0.8}; rho in {0,0.2,0.4}; M in {20,30,40,60,80};
// persons/factor in {4,6,8,12}; K in {2,3,4,5,6,7}, seven being the magic-number-seven
// ceiling of factors Brown would preliminarily extract). Inputs are snapped to the
// nearest grid value, leaning to the SAFER side (lower c, higher rho) so a
// recommendation is never optimistic. This is planning guidance to read with your
// topic in mind, not a recipe.

import { readFileSync } from 'node:fs';

const GRID = { c: [0.4, 0.6, 0.8], rho: [0.0, 0.2, 0.4], M: [20, 30, 40, 60, 80], perFactorN: [4, 6, 8, 12], K: [2, 3, 4, 5, 6, 7] };

function parseArgs(argv) {
  const o = {};
  for (let i = 0; i < argv.length; i++) if (argv[i].startsWith('--')) o[argv[i].slice(2)] = argv[i + 1];
  return o;
}
// snap to nearest grid value; `safe` picks the lower (for c) / higher (for rho) neighbour on ties-toward-caution
function snap(val, options, safer) {
  const sorted = [...options].sort((a, b) => a - b);
  let best = sorted[0], bd = Infinity;
  for (const g of sorted) { const d = Math.abs(g - val); if (d < bd - 1e-9) { bd = d; best = g; } }
  // conservative nudge: if strictly between two grid points, lean to the harder condition
  if (val > sorted[0] && val < sorted[sorted.length - 1]) {
    const below = [...sorted].reverse().find(g => g <= val);
    const above = sorted.find(g => g >= val);
    if (below !== above) best = safer === 'low' ? below : above;
  }
  return best;
}

function loadTable(url) {
  const text = readFileSync(url, 'utf8').trim().split('\n');
  const head = text[0].split(',');
  return text.slice(1).map(line => {
    const f = line.split(',');
    const row = {};
    head.forEach((h, i) => { row[h] = f[i]; });
    return row;
  });
}

function main() {
  const a = parseArgs(process.argv.slice(2));
  const cIn = a.c !== undefined ? +a.c : 0.6;
  const rhoIn = a.rho !== undefined ? +a.rho : 0.2;
  const MIn = a.M !== undefined ? +a.M : 40;
  const Kin = a.K !== undefined ? +a.K : 3;
  const target = a.target !== undefined ? +a.target : 0.90;
  const thr = a.threshold !== undefined ? a.threshold : '0.90';
  const mode = a.mode || 'forced';
  const retention = a.retention || 'pa';
  const nlf = a.nonloaders !== undefined ? a.nonloaders : '0';
  const recCol = ({ '0.85': 'rec85', '0.90': 'rec90', '0.95': 'rec95' })[thr] || 'rec90';

  // ---- input validation ----
  const errs = [];
  for (const [name, v, raw] of [['c', cIn, a.c], ['rho', rhoIn, a.rho], ['target', target, a.target]]) {
    if (!Number.isFinite(v) || v < 0 || v > 1) errs.push(`--${name} must be a number in [0, 1] (got ${raw})`);
  }
  if (!Number.isFinite(MIn) || MIn < 2) errs.push(`--M must be a positive integer (got ${a.M})`);
  if (!['0.85', '0.90', '0.95'].includes(thr)) errs.push(`--threshold must be one of 0.85, 0.90, 0.95 (got ${thr})`);
  if (!['forced', 'free'].includes(mode)) errs.push(`--mode must be forced or free (got ${mode})`);
  if (!['pa', 'oracle'].includes(retention)) errs.push(`--retention must be pa or oracle (got ${retention})`);
  if (!['0', '0.2'].includes(String(nlf))) errs.push(`--nonloaders must be 0 or 0.2 (got ${nlf})`);
  if (errs.length) { console.error('Invalid input:\n  ' + errs.join('\n  ')); process.exit(1); }
  if (!Number.isFinite(Kin) || Kin < 1) { console.error(`Invalid input:\n  --K must be a positive integer (got ${a.K})`); process.exit(1); }
  const K = snap(Math.round(Kin), GRID.K, 'high'); // expected viewpoint count, snapped to the simulated grid (2..7)
  const kCeiling = Math.round(Kin) > 7;            // beyond the magic-number-seven ceiling

  const M = snap(MIn, GRID.M, 'low');       // safer = fewer statements
  const cCons = snap(cIn, GRID.c, 'low');    // conservative = lower communality
  const rhoCons = snap(rhoIn, GRID.rho, 'high'); // conservative = higher (more overlapping) distinctness
  const cOpt = snap(cIn, GRID.c, 'high');
  const rhoOpt = snap(rhoIn, GRID.rho, 'low');
  const between = cCons !== cOpt || rhoCons !== rhoOpt;

  const table = loadTable(new URL('./results.csv', import.meta.url));
  // returns {pf, rec, lo, hi} for the smallest perFactorN hitting target, or {fail, best}
  function recommend(c, rho) {
    const rows = table.filter(r => +r.M === M && +r.c === c && +r.rho === rho
      && r.mode === mode && r.retention === retention && +r.K === K && r.nonLoaderFrac === String(nlf))
      .sort((x, y) => +x.perFactorN - +y.perFactorN);
    if (!rows.length) return null;
    const hit = rows.find(r => +r[recCol] >= target);
    // Wilson interval is stored only for the 0.90 threshold; do not pair it with rec85/rec95.
    if (hit) return { pf: +hit.perFactorN, rec: Math.round(+hit[recCol] * 100), lo: recCol === 'rec90' && hit.rec90_lo ? Math.round(+hit.rec90_lo * 100) : null, hi: recCol === 'rec90' && hit.rec90_hi ? Math.round(+hit.rec90_hi * 100) : null };
    const best = rows[rows.length - 1];
    return { fail: true, bestPf: +best.perFactorN, bestRec: Math.round(+best[recCol] * 100) };
  }

  console.log('Q-SYNTH design calculator');
  console.log('-------------------------');
  console.log(`Inputs: c=${cIn}, rho=${rhoIn}, M=${MIn} statements, K=${Kin} viewpoints, ${mode}, ${retention} retention, target recovery ${Math.round(target * 100)}% at |phi|>=${thr}.`);
  console.log(`Simulated grid: c in {0.4,0.6,0.8}, rho in {0,0.2,0.4}, M in {20,30,40,60,80}, persons/factor in {4,6,8,12}, K in {2,3,4,5,6,7}.`);
  if (M !== MIn) console.log(`(M snapped to ${M}.)`);
  if (kCeiling) console.log(`(K above the simulated ceiling; read at K=7, the magic-number-seven limit. Recovering more than seven distinct viewpoints is beyond what the table covers and is generally impractical.)`);
  else if (K !== Math.round(Kin)) console.log(`(K snapped to ${K}.)`);

  const printRec = (label, c, rho, r) => {
    if (!r) { console.log(`${label} (c=${c}, rho=${rho}): no matching condition.`); return; }
    if (r.fail) { console.log(`${label} (c=${c}, rho=${rho}): even ${r.bestPf}/factor reaches only ${r.bestRec}% -> design region is fragile.`); return; }
    const N = r.pf * K, recruit = Math.ceil(N / 0.8);
    console.log(`${label} (c=${c}, rho=${rho}): about ${r.pf} defining sorts/factor -> ${r.rec}%${r.lo != null ? ` [${r.lo},${r.hi}] 95% CI` : ''}; ${r.pf}x${K}=${N} defining sorts, recruit ~${recruit}+ participants (assuming ~80% of recruits define a factor).`);
  };

  console.log('');
  if (between) {
    console.log('Your c/rho fall between simulated grid points, so a range is shown:');
    printRec('  conservative', cCons, rhoCons, recommend(cCons, rhoCons));
    printRec('  optimistic  ', cOpt, rhoOpt, recommend(cOpt, rhoOpt));
    console.log('  Plan toward the conservative end unless you are confident your viewpoints are strong and distinct.');
  } else {
    printRec('RECOMMENDATION', cCons, rhoCons, recommend(cCons, rhoCons));
  }

  // statement-count note
  if (M < 40) console.log(`\nStatements: you specified ${MIn}; more statements make every by-person correlation steadier. Going to 40 statements typically lifts recovery noticeably at the margin.`);
  else console.log(`\nStatements: 40 is a solid count; each statement is an observation underlying every correlation, so do not go much lower.`);

  console.log('\nNote: this is planning guidance read against a simulated recovery table, not a fixed rule.');
  console.log('Estimate c (mean defining-sort loading squared) and rho (mean absolute correlation between factor arrays) from a pilot or a published study on a near topic. High communality lets you plan more modestly; overlapping viewpoints call for more.');
}

main();
