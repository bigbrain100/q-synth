// Ecological calibration — estimate realistic c (shared variance) and ρ (factor
// distinctness) from a PUBLISHED Q dataset, to show the sweep's parameter ranges
// (c∈{0.4,0.6,0.8}, ρ∈{0,0.4}) bracket reality. Public data only (firewall).
//
// Anchor: the classic Lipset study (Lipset 1963; Brown 1980), distributed as the
// `lipset` dataset in the qmethod R package. 9 Q-sorts × 33 statements, forced
// -4..+4, standard 3-factor solution.
import { readFileSync } from 'node:fs';
import { correlationMatrix, pcaExtract, varimax } from '../engine.mjs';
import { factorScores } from '../engine.mjs';
import { centroidExtract } from '../centroid.mjs';
import { zscore, tuckerPhi } from '../qsynth.mjs';

const { N, M, sorts } = JSON.parse(readFileSync(new URL('./lipset-sorts.json', import.meta.url)));
const K = 3;                                   // standard lipset solution
const sig = 2.58 / Math.sqrt(M);               // p<.01 significant loading threshold

function analyse(label, loadings) {
  const rot = varimax(loadings).loadings;      // N × K
  // flag each sort to its highest |loading| if significant (simple defining rule)
  const defs = Array.from({ length: K }, () => []);
  const communalities = [];
  rot.forEach((row, i) => {
    let bf = 0, bv = 0; row.forEach((v, f) => { if (Math.abs(v) > Math.abs(bv)) { bv = v; bf = f; } });
    communalities.push(row.reduce((s, v) => s + v * v, 0));
    if (Math.abs(bv) >= sig) defs[bf].push(bv * bv);          // defining loading²
  });
  // factor arrays (z-scored) and their pairwise distinctness
  const arrays = rot[0].map((_, f) => zscore(factorScores(sorts, rot.map(r => r[f]), 0)));
  const rho = [];
  for (let a = 0; a < K; a++) for (let b = a + 1; b < K; b++) rho.push(Math.abs(tuckerPhi(arrays[a], arrays[b])));
  const allDef = defs.flat();
  const meanC = allDef.reduce((s, v) => s + v, 0) / (allDef.length || 1);
  const meanComm = communalities.reduce((s, v) => s + v, 0) / N;
  console.log(`\n[${label}]`);
  console.log(`  defining sorts per factor : ${defs.map(d => d.length).join(', ')}  (sig |loading|≥${sig.toFixed(3)})`);
  console.log(`  c  (mean defining loading²): ${meanC.toFixed(3)}   [range ${Math.min(...allDef).toFixed(2)}–${Math.max(...allDef).toFixed(2)}]`);
  console.log(`  communality (mean, K=${K})   : ${meanComm.toFixed(3)}`);
  console.log(`  ρ  (|φ| between arrays)    : ${rho.map(v => v.toFixed(3)).join(', ')}   mean ${(rho.reduce((s, v) => s + v, 0) / rho.length).toFixed(3)}`);
}

const cm = correlationMatrix(sorts);
console.log(`Lipset ecological anchor — N=${N} sorts, M=${M} statements, K=${K}`);
analyse('PCA + varimax', pcaExtract(cm, K).loadings);
analyse('centroid + varimax', centroidExtract(cm, K).loadings);
