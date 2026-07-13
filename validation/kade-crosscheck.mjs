// KADE cross-check — validates q-synth's self-contained public engine
// (engine.mjs, centroid.mjs) against REAL KADE analysis exports. Auto-detects
// geometry (N sorts, M statements, extracted & retained factor counts) and the
// extraction method (PCA vs centroid), so one script covers KADE exports of
// different sizes and methods.
//
// IRB firewall: reads ONLY numeric matrices; reports ONLY engine-vs-KADE
// agreement metrics — never statement text or any participant viewpoint. The
// KADE .xlsx files are human-subjects data and are NOT committed; pass paths as
// arguments.
//
// Usage:  node kade-crosscheck.mjs <KADE_results.xlsx> [more.xlsx ...]
import { createRequire } from 'module';
// xlsx is an optional dependency (see package.json); install it with `npm install`.
const require = createRequire(import.meta.url);
let xlsx;
try { xlsx = require('xlsx'); }
catch { console.error("This cross-check needs the 'xlsx' reader. Run `npm install` first."); process.exit(1); }
import { correlationMatrix, pcaExtract, varimax } from '../engine.mjs';
import { centroidExtract } from '../centroid.mjs';

const files = process.argv.slice(2);
if (!files.length) { console.error('usage: node kade-crosscheck.mjs <KADE_results.xlsx> [...]'); process.exit(1); }

const V = (ws, R, C) => { const c = ws[xlsx.utils.encode_cell({ r: R, c: C })]; return c ? c.v : null; };
const rng = ws => xlsx.utils.decode_range(ws['!ref']);
const findRow = (ws, label, col = 0) => { const r = rng(ws); for (let R = r.s.r; R <= r.e.r; R++) if (String(V(ws, R, col)).trim() === label) return R; return -1; };
const maxAbsDiff = (A, B) => { let m = 0, at = ''; for (let i = 0; i < A.length; i++) for (let j = 0; j < A[0].length; j++) { const d = Math.abs(A[i][j] - B[i][j]); if (d > m) { m = d; at = `[${i}][${j}] ours=${A[i][j]} kade=${B[i][j]}`; } } return { m, at }; };

function crosscheck(file) {
  const wb = xlsx.readFile(file);
  const qs = wb.Sheets['Q sorts']; const qr = rng(qs); const hdr = findRow(qs, '参与者');
  const stmtCols = []; for (let C = 1; C <= qr.e.c; C++) if (/^s?\d+$/i.test(String(V(qs, hdr, C)).trim())) stmtCols.push(C);
  const M = stmtCols.length;
  const sorts = [];
  for (let R = hdr + 1; R <= qr.e.r; R++) { const id = V(qs, R, 0); if (id == null || id === '') break; sorts.push(stmtCols.map(C => Number(V(qs, R, C)))); }
  const N = sorts.length;

  const uw = wb.Sheets['未旋转的因子矩阵']; const uh = findRow(uw, 'Nm');
  const facCols = []; for (let C = 0; C <= rng(uw).e.c; C++) if (/因子/.test(String(V(uw, uh, C)))) facCols.push(C);
  const Kext = facCols.length;
  const fw = wb.Sheets['因子负荷']; const fh = findRow(fw, 'Nm');
  const loadCols = []; for (let C = 2; C <= rng(fw).e.c; C++) if (/因子/.test(String(V(fw, fh, C)))) loadCols.push(C);
  const Kret = loadCols.length;

  // [A] correlation matrix (KADE ×100 int)
  const cm = correlationMatrix(sorts);
  const cw = wb.Sheets['相关矩阵']; const ch = findRow(cw, '参与者');
  const kadeCorr = sorts.map((_, i) => sorts.map((__, j) => Number(V(cw, ch + 1 + i, 1 + j))));
  const A = maxAbsDiff(cm.map(r => r.map(v => Math.round(v * 100))), kadeCorr);

  // extraction method: PCA factor-1 SS is the maximal single-factor variance.
  // If KADE's factor-1 SS is materially below our PCA factor-1 SS → centroid/PAF.
  const pca = pcaExtract(cm, Kext);
  const eigRow = findRow(uw, '特征值', 1);
  const kadeEig = facCols.map(C => Number(V(uw, eigRow, C)));
  const method = Math.abs(pca.eigenvalues[0] - kadeEig[0]) < 1e-3 ? 'PCA' : (kadeEig[0] < pca.eigenvalues[0] ? 'centroid/PAF' : 'unknown');

  console.log(`\n### ${file.split(/[\/]/).pop()}  —  N=${N}, M=${M}, extracted=${Kext}, retained=${Kret}, method=${method}`);
  console.log(`[A] correlation ×100 int       maxAbsDiff=${A.m}   ${A.m === 0 ? 'EXACT ✓' : '✗ ' + A.at}`);

  if (method !== 'PCA') {
    // KADE used the traditional centroid method → cross-check q-synth's centroid
    // (public math; iterated communalities). Factor recovery is scale-invariant,
    // so the metric that matters is per-factor pattern agreement (column corr).
    const kadeUn = sorts.map((_, i) => facCols.map(C => Number(V(uw, uh + 1 + i, C))));
    const cen = centroidExtract(cm, Kext);
    const al = cen.loadings.map(r => r.slice());
    for (let j = 0; j < Kext; j++) { let d = 0; for (let i = 0; i < N; i++) d += cen.loadings[i][j] * kadeUn[i][j]; const s = d < 0 ? -1 : 1; for (let i = 0; i < N; i++) al[i][j] *= s; }
    let loadMax = 0; const colCorr = [];
    for (let j = 0; j < Kext; j++) {
      let sx = 0, sy = 0, sxy = 0, sxx = 0, syy = 0;
      for (let i = 0; i < N; i++) { const x = al[i][j], y = kadeUn[i][j]; loadMax = Math.max(loadMax, Math.abs(x - y)); sx += x; sy += y; sxy += x * y; sxx += x * x; syy += y * y; }
      colCorr.push((N * sxy - sx * sy) / Math.sqrt((N * sxx - sx * sx) * (N * syy - sy * sy)));
    }
    const eigD = Math.max(...cen.eigenvalues.map((v, j) => Math.abs(v - kadeEig[j])));
    const okPattern = Math.min(...colCorr) > 0.999;
    console.log(`[B] centroid loadings (pattern) per-factor col-corr = ${colCorr.map(v => v.toFixed(5)).join(', ')}   ${okPattern ? '✓ (structure reproduced)' : '✗'}`);
    console.log(`    loadings maxAbsDiff=${loadMax.toFixed(4)}   eigenvalues maxAbsDiff=${eigD.toFixed(4)} (~${(eigD / kadeEig[0] * 100).toFixed(1)}% of F1; communality-convergence convention, scale-invariant for recovery)`);
    console.log(`    our centroid eig : ${cen.eigenvalues.map(v => v.toFixed(4)).join(', ')}  (converged in ${cen.iterations} iters)`);
    console.log(`    KADE centroid eig: ${kadeEig.map(v => v.toFixed(4)).join(', ')}`);
    return;
  }

  // [B] unrotated PCA loadings + eigenvalues
  const kadeUn = sorts.map((_, i) => facCols.map(C => Number(V(uw, uh + 1 + i, C))));
  const al = pca.loadings.map(r => r.slice());
  for (let j = 0; j < Kext; j++) { let d = 0; for (let i = 0; i < N; i++) d += pca.loadings[i][j] * kadeUn[i][j]; const s = d < 0 ? -1 : 1; for (let i = 0; i < N; i++) al[i][j] *= s; }
  const B = maxAbsDiff(al.map(r => r.map(v => Number(v.toFixed(4)))), kadeUn);
  const eigD = Math.max(...pca.eigenvalues.map((v, j) => Math.abs(v - kadeEig[j])));
  console.log(`[B] unrotated PCA loadings 4dp maxAbsDiff=${B.m}   ${B.m <= 1e-4 ? '✓' : '✗ ' + B.at}`);
  console.log(`    eigenvalues                maxAbsDiff=${eigD.toExponential(2)}   ${eigD < 1e-4 ? '✓' : '✗'}`);

  // [C] varimax on RETAINED factors
  const rot = varimax(pcaExtract(cm, Kret).loadings).loadings;
  const kadeRot = sorts.map((_, i) => loadCols.map(C => Number(V(fw, fh + 1 + i, C))));
  const used = new Set(); const cmap = [];
  for (let jo = 0; jo < Kret; jo++) { let best = -1, bv = -1, bs = 1; for (let jk = 0; jk < Kret; jk++) { if (used.has(jk)) continue; let d = 0; for (let i = 0; i < N; i++) d += rot[i][jo] * kadeRot[i][jk]; if (Math.abs(d) > bv) { bv = Math.abs(d); best = jk; bs = d < 0 ? -1 : 1; } } used.add(best); cmap.push({ jo, jk: best, sign: bs }); }
  const oursRe = rot.map(row => { const o = new Array(Kret); cmap.forEach(m => o[m.jk] = Number((row[m.jo] * m.sign).toFixed(4))); return o; });
  const C = maxAbsDiff(oursRe, kadeRot);
  let over = 0; for (let i = 0; i < N; i++) for (let j = 0; j < Kret; j++) if (Math.abs(oursRe[i][j] - kadeRot[i][j]) > 1e-4) over++;
  console.log(`[C] varimax loadings 4dp       maxAbsDiff=${C.m.toFixed(4)}   ${C.m <= 5e-4 ? '✓ rounding-level' : '✗ ' + C.at}   (${over}/${N * Kret} cells >1e-4)`);
}
for (const f of files) {
  try { crosscheck(f); }
  catch (e) {
    console.error(`\nCould not process "${f}": ${e.code === 'ENOENT' ? 'file not found' : e.message}`);
    console.error('Provide the path to a KADE results .xlsx export you have on disk.');
    process.exitCode = 1;
  }
}
