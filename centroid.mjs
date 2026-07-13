// Classic centroid factor extraction with iterated communalities
// (Brown 1980, Political Subjectivity, Appendix; Thompson 1980) — the
// traditional Q extraction method, as implemented by PQMethod / KADE.
//
// Public mathematics only, self-contained alongside engine.mjs (which provides
// the PCA path). Its purpose is a robustness check: to show the M×N recovery
// standard holds under the traditional centroid method, not just PCA.
//
// Fidelity: cross-checked against a real KADE centroid export (N=24, M=40) in
// validation/kade-crosscheck.mjs. Reproduces KADE's factor structure to a
// per-factor column correlation of ≥0.9999; extracted SS ("eigenvalues") land
// within ~2% of KADE's (a communality-convergence convention difference).
// Because factor recovery is scale-invariant (it depends on the loading PATTERN,
// not absolute SS), the residual 2% does not affect recovery outcomes.
//
// Algorithm:
//   diagonal ← communality estimates (initial: max |r| in each row), then for
//   each outer iteration:
//     for each factor on the running residual matrix:
//       1. sign-reflect variables (one at a time, most-negative column first)
//          until every reflected column sum ≥ 0 — maximizes the grand sum
//       2. loading_j = sign_j · (reflected column sum_j) / sqrt(grand sum)
//       3. residualize: R[i][j] -= a_i · a_j
//     recompute communalities h²_i = Σ_f loading_if², put on the diagonal, repeat
//   until communalities converge (or maxIter).

/**
 * @param {number[][]} corr — N×N correlation matrix (off-diagonal correlations)
 * @param {number} nFactors — number of centroid factors to extract
 * @param {{maxIter?: number, tol?: number}} [opts]
 * @returns {{ loadings: number[][], eigenvalues: number[], communalities: number[], iterations: number }}
 */
export function centroidExtract(corr, nFactors, opts = {}) {
  const { maxIter = 50, tol = 1e-6 } = opts;
  const n = corr.length;

  // initial communality estimate: highest absolute off-diagonal correlation per row
  let h = new Array(n);
  for (let i = 0; i < n; i++) { let mx = 0; for (let j = 0; j < n; j++) if (j !== i) mx = Math.max(mx, Math.abs(corr[i][j])); h[i] = mx; }

  let loadings, eigenvalues, iterations = 0;
  for (let it = 0; it < maxIter; it++) {
    iterations = it + 1;
    const R = corr.map(r => r.slice());
    for (let i = 0; i < n; i++) R[i][i] = h[i];

    loadings = Array.from({ length: n }, () => new Array(nFactors).fill(0));
    eigenvalues = new Array(nFactors).fill(0);

    for (let f = 0; f < nFactors; f++) {
      const sign = new Array(n).fill(1);
      const colSum = j => { let s = 0; for (let i = 0; i < n; i++) s += sign[i] * sign[j] * R[i][j]; return s; };
      for (let guard = 0; guard < 10000; guard++) {
        let worst = -1, worstVal = 0;
        for (let j = 0; j < n; j++) { const s = colSum(j); if (s < worstVal) { worstVal = s; worst = j; } }
        if (worst === -1) break;
        sign[worst] *= -1;
      }
      let G = 0; const c = new Array(n);
      for (let j = 0; j < n; j++) { c[j] = colSum(j); G += c[j]; }
      // Degenerate residual: no common variance left to extract on a centroid
      // axis (grand sum ≤ 0). Emit zero loadings for this and remaining factors
      // rather than sqrt(negative) → NaN. Happens at small N with high shared
      // variance once the reduced matrix is near rank-deficient.
      if (!(G > 1e-9)) { for (let g = f; g < nFactors; g++) { for (let j = 0; j < n; j++) loadings[j][g] = 0; eigenvalues[g] = 0; } break; }
      const sqrtG = Math.sqrt(G);
      const a = new Array(n);
      for (let j = 0; j < n; j++) { a[j] = sign[j] * (c[j] / sqrtG); loadings[j][f] = a[j]; }
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) R[i][j] -= a[i] * a[j];
      let ss = 0; for (let j = 0; j < n; j++) ss += a[j] * a[j];
      eigenvalues[f] = ss;
    }

    // update communalities and test convergence; cap at 1.0 to prevent Heywood
    // overshoot (h² > 1) driving the reduced matrix non-positive-definite.
    let maxDelta = 0;
    const hNew = new Array(n);
    for (let i = 0; i < n; i++) { let s = 0; for (let f = 0; f < nFactors; f++) s += loadings[i][f] * loadings[i][f]; s = Math.min(1, s); hNew[i] = s; maxDelta = Math.max(maxDelta, Math.abs(s - h[i])); }
    h = hNew;
    if (maxDelta < tol) break;
  }

  return { loadings, eigenvalues, communalities: h, iterations };
}
