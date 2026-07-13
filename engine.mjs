// Self-contained public-math engine for Q-SYNTH.
//
// Clean-room implementations of the standard psychometric routines the recovery
// simulation needs: Pearson correlation, symmetric eigen-extraction (Jacobi) for
// principal components, Kaiser varimax rotation, Brown-weighted factor scores,
// and forced-distribution rendering. These are textbook algorithms, identical in
// substance to what KADE, PQMethod and qmethod compute; nothing here is specific
// to any closed product. Agreement with a real KADE analysis is verified in
// validation/ (correlation bit-exact, unrotated loadings and eigenvalues to
// floating point, varimax to a 4th-decimal rounding level).
//
// MIT-licensed. PCA uses a Jacobi symmetric eigensolver (not an SVD), so the file
// carries no third-party code.

// Banker's-neutral rounding matching the "round half away from zero" convention
// used across Q software for reported loadings and correlations.
export function evenRound(value, precision) {
  precision |= 0;
  const m = Math.pow(10, precision);
  value *= m;
  const sgn = (value > 0) | -(value < 0);
  const isHalf = value % 1 === 0.5 * sgn;
  const f = Math.floor(value);
  if (isHalf) value = f + (sgn > 0);
  return (isHalf ? value : Math.round(value)) / m;
}

// Pearson product-moment correlation; 0 if either vector is constant.
export function pearson(x, y) {
  const n = x.length;
  if (n === 0) return 0;
  let s1 = 0, s2 = 0, s1q = 0, s2q = 0, ps = 0;
  for (let i = 0; i < n; i++) { s1 += x[i]; s2 += y[i]; s1q += x[i] * x[i]; s2q += y[i] * y[i]; ps += x[i] * y[i]; }
  const num = ps - (s1 * s2) / n;
  const den = Math.sqrt((s1q - s1 * s1 / n) * (s2q - s2 * s2 / n));
  return den === 0 ? 0 : num / den;
}

// N x N by-person correlation matrix, correlations rounded to 5 decimals as Q
// software reports them.
export function correlationMatrix(sorts) {
  const n = sorts.length;
  const cm = Array.from({ length: n }, () => new Array(n));
  for (let i = 0; i < n; i++) {
    for (let k = i; k < n; k++) {
      const r = evenRound(pearson(sorts[i], sorts[k]), 5);
      cm[i][k] = r;
      if (k !== i) cm[k][i] = r;
    }
  }
  return cm;
}

export function matIdentity(n) {
  const I = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let i = 0; i < n; i++) I[i][i] = 1;
  return I;
}

// Classic cyclic Jacobi eigensolver for a real symmetric matrix. Returns
// eigenvalues and eigenvectors sorted by descending eigenvalue. eigenvectors[i][j]
// is the i-th component of the j-th eigenvector (columns are eigenvectors).
export function jacobi(A, maxIter = 200, tol = 1e-12) {
  const n = A.length;
  const S = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) S[i * n + j] = A[i][j];
  const V = new Float64Array(n * n);
  for (let i = 0; i < n; i++) V[i * n + i] = 1;

  for (let iter = 0; iter < maxIter; iter++) {
    let maxOff = 0, p = 0, q = 1;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      const a = Math.abs(S[i * n + j]);
      if (a > maxOff) { maxOff = a; p = i; q = j; }
    }
    if (maxOff < tol) break;
    const app = S[p * n + p], aqq = S[q * n + q], apq = S[p * n + q];
    let cos, sin;
    if (Math.abs(app - aqq) < tol * 1e-3) {
      const a = Math.PI / 4; cos = Math.cos(a); sin = Math.sin(a);
    } else {
      const tau = (aqq - app) / (2 * apq);
      const t = Math.sign(tau) / (Math.abs(tau) + Math.sqrt(1 + tau * tau));
      cos = 1 / Math.sqrt(1 + t * t); sin = t * cos;
    }
    for (let i = 0; i < n; i++) {
      if (i === p || i === q) continue;
      const sip = S[i * n + p], siq = S[i * n + q];
      S[i * n + p] = cos * sip - sin * siq; S[p * n + i] = S[i * n + p];
      S[i * n + q] = sin * sip + cos * siq; S[q * n + i] = S[i * n + q];
    }
    S[p * n + p] = cos * cos * app - 2 * sin * cos * apq + sin * sin * aqq;
    S[q * n + q] = sin * sin * app + 2 * sin * cos * apq + cos * cos * aqq;
    S[p * n + q] = 0; S[q * n + p] = 0;
    for (let i = 0; i < n; i++) {
      const vip = V[i * n + p], viq = V[i * n + q];
      V[i * n + p] = cos * vip - sin * viq;
      V[i * n + q] = sin * vip + cos * viq;
    }
  }
  const ev = new Array(n);
  for (let i = 0; i < n; i++) ev[i] = S[i * n + i];
  const idx = Array.from({ length: n }, (_, i) => i).sort((a, b) => ev[b] - ev[a]);
  const eigenvalues = idx.map(i => ev[i]);
  const eigenvectors = Array.from({ length: n }, () => new Array(n));
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) eigenvectors[i][j] = V[i * n + idx[j]];
  return { eigenvalues, eigenvectors };
}

// Principal-component extraction from a correlation matrix. loading[i][j] =
// eigenvector_ij * sqrt(eigenvalue_j). Sign of each factor is fixed by the
// standard inflection rule (flip if the directed sum of squares is negative), so
// the result is canonical regardless of the eigensolver's sign convention.
export function pcaExtract(corrMatrix, nFactors) {
  const n = corrMatrix.length;
  const { eigenvalues: allEigenvalues, eigenvectors: U } = jacobi(corrMatrix);
  const eigenvalues = allEigenvalues.slice(0, nFactors);
  const loadings = Array.from({ length: n }, () => new Array(nFactors));
  const inflection = [];
  for (let j = 0; j < nFactors; j++) {
    const sqrtEV = Math.sqrt(Math.max(0, eigenvalues[j]));
    let crit = 0;
    for (let i = 0; i < n; i++) {
      const load = evenRound(U[i][j] * sqrtEV, 8);
      loadings[i][j] = load;
      crit += evenRound(load * Math.abs(load), 8);
    }
    inflection.push(evenRound(crit, 8));
  }
  for (let j = 0; j < nFactors; j++) if (inflection[j] < 0) for (let i = 0; i < n; i++) loadings[i][j] = -loadings[i][j];
  const varianceExplained = eigenvalues.map(ev => ev / n);
  return { eigenvalues, allEigenvalues, loadings, varianceExplained };
}

// Kaiser varimax rotation with Kaiser (row) normalization. Returns rotated
// loadings and communalities.
export function varimax(loadings, maxIter = 100, tol = 1e-6) {
  const n = loadings.length, k = loadings[0].length;
  if (k < 2) {
    return { loadings: loadings.map(r => [...r]), communalities: loadings.map(r => r.reduce((s, v) => s + v * v, 0)) };
  }
  const comm = new Float64Array(n);
  for (let i = 0; i < n; i++) { let ss = 0; for (let j = 0; j < k; j++) ss += loadings[i][j] ** 2; comm[i] = Math.sqrt(ss); }
  const B = Array.from({ length: n }, (_, i) => { const h = comm[i]; return h === 0 ? new Array(k).fill(0) : loadings[i].map(v => v / h); });
  let prev = -Infinity;
  for (let iter = 0; iter < maxIter; iter++) {
    for (let p = 0; p < k - 1; p++) for (let q = p + 1; q < k; q++) {
      let a = 0, b = 0, c = 0, d = 0;
      for (let i = 0; i < n; i++) {
        const xi = B[i][p], yi = B[i][q];
        const ui = xi * xi - yi * yi, vi = 2 * xi * yi;
        a += ui; b += vi; c += ui * ui - vi * vi; d += 2 * ui * vi;
      }
      const num = d - 2 * a * b / n, den = c - (a * a - b * b) / n;
      const angle = (Math.abs(den) < 1e-15 && Math.abs(num) < 1e-15) ? 0 : Math.atan2(num, den) / 4;
      const cosA = Math.cos(angle), sinA = Math.sin(angle);
      for (let i = 0; i < n; i++) {
        const bp = B[i][p], bq = B[i][q];
        B[i][p] = cosA * bp + sinA * bq;
        B[i][q] = -sinA * bp + cosA * bq;
      }
    }
    let crit = 0;
    for (let j = 0; j < k; j++) {
      let sum = 0, sumSq = 0;
      for (let i = 0; i < n; i++) { const sq = B[i][j] ** 2; sum += sq; sumSq += sq * sq; }
      crit += sumSq / n - (sum / n) ** 2;
    }
    if (Math.abs(crit - prev) < tol) break;
    prev = crit;
  }
  const rotated = Array.from({ length: n }, (_, i) => B[i].map(v => v * comm[i]));
  return { loadings: rotated, communalities: Array.from(comm).map(h => h * h) };
}

// Brown-weighted, z-scored factor scores (the factor array) from one loading
// column. Weight of a defining sort is l/(1-l^2); scores standardized with the
// sample (N-1) standard deviation.
export function factorScores(sorts, loadings, threshold) {
  const P = sorts.length, M = sorts[0].length, t = threshold || 0;
  const scores = new Array(M).fill(0);
  const weights = [], def = [];
  for (let p = 0; p < P; p++) if (Math.abs(loadings[p]) >= t) { const l = loadings[p]; weights.push(l / (1 - l * l)); def.push(p); }
  if (def.length === 0) return scores;
  let absW = 0; for (const w of weights) absW += Math.abs(w);
  if (absW === 0) return scores;
  for (let m = 0; m < M; m++) { let w = 0; for (let i = 0; i < def.length; i++) w += weights[i] * sorts[def[i]][m]; scores[m] = w / absW; }
  let mean = 0; for (let m = 0; m < M; m++) mean += scores[m]; mean /= M;
  let ssq = 0; for (let m = 0; m < M; m++) ssq += (scores[m] - mean) ** 2;
  const sd = Math.sqrt(ssq / (M - 1));
  if (sd > 0) for (let m = 0; m < M; m++) scores[m] = (scores[m] - mean) / sd;
  return scores;
}

// Rank M scores into a forced distribution (ascending), returning grid positions.
export function generateCompositeSort(fScores, distribution) {
  const M = fScores.length;
  const composite = new Array(M);
  const indexed = fScores.map((score, idx) => ({ idx, score })).sort((a, b) => a.score - b.score);
  const K = distribution.length;
  const minVal = -Math.floor((K - 1) / 2);
  let si = 0;
  for (let col = 0; col < K; col++) {
    const val = minVal + col;
    for (let c = 0; c < distribution[col]; c++) { composite[indexed[si].idx] = val; si++; }
  }
  return composite;
}
