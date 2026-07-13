// Tests for the generator, centroid extraction, recovery, and Wilson interval.
// Run with: node --test
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { makeRng, buildDistribution, makeIdeals, makeParticipants, recover, wilson, tuckerPhi, zscore } from '../qsynth.mjs';
import { centroidExtract } from '../centroid.mjs';
import { correlationMatrix } from '../engine.mjs';

const approx = (a, b, eps = 1e-6) => Math.abs(a - b) <= eps;

test('makeRng: deterministic given a seed', () => {
  const a = makeRng(42), b = makeRng(42);
  for (let i = 0; i < 100; i++) assert.equal(a.u(), b.u());
});

test('buildDistribution: sums to M with every column at least 1', () => {
  for (const M of [20, 33, 40, 60]) {
    const d = buildDistribution(M);
    assert.equal(d.reduce((a, b) => a + b, 0), M);
    for (const c of d) assert.ok(c >= 1);
  }
});

test('makeParticipants: correct counts and valid forced sorts', () => {
  const rng = makeRng(1), K = 3, M = 40, pf = 6, dist = buildDistribution(M);
  const idealZ = makeIdeals(rng, K, M, 0.2, dist);
  assert.equal(idealZ.length, K);
  assert.equal(idealZ[0].length, M);
  const { sorts, trueType } = makeParticipants(rng, idealZ, pf, 0.6, 'forced', dist, 0);
  assert.equal(sorts.length, K * pf);
  assert.equal(trueType.length, K * pf);
  // every forced sort is a permutation into the same grid, so equal column-sum
  const sums = sorts.map(s => s.reduce((a, b) => a + b, 0));
  for (const s of sums) assert.ok(approx(s, sums[0], 1e-9), 'forced sorts share the grid total');
});

test('makeParticipants: non-defining sorts are added', () => {
  const rng = makeRng(2), dist = buildDistribution(40);
  const idealZ = makeIdeals(rng, 3, 40, 0, dist);
  const { sorts, trueType } = makeParticipants(rng, idealZ, 6, 0.6, 'forced', dist, 0.2);
  assert.ok(sorts.length > 18, 'has extra non-defining sorts');
  assert.ok(trueType.includes(-1), 'non-defining sorts flagged as -1');
});

test('centroidExtract: finite non-negative eigenvalues; communalities in [0,1]', () => {
  // Note: unlike PCA, the centroid method does not guarantee strictly descending
  // factor sums of squares, so we do not assert ordering here.
  const rng = makeRng(3), dist = buildDistribution(40);
  const idealZ = makeIdeals(rng, 3, 40, 0.3, dist);
  const { sorts } = makeParticipants(rng, idealZ, 8, 0.7, 'forced', dist, 0);
  const cm = correlationMatrix(sorts);
  const { loadings, eigenvalues, communalities } = centroidExtract(cm, 3);
  for (const row of loadings) for (const v of row) assert.ok(Number.isFinite(v), 'loading finite (no NaN)');
  for (const e of eigenvalues) assert.ok(Number.isFinite(e) && e >= -1e-9, 'eigenvalue finite and >=0');
  for (const h of communalities) assert.ok(h >= -1e-9 && h <= 1 + 1e-9, 'communality in [0,1]');
});

test('centroid stays finite in the small-N high-c corner (regression for the NaN guard)', () => {
  const rng = makeRng(99), dist = buildDistribution(40);
  for (let t = 0; t < 20; t++) {
    const idealZ = makeIdeals(rng, 3, 40, 0.4, dist);
    const { sorts } = makeParticipants(rng, idealZ, 4, 0.8, 'forced', dist, 0); // N=12, c=0.8
    const { loadings, eigenvalues } = centroidExtract(correlationMatrix(sorts), 3);
    assert.ok(loadings.flat().every(Number.isFinite), 'no NaN loadings');
    assert.ok(eigenvalues.every(Number.isFinite), 'no NaN eigenvalues');
  }
});

test('recover: returns estK factors and congruences in [0,1]', () => {
  const rng = makeRng(4), K = 3, M = 40, dist = buildDistribution(M);
  const idealZ = makeIdeals(rng, K, M, 0.2, dist);
  const { sorts } = makeParticipants(rng, idealZ, 8, 0.7, 'forced', dist, 0);
  for (const ext of ['pca', 'centroid']) {
    const { estK, phis } = recover(sorts, idealZ, K, ext);
    assert.equal(estK, K);
    assert.equal(phis.length, K);
    for (const p of phis) assert.ok(p >= 0 && p <= 1 + 1e-9, `${ext} congruence in [0,1]`);
  }
});

test('recover: strong distinct viewpoints are recovered well', () => {
  const rng = makeRng(5), K = 3, M = 40, dist = buildDistribution(M);
  const idealZ = makeIdeals(rng, K, M, 0.0, dist);
  const { sorts } = makeParticipants(rng, idealZ, 8, 0.8, 'forced', dist, 0);
  const { phis } = recover(sorts, idealZ, K, 'pca');
  assert.ok(Math.min(...phis) > 0.9, 'high communality + distinct -> recovered');
});

test('tuckerPhi and zscore basics', () => {
  const a = [1, 2, 3, 4], b = [2, 4, 6, 8];
  assert.ok(approx(tuckerPhi(a, b), 1), 'proportional vectors -> congruence 1');
  const z = zscore([1, 2, 3, 4, 5]);
  const mean = z.reduce((x, y) => x + y, 0) / z.length;
  assert.ok(approx(mean, 0, 1e-9), 'zscore mean 0');
});

test('wilson: interval brackets the point estimate within [0,1]', () => {
  for (const [x, n] of [[0, 500], [250, 500], [500, 500], [1, 500], [499, 500]]) {
    const [lo, hi] = wilson(x, n);
    const p = x / n;
    assert.ok(lo >= 0 && hi <= 1, 'within [0,1]');
    assert.ok(lo <= p + 1e-9 && p <= hi + 1e-9, 'brackets p');
  }
});
