// Tests for the self-contained numerical engine.
// Run with: node --test
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { correlationMatrix, pearson, pcaExtract, varimax, factorScores, generateCompositeSort } from '../engine.mjs';

const approx = (a, b, eps = 1e-6) => Math.abs(a - b) <= eps;

test('pearson: perfect correlation and anti-correlation', () => {
  assert.ok(approx(pearson([1, 2, 3, 4], [2, 4, 6, 8]), 1));
  assert.ok(approx(pearson([1, 2, 3, 4], [4, 3, 2, 1]), -1));
  assert.equal(pearson([1, 1, 1], [1, 2, 3]), 0); // zero variance -> 0
});

test('correlationMatrix: symmetric with unit diagonal', () => {
  const sorts = [[2, -1, 0, 1, -2], [1, 0, -1, 2, -2], [-2, 1, 0, -1, 2]];
  const cm = correlationMatrix(sorts);
  for (let i = 0; i < 3; i++) {
    assert.ok(approx(cm[i][i], 1, 1e-9));
    for (let j = 0; j < 3; j++) assert.ok(approx(cm[i][j], cm[j][i], 1e-12));
  }
});

test('pcaExtract: eigenvalues sum to N and equal column sums of squared loadings', () => {
  const sorts = Array.from({ length: 6 }, (_, i) =>
    [0, 1, 2, 3, 4, 5, 6, 7].map(j => Math.sin(i + j) + 0.3 * Math.cos(2 * j - i)));
  const cm = correlationMatrix(sorts);
  const { eigenvalues, allEigenvalues, loadings } = pcaExtract(cm, 6);
  const N = sorts.length;
  assert.ok(approx(allEigenvalues.reduce((a, b) => a + b, 0), N, 1e-4), 'trace = N');
  // eigenvalues descending
  for (let j = 1; j < eigenvalues.length; j++) assert.ok(eigenvalues[j] <= eigenvalues[j - 1] + 1e-9);
  // column sum of squared loadings == eigenvalue (PCA property)
  for (let j = 0; j < 6; j++) {
    const ss = loadings.reduce((s, row) => s + row[j] * row[j], 0);
    assert.ok(approx(ss, eigenvalues[j], 1e-3), `col ${j} SS = eigenvalue`);
  }
});

test('varimax: preserves per-row communalities', () => {
  const sorts = Array.from({ length: 7 }, (_, i) =>
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9].map(j => Math.cos(i * 0.7 + j) + 0.2 * (j % 3)));
  const cm = correlationMatrix(sorts);
  const { loadings } = pcaExtract(cm, 3);
  const { loadings: rot } = varimax(loadings);
  for (let i = 0; i < loadings.length; i++) {
    const hBefore = loadings[i].reduce((s, v) => s + v * v, 0);
    const hAfter = rot[i].reduce((s, v) => s + v * v, 0);
    assert.ok(approx(hBefore, hAfter, 1e-6), `row ${i} communality preserved`);
  }
});

test('factorScores: returns a z-scored array (mean 0, sd 1)', () => {
  const sorts = [[4, 3, 1, -1, -2, -3, -1, 0, 1, -2], [3, 4, 0, -2, -1, -3, 0, 1, 1, -3],
    [-3, -2, 0, 2, 3, 4, 1, 0, -1, 1], [-2, -3, 1, 3, 2, 4, 0, -1, -1, 0]];
  const loadings = [0.8, 0.75, -0.7, -0.72];
  const fs = factorScores(sorts, loadings, 0);
  const M = fs.length;
  const mean = fs.reduce((a, b) => a + b, 0) / M;
  const sd = Math.sqrt(fs.reduce((s, v) => s + (v - mean) ** 2, 0) / (M - 1));
  assert.ok(approx(mean, 0, 1e-9), 'mean 0');
  assert.ok(approx(sd, 1, 1e-9), 'sd 1');
});

test('generateCompositeSort: output respects the forced distribution', () => {
  const dist = [1, 2, 3, 2, 1]; // sums to 9, grid -2..+2
  const scores = [0.5, -0.3, 1.2, -1.1, 0.0, 0.9, -0.7, 0.2, -0.9];
  const comp = generateCompositeSort(scores, dist);
  assert.equal(comp.length, 9);
  // count of each grid value must match the distribution
  const counts = {};
  for (const v of comp) counts[v] = (counts[v] || 0) + 1;
  assert.deepEqual([counts[-2], counts[-1], counts[0], counts[1], counts[2]], dist);
});
