# Q-SYNTH

[![DOI](https://zenodo.org/badge/1298889650.svg)](https://doi.org/10.5281/zenodo.21429739)

**A synthetic Q-sort recovery simulator, and a G\*Power-style design calculator, for evidence-based sample-size (N persons) and item-count (M statements) guidance in Q methodology.**

Q-SYNTH plants known viewpoints as factor arrays, generates synthetic Q sorts under controlled communality and viewpoint distinctness, runs a standard by-person factor analysis (principal components or centroid, with varimax rotation, factor arrays computed from auto-flagged defining sorts exactly as mainstream Q software computes them), and measures how reliably the planted structure is recovered. It turns Q's untested rules of thumb ("40–80 statements", "4–5 defining sorts per factor") into condition-specific, confidence-interval-backed guidance, and its calculator turns that guidance into a concrete recommendation for a planned study.

It is a self-contained research tool with no dependency on any closed analysis engine. Its numerical routines are open implementations of the standard Q factor-analysis algorithms, cross-validated against KADE (see `validation/`).

## Requirements

- [Node.js](https://nodejs.org) version 18 or newer. No other runtime dependency is required for the simulation or the calculator. (`xlsx` is an optional dependency used only by the KADE cross-check in `validation/`.)

## Installation

```bash
git clone https://github.com/bigbrain100/q-synth.git
cd q-synth
```

The simulation, the calculator, and the tests run with Node alone: **no dependencies and no build step**. (The optional KADE cross-check needs a spreadsheet reader; see below.)

## Usage

### Design calculator (plan a study)

**In a browser (no install):** use the live tool at <https://bigbrain100.github.io/q-synth/calculator.html>, or open `calculator.html` locally in any web browser. A step-by-step tutorial written for researchers with no factor-analysis background is in [TUTORIAL.md](TUTORIAL.md) (English) and [TUTORIAL.zh-Hant.md](TUTORIAL.zh-Hant.md) (繁體中文). It provides a G*Power-style form (expected communality, viewpoint overlap, statements, target recovery) and returns the recommended defining sorts per factor and participant count, with a recovery-rate curve. The page is self-contained (the lookup table is embedded), so it works by double-clicking the file, and it can be deployed (e.g. GitHub Pages) as a live tool.

**On the command line:** given the kind of viewpoints you expect, get the sample size and item count needed for the factors to be recovered:

```bash
node calculate.mjs --c 0.6 --rho 0.2 --M 40 --K 3 --target 0.90
```

`--c` expected communality of defining sorts (0–1), `--rho` expected correlation between viewpoints (0–1; higher means viewpoints overlap more and are harder to separate), `--M` planned number of statements, `--target` desired recovery rate. Estimate `c` as the mean squared loading of defining sorts, and `rho` as the mean absolute correlation between factor arrays, from a pilot or a published study on a near topic.

The calculator reads a **coarse simulated lookup table** (c ∈ {0.4, 0.6, 0.8}, rho ∈ {0, 0.2, 0.4}, M ∈ {20, 30, 40, 60, 80}, defining sorts per factor ∈ {4, 6, 8, 12}, K ∈ {2..7} expected viewpoints). Inputs are snapped to the nearest grid point, leaning to the safer side; when they fall between grid points the calculator prints both a conservative and an optimistic recommendation. It is planning guidance to read against your topic, not a continuous power function.

### Run the simulation sweep

```bash
node run.mjs 1000             # parallel (worker_threads), 1000 reps/cell -> results.csv + figures/
node robustness.mjs 1000      # PCA vs centroid, paired -> robustness.csv
node seeding/seed-calibrate.mjs   # ecological calibration on the public Lipset data
```

All runs are deterministic (seeded). Swept factors: statements (M), defining sorts per factor (→ N), communality (c), distinctness (rho), forced/free distribution, share of non-defining sorts, and factor-retention rule. Factor arrays are computed from **auto-flagged defining sorts** (|loading| > 1.96/√M and a majority of the sort's common variance — the auto-flag rule of the `qmethod` package and KADE); recovery means every planted viewpoint is matched by a distinct recovered factor (optimal one-to-one assignment) at |Tucker φ| ≥ 0.85 / 0.90 / 0.95, reported with Wilson 95% confidence intervals.

### Validate against KADE (optional)

```bash
npm install xlsx                              # spreadsheet reader, only for this step
node validation/kade-crosscheck.mjs <KADE_results.xlsx>
```

`xlsx` is not a project dependency (a normal install pulls nothing); install it only when you want to run this optional cross-check. You supply your own KADE export. The correlation matrix, unrotated PCA loadings, eigenvalues, and varimax loadings are compared to KADE's; centroid structure is compared by per-factor congruence.

## Testing

```bash
npm test          # node --test
```

The tests pin down the numerical behavior of the engine (correlation, PCA, varimax, factor scores, forced rendering), the centroid extractor (including the small-N high-communality regression), the generator, the recovery metric, and the Wilson interval.

## What Q-SYNTH is for

The exploratory-factor-analysis sample-size literature (MacCallum et al., 1999) established that required sample size is not fixed: it depends on communality and how well each factor is determined. That paradigm had never been transferred to Q's by-person factoring, where the roles invert: **M statements** are the paired observations that steady each person-to-person correlation, and **N persons** are the variables being factored (so defining sorts per factor is the overdetermination, Brown's "4–5 per factor"). In R-mode designs this kind of a-priori planning is served by tools like G\*Power; Q had no equivalent. Q-SYNTH provides one.

## Documentation

- `TUTORIAL.md` / `TUTORIAL.zh-Hant.md` — step-by-step calculator tutorial (English / 繁體中文), no factor-analysis background assumed.
- `DATA_GENERATION.md` — how the synthetic data are generated and how to reproduce every result.
- `ROBUSTNESS.md` — PCA vs centroid robustness of the design standard.
- `PRIOR_ART.md` — how Q-SYNTH relates to existing Q tools and the sample-size literature.
- `validation/README.md` — the KADE cross-check.
- `paper.md` — the accompanying software paper.

## Citation

If you use Q-SYNTH, please cite it (all-versions DOI: [10.5281/zenodo.21429739](https://doi.org/10.5281/zenodo.21429739)):

> Lee, M.-S., Wang, W.-H., & Tsai, Y.-J. (2026). *Q-SYNTH: A recovery-simulation tool for sample-size and item-count guidance in Q methodology* (Version 1.1.0) [Computer software]. https://doi.org/10.5281/zenodo.21429739

or use GitHub's **Cite this repository** button. An accompanying software paper is included in the repository (`paper.md`).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) and the [Code of Conduct](CODE_OF_CONDUCT.md).

## License

MIT. See [LICENSE](LICENSE).
