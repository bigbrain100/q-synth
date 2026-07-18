---
title: 'Q-SYNTH: A recovery-simulation tool for sample-size and item-count guidance in Q methodology'
tags:
  - Q methodology
  - factor analysis
  - Monte Carlo simulation
  - research design
  - sample size
  - JavaScript
authors:
  - name: Ming-Shinn Lee
    orcid: 0000-0002-2932-7897
    corresponding: true
    affiliation: 1
  - name: Wan-Hsiang Wang
    orcid: 0000-0002-3011-5699
    affiliation: "1, 2"
  - name: Yi-Jung Tsai
    orcid: 0009-0009-4124-0493
    affiliation: 1
affiliations:
  - name: Department of Education and Human Potentials Development, Hua-Shih College of Education, National Dong Hwa University, Hualien, Taiwan
    index: 1
  - name: Department of Nursing, Hualien Tzu Chi Hospital, Hualien, Taiwan
    index: 2
date: 16 July 2026
bibliography: paper.bib
---

# Summary

Q methodology studies human subjectivity by having participants rank-order a set of statements into a forced quasi-normal grid, a Q sort, and then correlating and factor-analyzing the *people* rather than the items [@Stephenson1953; @Brown1980; @McKeownThomas2013]. The factors that emerge are shared points of view. A perennial and practical question in designing such a study is how many statements and how many participants it needs. In practice this is settled by rules of thumb, roughly 40 to 80 statements and 4 to 5 clearly loading sorts per factor [@WattsStenner2012; @Brown1980], whose adequacy has never been tested against the criterion that matters most: whether the intended shared viewpoints are actually recovered.

`Q-SYNTH` is a small, self-contained simulator that supplies that missing test. It plants known viewpoints as factor arrays, generates synthetic Q sorts under a controlled level of shared variance (communality) and viewpoint distinctness, runs a standard by-person factor analysis, and measures how often the planted structure is recovered. It sweeps the design space, the number of statements, the number of defining sorts per factor, communality, distinctness, the forced or free distribution, the share of non-defining participants, and the factor-retention rule, and reports recovery rates with Wilson confidence intervals. This turns untested conventions into condition-specific, evidence-based guidance for planning a Q study.

# Statement of need

Existing Q software addresses a different, later stage of the workflow, and `Q-SYNTH` is meant to complement rather than replace it. `PQMethod` [@Schmolck2014], KADE [@Banasick2019], and the R package `qmethod` [@Zabala2014] analyze data that have already been collected: they extract and rotate factors and compute factor arrays. Bootstrap resampling within `qmethod` [@ZabalaPascual2016] quantifies the stability of a solution once the data are in hand. These are diagnostics on data one already has. They do not speak to the question a researcher faces *before* collecting anything: given the kind of viewpoints I expect, how many statements and participants should I plan for so that those viewpoints will be recoverable?

The equivalent question in ordinary (R-mode) factor analysis was long ago answered by simulation, showing that the required sample size is not a fixed number but depends on the communality of the variables and how well each factor is determined [@MacCallum1999]. In R-mode designs, a-priori planning of this kind is routinely supported by dedicated power-analysis software such as `G*Power` [@Faul2007]; Q methodology has had no equivalent. `Q-SYNTH` provides one, transposed to Q's by-person matrix, where the participants are the variables and the statements are the paired observations underlying every correlation between two sorts. In the spirit of `G*Power`, its calculator takes the communality and viewpoint overlap a researcher expects and returns the defining sorts per factor and the participant count needed to reach a target recovery rate. It is a planning aid backed by a simulated recovery table rather than a continuous power function: it reads a coarse grid (communality, distinctness, statement count, defining sorts per factor, and the number of viewpoints from two up to the magic-number-seven ceiling), snaps inputs to the nearest simulated condition, and reports a conservative and an optimistic recommendation when they fall between grid points. It is intended for Q researchers planning studies, for methodologists studying the behavior of Q factor solutions, and for teaching, and it makes the design trade-offs concrete and reproducible.

# State of the field

Two things distinguish `Q-SYNTH` from the tools above and explain why it is a new package rather than a contribution to one of them. First, it targets the design stage, not the analysis stage: it generates and analyzes *synthetic* sorts with a known planted structure, which lets it measure recovery against ground truth, something no analysis of a single real study can do. Second, it is deliberately self-contained. Rather than depend on `PQMethod`, KADE, or `qmethod`, it re-implements the standard extraction and rotation routines as open code and validates them against KADE, so the simulator has no external analysis dependency, runs anywhere Node.js runs, and can later be ported to the browser for interactive use. The `qmethod` ecosystem remains the natural home for *analyzing* Q data; `Q-SYNTH` fills the separate, previously empty niche of a-priori design planning.

# Functionality

`Q-SYNTH` runs on Node.js with no dependency on any closed or external analysis engine. Its numerical core (`engine.mjs`) provides open implementations of the standard routines: Pearson correlation, a Jacobi symmetric eigensolver for principal-component extraction, Kaiser varimax rotation, Brown-weighted factor scores, and forced-distribution rendering; `centroid.mjs` adds classic centroid extraction with iterated communalities. Both are verified against a real analysis exported from KADE [@Banasick2019]. For principal-component extraction the agreement is exact: the correlation matrix matches to the last reported digit, and the unrotated loadings, eigenvalues, and varimax loadings match to floating-point or fourth-decimal rounding. For centroid extraction, which is convention-dependent across software, `Q-SYNTH` reproduces KADE's factor *structure* to a per-factor congruence of at least 0.9995 (loadings within about 0.024, sums of squares within about two percent). Because factor recovery is scale-invariant, this is close enough for the recovery comparison, and it means simulated results speak to what mainstream Q software would produce.

For planning a study, `calculate.mjs` reads the simulation table and returns the smallest design (defining sorts per factor, and the participant count to recruit) that reaches a target recovery rate for a researcher's expected communality, viewpoint overlap, and statement count, showing a conservative and an optimistic recommendation when the inputs fall between simulated conditions. A self-contained browser version (`calculator.html`) offers the same calculation through a G*Power-style form with a recovery-rate plot, and requires no installation. The tool also ships a parameterized sweep (`run.mjs`) that writes recovery rates and confidence intervals to CSV with accompanying figures, parallelized across CPU cores with `worker_threads` and seeded per design cell so that results are bit-identical regardless of the number of workers or the order of execution, a paired comparison of principal-component versus centroid extraction (`robustness.mjs`), and an ecological calibration (`seeding/`) that estimates realistic communality and distinctness from a published dataset (the classic Lipset data distributed with `qmethod` [@Lipset1963; @Brown1980; @Zabala2014]). Recovery is judged by matching each recovered factor array to its planted counterpart with Tucker's coefficient of congruence [@Tucker1951], and factor retention can be studied under an oracle, under parallel analysis [@Horn1965], or against the eigenvalue-over-one rule that the Q literature cautions against, where Brown instead advised preliminarily retaining as many as seven factors before rotation, the non-statistical magic number seven [@Brown1980, pp. 222-223; @Held2016]. A companion methodological paper reports the substantive design findings; this software is the instrument behind them.

# Software design and research impact

The package is small and dependency-light by design (the simulation and the calculator use only the Node.js standard library; a spreadsheet reader is an optional dependency used solely by the KADE cross-check), which keeps it easy to install, test, and audit (its pure-Node core could later run in the browser). Reproducibility is built in rather than promised: the sweep seeds its generator per design cell from the cell's own parameters instead of threading one global stream through the loop, so a result depends only on the parameters and not on execution order or machine. This makes the whole study bit-reproducible and, at the same time, lets it run in parallel across CPU cores through `worker_threads` without changing a single output value, which keeps a full-precision sweep tractable on a laptop. Its intended impact is to give Q researchers a reproducible, quantitative basis for the sample-size and item-count decisions that have until now rested on rules of thumb, and to give reviewers and instructors a concrete way to reason about those decisions. The accompanying methodological paper, which reports the design findings the tool produced, is its first application.

# AI usage disclosure

Generative AI (Claude Opus 4.8, from Anthropic) was used to help write and refactor code, tests, and documentation, and to draft prose, under the authors' direction and review; the authors are responsible for all content. The synthetic data are produced by a deterministic, seeded Monte Carlo simulation and not by any generative model.

# Acknowledgements

We thank the developers of KADE, PQMethod, and `qmethod`, whose openly documented conventions made the cross-validation of `Q-SYNTH` possible.

# References
