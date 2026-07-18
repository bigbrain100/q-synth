# Q-SYNTH calculator tutorial

*A step-by-step guide to the [Q-SYNTH design calculator](https://bigbrain100.github.io/q-synth/calculator.html), written for Q researchers with no factor-analysis background. 繁體中文版請見 [TUTORIAL.zh-Hant.md](TUTORIAL.zh-Hant.md).*

## What this tool is

Quantitative researchers use G\*Power for a priori planning: give it an effect size and the power you want, and it returns the sample size you need. Q methodology never had an equivalent. Decisions about how many participants to recruit and how many statements to write have rested on rules of thumb — "40–80 statements", "4–5 defining sorts per factor" — that were never systematically tested.

Q-SYNTH tests them. It plants known viewpoints in synthetic data, generates thousands of simulated Q sorts under controlled conditions, runs a standard Q factor analysis, and measures how often the planted viewpoints are recovered. The calculator reads the resulting lookup table: you describe the study you are planning, and it returns how many defining sorts per factor — and hence how many participants — that design needs.

The correspondence with G\*Power:

| G\*Power | Q-SYNTH calculator |
|---|---|
| Effect size | Expected communality *c* (how strong the viewpoint signal is) |
| Power (1 − β) | Target recovery (how sure you want to be) |
| Output: required *N* | Output: defining sorts per factor → total participants |

The calculator runs entirely in your browser — no installation, no registration, works on a phone:

**<https://bigbrain100.github.io/q-synth/calculator.html>**

## The five inputs, one at a time

### 1. Expected communality (*c*) — how firmly people hold their view

This is the input that stops most people, and the idea is plain: what proportion of a person's sort is explained by the viewpoint they share? People with sharp, settled opinions produce sorts almost fully determined by their viewpoint — high *c*. People whose sorting wavers produce low *c*.

**Where to get it.** Find a published Q study on a nearby topic (or run a small pilot) and open its factor-loading table — KADE, PQMethod, and the `qmethod` R package all print one:

1. Keep only the rows flagged as defining sorts.
2. Square each loading.
3. Average the squares. That average is *c*.

Worked once: the defining sorts on a factor load .75, .80, .70. Squared: .56, .64, .49. Mean ≈ 0.56 → enter **0.6**.

**No reference study?** Enter 0.6 — the typical value in Q studies. Reserve 0.8 for participants who are seasoned stakeholders with firm, articulate positions.

### 2. Viewpoint overlap (*ρ*) — how alike the viewpoints are

Two viewpoints that agree on most statements and differ on a few are heavily overlapping — hard to separate statistically, so more sorts are needed. Viewpoints that genuinely diverge are easy to separate.

**Where to get it.** In the same reference study, find the small "correlations between factors" table (all three programs report it). Average the absolute values. Example: correlations .25, .31, .18 → mean ≈ 0.25 → enter 0.2, or 0.4 to be safe.

**No reference study?** Enter 0.4. Opposing camps on contested topics usually agree on more statements than you expect.

### 3. Viewpoints (*K*) — how many voices you expect

How many shared viewpoints do you expect the topic to yield? Two to four is typical in Q. If unsure, enter the largest number you plausibly expect — over-preparing is cheap, under-recruiting is not.

(If planning around an expected factor count sounds contrary to Q's philosophy, see the FAQ on the calculator page — short answer: Q already makes this assumption implicitly whenever it plans "4–5 sorts per factor"; the calculator just makes it explicit and testable.)

### 4. Statements (*M*) — the size of your Q set

How many statements participants will sort. This one is your design choice, drawn from your concourse; the calculator tells you whether it is enough. The simulations are unambiguous on one point: 40+ statements are meaningfully better than 20, because every statement is one more observation steadying each person-to-person correlation.

### 5. Target recovery — how sure you want to be

The analogue of statistical power: "if this study were rerun a hundred times, in how many should *all* the viewpoints come back cleanly?" The default 0.90 (nine in ten) is a sound choice; keep it.

Leave the **Advanced options** (congruence threshold, distribution, retention rule, non-defining sorts) at their defaults — they exist for methodologists.

## Reading the result

Suppose the calculator returns **8 defining sorts per factor** and you expect 3 viewpoints:

- 8 × 3 = 24 defining sorts.
- Not every participant loads cleanly on a factor (about 80% do, in practice).
- 24 ÷ 0.8 = 30 → **recruit about 30 participants.**

If your *c* and *ρ* fall between simulated conditions, you get a **conservative** and an **optimistic** recommendation plus a recovery-rate curve. Plan toward the conservative one unless you are confident your viewpoints are strong and distinct.

If the calculator warns that even 12 sorts per factor cannot reach your target, do not just recruit harder — that design region is fragile. Sharpen the Q set (raise *c*), plan for fewer viewpoints, or add statements.

The fastest way to learn the tool: press **Load example**, then **Calculate**, and trace the numbers above.

## Reporting it in your methods section

> Sample size was determined a priori with the Q-SYNTH recovery-simulation calculator (Lee, Tsai, & Wang, 2026). Assuming communality 0.6, viewpoint overlap 0.4, three expected viewpoints, and 40 statements, 8 defining sorts per factor were required for a .90 recovery rate, so 30 participants were recruited.

When a reviewer asks where your *N* came from, this is a far stronger answer than "similar studies used about 30."

## Citing Q-SYNTH

> Lee, M.-S., Tsai, Y.-J., & Wang, W.-H. (2026). *Q-SYNTH: A recovery-simulation tool for sample-size and item-count guidance in Q methodology* [Computer software]. https://github.com/bigbrain100/q-synth

Or use GitHub's **Cite this repository** button (top right of the repo page).
