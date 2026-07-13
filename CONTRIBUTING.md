# Contributing to Q-SYNTH

Contributions, bug reports, and feature requests are welcome.

## Reporting bugs or asking questions

Please open an issue on the project's issue tracker. For a bug, include:

- what you ran (the exact command and options),
- what you expected and what happened,
- your Node.js version (`node --version`) and operating system.

## Suggesting features

Open an issue describing the use case. Because Q-SYNTH aims to give trustworthy
design guidance, proposals that add or change the simulation or the recovery
criterion should explain the methodological rationale and, where possible, cite
supporting literature.

## Contributing code

1. Fork the repository and create a branch for your change.
2. Keep the code dependency-free: the simulation, calculator, and tests use only
   the Node.js standard library. The only external package is `xlsx`, which is not
   a project dependency; install it manually (`npm install xlsx`) only to run the
   optional KADE cross-check.
3. Add or update tests in `tests/` (run with `npm test`, i.e. `node --test`).
   New numerical behavior should come with a test that pins it down.
4. If you change extraction, rotation, or the recovery metric, re-run the
   validation in `validation/` against a KADE export where relevant, and note
   the result in your pull request.
5. Open a pull request describing the change and its motivation.

## Running the checks locally

```bash
npm test            # unit tests (node --test)
node calculate.mjs --c 0.6 --rho 0.2 --M 40   # smoke-test the calculator
```

## Code of conduct

Participation in this project is governed by the [Code of Conduct](CODE_OF_CONDUCT.md).
