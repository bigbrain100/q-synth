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

## Release checklist (maintainers)

1. All tests pass (`npm test`) and the calculator smoke-test runs.
2. Bump `version` in `package.json` and `CITATION.cff` (keep `date-released` current).
3. Create the GitHub release. Zenodo archives it automatically via the webhook.
4. **Fix the Zenodo record by hand — every release.** Zenodo's GitHub integration
   currently ignores `.zenodo.json` and `CITATION.cff`, so each auto-archived
   version is created with the GitHub account name as the only author and the
   repository slug as the title. On zenodo.org, open the new version → Edit →
   set the title to the CITATION.cff title and the creators to the three authors
   (Lee, Ming-Shinn; Wang, Wan-Hsiang; Tsai, Yi-Jung — with ORCIDs and
   affiliations, in that order) → Publish. The DOI does not change.
5. Verify with the API that the record is correct:
   `curl -s https://zenodo.org/api/records/<new-version-id>` — check `title`
   and `creators`. The all-versions (concept) DOI stays 10.5281/zenodo.21429739.

## Code of conduct

Participation in this project is governed by the [Code of Conduct](CODE_OF_CONDUCT.md).
