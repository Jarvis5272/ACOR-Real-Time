# REAL10 Regression

The release gate runs the frozen package and this source package on ten real
sources at P1/P4/P16/P20/P32 with seed
`6c2f91a47bd803e52a9c14f8d730be65` (100 runs). Each source also performs one
OS-random P20 run and a fresh `--seeds` replay of the recorded seed (20 runs).

The hard gate is zero missing, extra, duplicate, or mismatched reconstructed
sequences and zero difference in exact, accuracy, and mean edit distance.
Schema-only raw-pred changes are audited separately.

## Release result

Version 0.1.0 closed all 120 actual tasks:

- 50/50 frozen-old versus release-new fixed-seed pairs passed;
- 0 missing, extra, duplicate, or mismatched reconstructed sequences;
- exact, accuracy, and mean edit distance were identical in every pair;
- each source produced one canonical prediction across P1/P4/P16/P20/P32;
- `runner_workers` and the recorded worker CPU list matched the requested
  thread count in all 50 release runs;
- 10/10 random-seed P20 runs reproduced after regenerating, rather than
  copying, the order from the recorded seed.

The Zenodo source contains 680 truth rows for clusters with no input reads.
The regression harness excludes only those truth-only rows and records both
the original and effective truth SHA. The old and new package are evaluated
against the same input-closed truth view; no source reads or original truth
files are modified. All other sources required no truth-scope adjustment.

The frozen evaluator did not expose integer count fields in the old package
summary. `exact_count`, `total_ed`, and `total_truth_bases` for those 50 rows
are copied from the paired strict evaluator only after same-order canonical
prediction identity is established. Their provenance is explicitly
`PAIRED_NEW_EVALUATOR_AFTER_CANONICAL_IDENTITY`.

Machine-readable evidence is included under `regression/`:

- `REAL10_DATA_MANIFEST.tsv`
- `REAL10_RUN_LEDGER.tsv`
- `REAL10_REGRESSION_LEDGER.tsv`
- `EXPECTED_PREDICTION_HASHES.tsv`
- `RANDOM_REPLAY_LEDGER.tsv`
- `THREAD_IDENTITY_AUDIT.tsv`
- `FINAL_STATUS.json`
