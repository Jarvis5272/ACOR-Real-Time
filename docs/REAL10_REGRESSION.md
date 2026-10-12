# REAL10 regression results

The v0.1.0 comparison used ten real sources at P1/P4/P16/P20/P32 with seed
`6c2f91a47bd803e52a9c14f8d730be65`, giving 100 runs across the reference
and release packages. Each source also had a random-seed P20 run and a
fresh replay, giving 20 additional runs.

The comparison checked missing, extra, duplicate, and mismatched sequences,
together with Exact, Accuracy, and Mean ED. Prediction-schema changes were
recorded separately from sequence differences.

## Results

All 120 tasks completed:

- All 50 fixed-seed reference/release pairs had identical reconstructed
  sequences, with no missing, extra, duplicate, or mismatched outputs.
- Exact, Accuracy, and Mean ED were identical in every pair.
- Each source produced one canonical prediction across P1/P4/P16/P20/P32.
- `runner_workers` and the recorded worker CPU list matched the requested
  count in all 50 release runs.
- All 10 random-seed P20 runs reproduced when their read orders were
  regenerated from the recorded seeds.

## Evaluation scope

The Zenodo source contained 680 reference rows for clusters without input
reads. Both packages were evaluated against the same reference view,
excluding these rows. The original and effective reference hashes are
recorded in the regression files.

The reference evaluator omitted integer count fields from its summary.
For its 50 rows, `exact_count`, `total_ed`, and `total_truth_bases`
were taken from the paired evaluator after canonical prediction identity
was established. These fields carry the provenance label
`PAIRED_NEW_EVALUATOR_AFTER_CANONICAL_IDENTITY`.

## Supporting files

The records are included under `regression/`:

- `REAL10_DATA_MANIFEST.tsv`
- `REAL10_RUN_LEDGER.tsv`
- `REAL10_REGRESSION_LEDGER.tsv`
- `EXPECTED_PREDICTION_HASHES.tsv`
- `RANDOM_REPLAY_LEDGER.tsv`
- `THREAD_IDENTITY_AUDIT.tsv`
- `FINAL_STATUS.json`
