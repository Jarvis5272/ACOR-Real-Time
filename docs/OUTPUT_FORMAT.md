# Output format

## Session and round files

Each completed session contains `RUN_MANIFEST.json`, `RUN_STATUS.json`,
and `SUMMARY.tsv`. A result directory is created for each thread/round,
containing `pred.tsv`, `summary.tsv`, and, when truth is available,
`evaluation.tsv`. Order files are shared under the session's `.orders/`;
runner artifacts are stored within each round's `.internal/`.

Round outputs are published after successful completion. Failed or
interrupted sessions use `.failed` or `.interrupted` suffixes, with
the corresponding state recorded in `RUN_STATUS.json`.

## Predictions and timing

`pred.tsv` contains `row_key`, `dataset_id`, `cluster_id`,
`legacy_prefix_alias`, `stop_index`, `consumed_reads`, `total_reads`,
`reconstructed_sequence`, and `status`.
In `NO_STOP` mode, `stop_index`, `consumed_reads`, and `total_reads`
are equal.

The terminal `time` value is reconstruction engine wall time.
Summary and manifest files also retain round and session wall times.
Saving and trigger fields are omitted or `NOT_APPLICABLE` in this release.
