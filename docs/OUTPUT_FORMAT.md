# Output Format

Each committed round contains `pred.tsv`, `summary.tsv`, and, when truth is
available, `evaluation.tsv`. Internal order and runner artifacts live below
`.internal/`.

`pred.tsv` records `cluster_id`, `reconstructed_sequence`, real `stop_index`,
`consumed_reads`, `total_reads`, and the explicitly named
`legacy_prefix_alias`. In `NO_STOP`, stop/consumed/total are equal.

The terminal `time` value is reconstruction engine wall time. Summary and
manifest files also retain round and session wall times.
