# Input format

## Assigned reads

`reads.tsv` contains exactly seven tab-separated columns in this order:

```text
dataset_id	cluster_id	read_id	read_sequence	original_read_sequence	read_quality	design_length
```

The first four fields must be non-empty. Each file contains one non-empty
dataset; each cluster occupies one contiguous block and has a unique
`cluster_id`. Read IDs are unique within their cluster.

`design_length` is constant within a cluster and ranges from the selected
`k` to 1,000,000 bases. The optional `--length` argument overrides this
field for reconstruction.

Only `read_sequence` enters the reconstruction model.
`original_read_sequence` and `read_quality` retain input provenance.

## Reference sequences

The optional `truth.tsv` accepts either of these tab-separated schemas:

```text
cluster_id	sequence
```

```text
row_key	dataset_id	reference_sequence
```

References are loaded after a prediction exists and are used for evaluation.
The `--aim` option requires `truth.tsv`; an existing reference file is
also evaluated with `--no-aim`.
