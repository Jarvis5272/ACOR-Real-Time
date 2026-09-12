# Input Format

`reads.tsv` must contain exactly these seven tab-separated columns:

```text
dataset_id cluster_id read_id read_sequence original_read_sequence read_quality design_length
```

The first four fields must be non-empty. Clusters must be contiguous and
globally unique; read IDs must be unique within a cluster; design length must
be constant within a cluster. `original_read_sequence` and `read_quality` are
carried as provenance only and never enter reconstruction.

`truth.tsv` is optional. Accepted schemas are `cluster_id,sequence` or the
legacy `row_key,dataset_id,reference_sequence`. Truth is loaded only after a
prediction exists and is used only for evaluation.
