# ACOR v0.2.0 performance and compatibility audit

## Protocol

The v0.1.0 tag and the v0.2.0 release candidate were rebuilt and executed on
the same server, inputs, fixed seed, and default `k=9` configuration. The host
had two Intel Xeon Silver 4316 processors (20 physical cores per socket, two
hardware threads per core), 80 logical CPUs, two NUMA nodes, and
540,671,111,168 bytes of installed memory. The attached A800 accelerators were
not used by ACOR.

The eight datasets were Oligo0, Nature Communications 2019, binned nanopore,
Microsoft CNR, Trellis BMA, DNAFormer nanopore, DNAFormer Illumina, and Zenodo
10943282. Each version was tested at P1, P20, and P40. This is a deterministic
release audit, not a substitute for repeated inferential performance analysis.
Raw measurements are in
[`benchmarks/v0.2/v010_v020_a800.tsv`](../benchmarks/v0.2/v010_v020_a800.tsv).

## Default k=9 compatibility

All 48 runs completed with `PASS`. For every dataset, the v0.1.0 and v0.2.0
prediction files had identical SHA-256 values at P1, P20, and P40. Exact
reconstruction rate, base accuracy, and mean edit distance were therefore also
identical. The new scheduler and state representation do not change the paper
version's default output.

At P1, v0.2.0 was about 4.8% slower at the median because its release path adds
configuration and reusable-state machinery. At P20, v0.2.0 was faster on six
of eight datasets, essentially tied on the two largest datasets, and reduced
engine time by a median factor of 1.17. Representative P20 measurements are:

| Dataset | v0.1.0 (s) | v0.2.0 (s) | old/new |
|---|---:|---:|---:|
| Oligo0 | 0.322 | 0.222 | 1.45x |
| Nature Communications 2019 | 0.332 | 0.239 | 1.39x |
| Binned nanopore | 0.579 | 0.493 | 1.17x |
| Microsoft CNR | 0.589 | 0.507 | 1.16x |
| Trellis BMA | 0.410 | 0.316 | 1.30x |
| DNAFormer nanopore | 5.074 | 5.090 | 1.00x |
| DNAFormer Illumina | 3.431 | 3.403 | 1.01x |
| Zenodo 10943282 | 5.078 | 5.039 | 1.01x |

Seven of eight one-pass P40 measurements also favored v0.2.0. The first Oligo0
P40 v0.2.0 measurement was an isolated timing outlier; a five-order paired
rerun produced mean engine times of 0.377 s for v0.1.0 and 0.218 s for v0.2.0
(1.73x). This rerun is reported separately rather than replacing the raw
one-pass record.

## High-k memory audit

The decisive v0.2.0 change appears at `k>=12`. The former dense representation
allocated counters proportional to `4^k` for every worker. v0.2.0 stores only
observed nodes and edges, so the ledger contribution is
`O(|V_obs| + |E_obs|)`.

The eight datasets were rerun at `k=12` with P1 and P40. All 16 predictions
were byte-identical to their saved dense references. P40 peak RSS was
41.84--658.97 MiB, compared with 22.53--23.14 GiB for the historical dense P20
runs: a 36.0--551.5-fold reduction despite using twice as many workers.

| Dataset | dense k12 P20 | sparse k12 P40 | reduction |
|---|---:|---:|---:|
| Oligo0 | 22.55 GiB | 123.08 MiB | 187.6x |
| Nature Communications 2019 | 22.53 GiB | 50.76 MiB | 454.5x |
| Binned nanopore | 22.56 GiB | 73.87 MiB | 312.7x |
| Microsoft CNR | 22.56 GiB | 83.04 MiB | 278.2x |
| Trellis BMA | 22.53 GiB | 41.84 MiB | 551.5x |
| DNAFormer nanopore | 22.99 GiB | 540.50 MiB | 43.6x |
| DNAFormer Illumina | 23.14 GiB | 658.97 MiB | 36.0x |
| Zenodo 10943282 | 23.00 GiB | 538.29 MiB | 43.8x |

Across these k12 runs, aggregate compute time decreased from 293.744 s at P1
to 8.568 s at P40 (34.28x). Per-dataset speedup ranged from 25.53x to 38.87x;
the median was 35.03x, or 87.6% median parallel efficiency. The complete
32-cell formula-k/k12 matrix is in
[`benchmarks/v0.2/high_k_sparse_matrix.tsv`](../benchmarks/v0.2/high_k_sparse_matrix.tsv).

## Interpretation

v0.2.0 is not a new reconstruction objective. It is an output-compatible
execution and memory release: default short-oligo quality is unchanged,
parallel throughput improves most clearly for many-cluster workloads, and the
structural high-k memory failure is removed. Optimal worker count remains
workload-dependent; very small jobs may saturate before all physical cores are
used.
