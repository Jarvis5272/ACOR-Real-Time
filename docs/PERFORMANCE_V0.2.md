# ACOR v0.2.0 performance and compatibility

## Comparison protocol

The v0.1.0 tag and v0.2.0 release candidate were rebuilt on the same server
and compared using identical inputs, a fixed seed, and the default `k=9`
configuration. The host had two Intel Xeon Silver 4316 processors, with
20 physical cores per socket and two hardware threads per core. It provided
80 logical CPUs, two NUMA nodes, and 540,671,111,168 bytes of reported memory.
ACOR ran on the CPUs.

The eight datasets were Oligo0, ncomms19 365-dishes, binned nanopore,
Microsoft CNR, Trellis BMA, DNAFormer nanopore, DNAFormer Illumina, and
Zenodo 10943282. Each version was tested at P1, P20, and P40, where P denotes
worker count. Each configuration has one recorded timing observation.

The `engine_time_s` field reports runner wall time, including input loading
and prediction writing. The `compute_time_s` field reports reconstruction
compute time separately. Measurements are available in
[`benchmarks/v0.2/v010_v020_a800.tsv`](../benchmarks/v0.2/v010_v020_a800.tsv).

## Default k=9 results

All 48 runs completed with `PASS`. The 24 matched version/dataset/worker
pairs had identical prediction-file SHA-256 values, Exact, Accuracy, and
Mean ED. These comparisons establish output compatibility on the tested
inputs and configurations.

At P1, v0.2.0 runner times were higher in all eight datasets, with a median
increase of about 4.8%. At P20, the median old/new runner-time ratio was
1.17. The three largest workloads had ratios close to 1.00.
The P20 measurements are:

| Dataset | v0.1.0 (s) | v0.2.0 (s) | old/new |
|---|---:|---:|---:|
| Oligo0 | 0.322 | 0.222 | 1.45x |
| ncomms19 365-dishes | 0.332 | 0.239 | 1.39x |
| Binned nanopore | 0.579 | 0.493 | 1.17x |
| Microsoft CNR | 0.589 | 0.507 | 1.16x |
| Trellis BMA | 0.410 | 0.316 | 1.30x |
| DNAFormer nanopore | 5.074 | 5.090 | 1.00x |
| DNAFormer Illumina | 3.431 | 3.403 | 1.01x |
| Zenodo 10943282 | 5.078 | 5.039 | 1.01x |

Seven of eight P40 runner measurements favored v0.2.0. The recorded Oligo0
P40 row favored v0.1.0. A separate five-order paired comparison reported
mean runner times of 0.377 s and 0.218 s, respectively, giving a 1.73x ratio.
That comparison remains separate from the one-observation matrix.

Compute times were higher for v0.2.0 at every tested default-k worker
configuration. Runner and compute measurements therefore describe different
aspects of version performance.

## High-k storage and scaling

For `k>=12`, v0.2.0 replaces dense counter allocation with storage for
observed contexts and transitions. The graph-ledger storage contribution is
`O(|V_obs| + |E_obs|)`. Input reads, length statistics, raw-read frequencies,
and worker state also contribute to process memory.

The eight datasets were evaluated at `k=12` with P1 and P40. All
16 prediction files were byte-identical to their saved dense references.
P40 peak RSS was 41.84–658.97 MiB, compared with 22.53–23.14 GiB for
historical dense P20 runs. The recorded configuration ratios ranged from
36.0 to 551.5. The worker counts differ between these storage comparisons.

| Dataset | dense k12 P20 | sparse k12 P40 | reduction |
|---|---:|---:|---:|
| Oligo0 | 22.55 GiB | 123.08 MiB | 187.6x |
| ncomms19 365-dishes | 22.53 GiB | 50.76 MiB | 454.5x |
| Binned nanopore | 22.56 GiB | 73.87 MiB | 312.7x |
| Microsoft CNR | 22.56 GiB | 83.04 MiB | 278.2x |
| Trellis BMA | 22.53 GiB | 41.84 MiB | 551.5x |
| DNAFormer nanopore | 22.99 GiB | 540.50 MiB | 43.6x |
| DNAFormer Illumina | 23.14 GiB | 658.97 MiB | 36.0x |
| Zenodo 10943282 | 23.00 GiB | 538.29 MiB | 43.8x |

Summed compute times across these k12 runs were 293.744 s at P1 and
8.567 s at P40, giving a 34.29x ratio. Per-dataset ratios ranged from
25.53x to 38.87x, with a median of 35.03x and median efficiency of 87.6%.

The complete 32-configuration matrix includes k12 and dataset-specific
context settings in
[`benchmarks/v0.2/high_k_sparse_matrix.tsv`](../benchmarks/v0.2/high_k_sparse_matrix.tsv).

## Practical implications

Version 0.2.0 preserves the tested default predictions and reduces the
memory reserved for unobserved high-k graph keys. Worker choice depends
on workload size and cluster distribution. Small jobs can reach their
best recorded time before all physical cores are used. The measurements
describe the reported builds, inputs, and execution settings.
