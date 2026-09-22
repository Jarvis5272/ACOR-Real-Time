# Limitations

- Execution is causal offline stream replay over a complete static dataset,
  not a persistent live sequencing interface.
- This release intentionally contains no learned early-stop module and reports
  no saving or trigger claim.
- ACOR expects reads to be assigned to clusters before reconstruction; online
  clustering is outside this release.
- `k` is user-selected. The software records it but does not yet infer an
  optimal value from read length, coverage, and error rate.
- Sparse storage removes the dense high-k allocation, but peak memory still
  grows with the number of distinct observed nodes/edges and concurrent worker
  states.
- Linux, C++17, zlib, pthreads, and CPU affinity support are required.
- The evaluator is byte- and case-sensitive and does not canonicalize input.
- `-march=native` is available only in the frozen-compatible profile and is not
  used by the portable release.
