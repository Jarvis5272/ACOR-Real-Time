# Limitations

- Execution is causal offline stream replay over a complete static dataset,
  not a persistent live sequencing interface.
- This release intentionally contains no learned early-stop module and reports
  no saving or trigger claim.
- Linux, C++17, zlib, pthreads, and CPU affinity support are required.
- The evaluator is byte- and case-sensitive and does not canonicalize input.
- `-march=native` is available only in the frozen-compatible profile and is not
  used by the portable release.
