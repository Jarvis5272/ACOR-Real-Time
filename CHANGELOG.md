# Changelog

## 0.2.0 - 2026-09-22

- Preserve the v0.1.0 `NO_STOP` reconstruction semantics and default `k=9`
  output while replacing per-cluster worker setup with a persistent,
  NUMA-aware adaptive worker pool.
- Reuse worker-local reconstruction state and deterministically split only
  genuinely large clusters, improving parallel throughput without changing
  prediction identity.
- Use an observed-key sparse evidence ledger automatically for `k>=12`,
  avoiding the exponential dense allocation that previously made high-k and
  longer-sequence runs impractical.
- Add the public `--k 5..15` option, exact ledger-merge regression coverage,
  and auditable configuration/backend fields in run manifests.

## 0.1.0 - 2026-08-22

- Initial source release of the ACOR reconstruction backbone.
- Frozen `NO_STOP` reconstruction semantics and reproducible per-cluster shuffling.
- Transactional Python orchestration, strict input/evaluator closure, and auditable manifests.
- CMake frozen-compatible and portable-release build profiles.
