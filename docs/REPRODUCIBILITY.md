# Reproducibility

## Read order

An explicit seed derives a 128-bit per-cluster entropy value from
`sha256(acor|minimal|seed|dataset|cluster_id)[:16]`. Python's
`random.Random(entropy).shuffle(permutation)` produces the read order.

Without `--seeds`, the CLI obtains a 32-hex-character seed from the
operating system and records it before reconstruction. A later invocation
with `--seeds` regenerates the order from that value.

Worker configurations within the same session share one order per
seed/round. The session and engine manifests record the configuration
and execution details needed to identify the run.

## Build profiles

With GCC or Clang, `portable-release` uses `-O3 -DNDEBUG`, while
`frozen-compatible` uses the reference flags `-O2 -DNDEBUG`.
MSVC uses `/O2 /DNDEBUG` for both profiles. Both profiles select the
Release configuration and require identical fixture predictions.

CMake disables interprocedural optimization and rejects architecture-specific
or LTO flags in both profiles, including inherited compile and link flags.
The `frozen-compatible` name identifies a compiler profile; the rebuilt
source is evaluated through prediction compatibility.

See the [regression results](REAL10_REGRESSION.md) for the recorded seed,
worker configurations, prediction comparisons, and supporting files.
