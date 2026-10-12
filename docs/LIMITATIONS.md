# Requirements and limitations

## Input and execution

ACOR processes assigned clusters with known design lengths. The CLI uses
`CAUSAL_OFFLINE_STREAM_REPLAY` over a complete static dataset. Live
sequencing connections and online clustering require separate integration.

In `NO_STOP` mode, all reads are consumed before the final decode.
This release has no learned early-stop module; saving and trigger metrics
are unavailable.

The user selects `k`; automatic selection from read length, coverage,
and error rate is outside the current release.

## Memory and reconstruction

Sparse storage allocates graph counters for observed contexts and
transitions. Peak memory also depends on read storage, length statistics,
distinct raw-read strings, and concurrent worker states. New observed keys
can increase memory use as more reads arrive.

Beam pruning can discard the correct sequence path during search.
Reconstruction quality varies with the chosen context length and supplied
reads.

## Platform and evaluation

The supported platform is Linux x86_64 with C++17, zlib, pthreads, and
CPU-affinity support. Python and CMake requirements are listed in the
[README](../README.md).

The evaluator compares sequences byte-for-byte, including case.
Nucleotide canonicalization applies within reconstruction.

Both build profiles forbid architecture-specific and LTO flags.
This includes `-march`, `-mtune`, and `-flto` with GCC/Clang, plus
`/arch` and `/GL` with MSVC. See [reproducibility](REPRODUCIBILITY.md)
for the compiler profiles.
