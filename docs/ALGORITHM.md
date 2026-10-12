# ACOR algorithm

ACOR reconstructs one sequence from a non-empty cluster of assigned reads
with a supplied design length `L`. Its default configuration uses `k=9`,
beam width `16`, a lognormal length prior, and `FinalOnlyExact` decoding.
The release mode is `NO_STOP`; every read contributes before decoding once.

## Evidence accumulation

Each read is normalized to uppercase A, C, G, and T, with other symbols
mapped to `N`. Valid k-mer contexts and observed single-base transitions
update the cluster ledger. Windows containing noncanonical symbols are
skipped. The ledger also records the first and last valid contexts, read
length statistics, and normalized read frequencies.

The public CLI accepts `k=5..15`. Contexts of size `k>=12` use observed-key
sparse counters; smaller contexts use dense counters. Both representations
store the same integer evidence and use the same final decoding rules.

## Closure-time decoding

Starting contexts are ranked by start support, with overall context counts
as a fallback. A beam search follows observed transitions and retains
candidates at design-window and empirical-window lengths. The maximum search
length `M` is the larger of `L+5` and the empirical upper length bound.

The first non-empty candidate family supplies the output:

1. Beam candidates of exactly `L` bases.
2. Normalized raw reads of exactly `L` bases.
3. Beam candidates at `L-d` and `L+d`, tested jointly for increasing
   `d=1..5`. Lengths below `k` are excluded.
4. Candidates from the empirical length window, named `arm0` in the code.
5. The highest-frequency normalized raw read, regardless of its length.

Beam pruning uses cumulative path scores. Complete candidates use normalized
path support, the length prior, and terminal support. Finite scores are
ordered through integer rounding at scale `10^12`, followed by lexical
tie-breaking. Raw reads are ordered by frequency and then lexical sequence.

## Execution

Independent clusters run through a persistent pinned worker pool, reusing
worker-local state. Large clusters can be split into evidence blocks and
merged in a fixed order before one final decode. Reported regression checks
compare predictions across worker counts and execution versions.

The CLI uses `CAUSAL_OFFLINE_STREAM_REPLAY`: it loads a complete static
dataset and feeds each cluster's reads in a reproducible shuffled order.
Live sequencing and upstream read assignment are described in the
[requirements and limitations](LIMITATIONS.md).
