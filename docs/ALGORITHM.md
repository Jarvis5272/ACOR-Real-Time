# ACOR Reconstruction Backbone

ACOR consumes each cluster in a reproducible shuffled read order and updates a
k-mer graph after every arrival. The default remains `k=9`, beam width `16`,
the lognormal length prior, and `FinalOnlyExact` decoding. The public CLI
accepts `k=5..15`; configurations with `k>=12` use an observed-key sparse
ledger while preserving the same integer evidence and final decoding rules.

The final candidate priority is frozen as: exact design-length beam, exact
length raw-read fallback, `M±1...±5` beam, arm0 consensus, then the
highest-frequency raw read. Ties use the frozen lexicographic path/sequence
fields. The release mode is `NO_STOP`; every read is consumed.

Independent clusters are dispatched through a persistent pinned worker pool.
Worker-local state is reused across clusters. A cluster is split into
independent evidence blocks only when it is large enough to amortize the merge;
block ledgers are combined in a fixed order before one final decode. Scheduling
therefore changes resource use, not candidate scores or output ordering.

This package implements causal offline stream replay. It scans a complete
static dataset, derives one per-cluster permutation, then feeds reads to the
incremental model in that order. It is not a direct sequencing-instrument API.
