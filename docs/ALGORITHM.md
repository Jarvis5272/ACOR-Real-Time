# ACOR Reconstruction Backbone

ACOR consumes each cluster in a reproducible shuffled read order and updates a
direct-indexed k-mer graph after every arrival. The frozen release uses `k=9`,
beam width `16`, the lognormal length prior, and `FinalOnlyExact` decoding.

The final candidate priority is frozen as: exact design-length beam, exact
length raw-read fallback, `M±1...±5` beam, arm0 consensus, then the
highest-frequency raw read. Ties use the frozen lexicographic path/sequence
fields. The release mode is `NO_STOP`; every read is consumed.

This package implements causal offline stream replay. It scans a complete
static dataset, derives one per-cluster permutation, then feeds reads to the
incremental model in that order. It is not a direct sequencing-instrument API.
