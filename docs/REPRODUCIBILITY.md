# Reproducibility

An explicit seed deterministically derives a 128-bit per-cluster entropy value
from `acor|minimal|seed|dataset|cluster_id`. Python's
`random.Random(entropy).shuffle(permutation)` then produces the read order.

Without `--seeds`, the CLI obtains a 32-hex-character seed from the operating
system. The manifest records that seed before runner launch. Replaying it with
`--seeds` regenerates the order; the first order file is not reused.

The frozen-compatible build uses the original server flags. The portable
release omits `-march=native` and LTO. Both profiles are required to produce
identical fixture predictions.
