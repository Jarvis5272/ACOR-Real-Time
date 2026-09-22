# ACOR Reconstruction Backbone

ACOR reconstructs a DNA sequence from a cluster of noisy reads. Version 0.2.0
keeps the `NO_STOP` reconstruction semantics of v0.1.0 while adding an
adaptive persistent worker pool and a sparse high-k evidence ledger. Every
read is consumed and no learned early-stop model is included.

## Scope

- Release mode: `NO_STOP`
- Execution model: `CAUSAL_OFFLINE_STREAM_REPLAY`
- Default configuration: `kmc_k9_b16_lognormal`
- Supported k-mer sizes: `5` through `15`
- Execution backend: persistent NUMA-aware adaptive worker pool
- License: MIT
- Author: Jarvis5272

The CLI scans a complete static dataset, generates a reproducible random read
permutation for every cluster, and feeds reads to the incremental model in that
order. This is a high-fidelity offline streaming simulation, not a persistent
interface connected directly to a sequencing instrument.

## Requirements

- Linux x86_64 with CPU-affinity support
- Python 3.10 or newer
- CMake 3.16 or newer
- A C++17 compiler
- zlib and pthreads

## Build

Portable release build:

```bash
./scripts/build.sh portable-release
```

Frozen-compatible build using the original server flags:

```bash
./scripts/build.sh frozen-compatible
```

Equivalent explicit commands are:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DACOR_BUILD_PROFILE=portable-release
cmake --build build --parallel
cmake --install build --prefix .
```

## Test

```bash
./scripts/test.sh
```

The tests include deterministic order golden files, CLI modes, random-seed
replay by regenerating the order, thread argument propagation, truth and hidden
column isolation, failure cleanup, evaluator boundaries, and comparator
strict-weak-order/frontier identity. A C++ regression additionally verifies
that independently accumulated evidence ledgers merge to the serial result.

## Run

```bash
python3 acor.py run \
  --data /path/to/dataset_directory \
  [--aim | --no-aim] \
  [--length N] \
  [--k 5..15] \
  [--threads 1,4,16,20,32] \
  [--seeds seed1,seed2 | --rounds 3] \
  --results /path/to/results
```

`reads.tsv` is required. If `truth.tsv` exists it is evaluated automatically,
including with `--no-aim`; `--aim` requires truth. Without `--seeds` or
`--rounds`, the operating system supplies one random 32-hex-character seed.
The manifest records it before runner launch so it can be replayed with
`--seeds`; replay regenerates the order and never reuses the first order file.

Each successfully committed session contains `RUN_MANIFEST.json`,
`RUN_STATUS.json`, `SUMMARY.tsv`, one shared order per seed/round, and one
transactionally published directory per thread/round. Failed or interrupted
sessions are explicitly named and never expose a formal round prediction.

The terminal `time` value is `reconstruction_engine_wall`. Round and session
wall times remain in summaries/manifests. Because early stopping is absent,
saving/trigger fields are omitted or `NOT_APPLICABLE`.

The default `--k 9` is backward-compatible with v0.1.0. For `k>=12`, ACOR
automatically stores only observed graph counters instead of allocating the
full direct-indexed counter space. The selected k-mer size, effective
configuration, and execution backend are recorded in `RUN_MANIFEST.json` and
each engine `COMMAND.json`.

## Data

No research dataset is included in this repository. See
[`docs/INPUT_FORMAT.md`](docs/INPUT_FORMAT.md) for the seven-column reads
contract and optional truth schemas. The small `tests/fixtures/` files are
test-only examples.

## Reproducibility and limits

- [`docs/ALGORITHM.md`](docs/ALGORITHM.md)
- [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md)
- [`docs/LIMITATIONS.md`](docs/LIMITATIONS.md)
- [`docs/REAL10_REGRESSION.md`](docs/REAL10_REGRESSION.md)
- [`docs/PERFORMANCE_V0.2.md`](docs/PERFORMANCE_V0.2.md)

The evaluator is byte- and case-sensitive. It does not canonicalize sequences.
The algorithm canonicalization rules apply only inside reconstruction and are
frozen in `provenance/ALGORITHM_LOCK.json`.
