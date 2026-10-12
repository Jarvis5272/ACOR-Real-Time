# ACOR

ACOR reconstructs DNA-storage sequences by accumulating evidence from noisy
reads and decoding each cluster after its final read. Version 0.2.0 adds a
persistent worker pool and sparse storage for larger k-mer contexts. The
default reconstruction configuration remains compatible with v0.1.0 in the
reported comparisons.

## Scope

- Release mode: `NO_STOP`, consuming every read before one final decode
- Execution model: `CAUSAL_OFFLINE_STREAM_REPLAY`
- Default configuration: `kmc_k9_b16_lognormal`
- Supported k-mer sizes: `5` through `15`
- Execution backend: persistent NUMA-aware adaptive worker pool
- License: MIT
- Author: Jarvis5272

The CLI reads a complete dataset of assigned clusters and replays a
reproducible random permutation of each cluster's reads. The reconstruction
model updates its evidence incrementally during replay. Live sequencing
connections and online read assignment require integration outside this CLI.

## Requirements

- Linux x86_64 with CPU-affinity support
- Python 3.10 or newer
- CMake 3.16 or newer
- A C++17 compiler
- zlib and pthreads

## Build

Build the default portable profile:

```bash
./scripts/build.sh portable-release
```

For the reference compiler profile:

```bash
./scripts/build.sh frozen-compatible
```

With GCC or Clang, these profiles use `-O3 -DNDEBUG` and `-O2 -DNDEBUG`,
respectively. Both profiles disable LTO and reject architecture-specific
compiler flags. See [build reproducibility](docs/REPRODUCIBILITY.md).

The portable build can also be configured directly:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DACOR_BUILD_PROFILE=portable-release
cmake --build build --parallel
cmake --install build --prefix .
```

## Test

```bash
./scripts/test.sh
```

Tests cover read-order reproducibility, CLI options, thread counts, reference
isolation, interrupted runs, evaluation boundaries, and score ordering.
A C++ regression compares merged block evidence with serial accumulation.

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

`reads.tsv` supplies the assigned reads and design lengths. If `truth.tsv`
exists, ACOR evaluates the predictions automatically, including with
`--no-aim`; `--aim` requires that file. References are used only for
evaluation after reconstruction. See [input formats](docs/INPUT_FORMAT.md).

Use `--seeds` to specify reproducible read orders or `--rounds` to generate
multiple random rounds. Without either option, one random 32-hex-character
seed is generated and recorded before reconstruction. Reusing that seed
regenerates the order in a new session.

For `k>=12`, ACOR stores observed graph counters in a sparse ledger.
The selected k-mer size, configuration, and execution backend are recorded in
`RUN_MANIFEST.json` and each engine `COMMAND.json`.

## Results

A completed session contains `RUN_MANIFEST.json`, `RUN_STATUS.json`,
`SUMMARY.tsv`, one shared read order per seed/round, and a result directory
for each thread/round. Round outputs are published after successful
completion. Failed or interrupted sessions carry `.failed` or
`.interrupted` suffixes and retain their status records.

The terminal `time` value reports reconstruction engine wall time
(`reconstruction_engine_wall`). Summaries also record round and session
wall times. Saving and trigger fields are omitted or `NOT_APPLICABLE` in
`NO_STOP` mode. See [output formats](docs/OUTPUT_FORMAT.md).

## Data and documentation

Research datasets are obtained separately; the repository includes small
test examples in `tests/fixtures/`. Input columns and optional reference
schemas are described in [input formats](docs/INPUT_FORMAT.md).

- [Algorithm](docs/ALGORITHM.md)
- [Reproducibility](docs/REPRODUCIBILITY.md)
- [Requirements and limitations](docs/LIMITATIONS.md)
- [v0.1.0 regression results](docs/REAL10_REGRESSION.md)
- [v0.2.0 performance and compatibility](docs/PERFORMANCE_V0.2.md)

Reconstruction maps canonical nucleotide symbols to uppercase and other
symbols to `N`. Evaluation compares sequences byte-for-byte, including
case. The reconstruction settings are recorded in
`provenance/ALGORITHM_LOCK.json`.
