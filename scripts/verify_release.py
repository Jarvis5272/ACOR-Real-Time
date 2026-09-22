#!/usr/bin/env python3
"""Fail-closed ACOR source and binary release verifier.

The verifier never extracts an archive.  Its allowlists are anchored in the
checked-out worktree, and every unregistered worktree or archive member is an
error.  It intentionally has no third-party Python dependency.
"""

from __future__ import annotations

import sys

# A release audit must not import `json.py`, `tarfile.py`, etc. from the
# attacker-controlled publication tree before it gets a chance to inspect that
# tree.  `sys` is built in, so this guard runs before every shadowable import.
# All supported command-line invocations therefore use `python3 -I`.
if __name__ == "__main__" and not sys.flags.isolated:
    print("RELEASE_AUDIT=FAIL reason=invoke verifier with python3 -I", file=sys.stderr)
    raise SystemExit(2)

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import stat
import tarfile
import tempfile
import unicodedata
import zlib
from pathlib import Path, PurePosixPath
from typing import Dict, Iterable, List, Mapping, Sequence, Set, Tuple


VERSION = "0.2.0"
TOP = f"acor-reconstruction-backbone-{VERSION}"
RELEASE_MTIME = 1_704_067_200
SOURCE_MANIFEST = "scripts/release_source_members.txt"
BINARY_MANIFEST = "scripts/release_binary_members.txt"

REQUIRED_SOURCE: Set[str] = {
    ".github/workflows/ci.yml",
    ".gitignore",
    "CHANGELOG.md",
    "CITATION.cff",
    "CMakeLists.txt",
    "LICENSE",
    "README.md",
    "THIRD_PARTY_NOTICES.md",
    "VERSION",
    "acor.py",
    "benchmarks/v0.2/high_k_sparse_matrix.tsv",
    "benchmarks/v0.2/v010_v020_a800.tsv",
    "docs/ALGORITHM.md",
    "docs/INPUT_FORMAT.md",
    "docs/LIMITATIONS.md",
    "docs/OUTPUT_FORMAT.md",
    "docs/PERFORMANCE_V0.2.md",
    "docs/REAL10_REGRESSION.md",
    "docs/REPRODUCIBILITY.md",
    "provenance/ALGORITHM_LOCK.json",
    "provenance/FROZEN_REFERENCE.json",
    "provenance/PATCH_LEDGER.tsv",
    "provenance/SOURCE_PROVENANCE.tsv",
    "regression/EXPECTED_PREDICTION_HASHES.tsv",
    "regression/FINAL_STATUS.json",
    "regression/RANDOM_REPLAY_LEDGER.tsv",
    "regression/REAL10_DATA_MANIFEST.tsv",
    "regression/REAL10_REGRESSION_LEDGER.tsv",
    "regression/REAL10_RUN_LEDGER.tsv",
    "regression/THREAD_IDENTITY_AUDIT.tsv",
    "scripts/build.sh",
    "scripts/make_release.sh",
    "scripts/regression_real10.py",
    SOURCE_MANIFEST,
    BINARY_MANIFEST,
    "scripts/test.sh",
    "scripts/verify_release.py",
    "src/acor_kernel.cpp",
    "src/acor_kernel.hpp",
    "src/acor_runner.cpp",
    "src/ed_pairs.cpp",
    "tests/common.py",
    "tests/fixtures/tiny_aim/reads.tsv",
    "tests/fixtures/tiny_aim/truth.tsv",
    "tests/fixtures/tiny_noaim/reads.tsv",
    "tests/fixtures/tiny_tie/frontier_reads.tsv",
    "tests/fixtures/tiny_tie/frontier_tiny_ordered.tsv",
    "tests/golden/frontier_tie.tsv",
    "tests/golden/frontier_tie_final.tsv",
    "tests/golden/frontier_tiny.tsv",
    "tests/golden/frontier_tiny_final.tsv",
    "tests/golden/tiny_fixed_order.jsonl.gz",
    "tests/golden/tiny_fixed_pred.tsv",
    "tests/test_cli.py",
    "tests/test_comparator.py",
    "tests/test_cpp_boundaries.py",
    "tests/test_evaluator.py",
    "tests/test_failure_cleanup.py",
    "tests/test_ledger_reduction.cpp",
    "tests/test_seed_replay.py",
    "tests/test_thread_identity.py",
    "tests/test_truth_isolation.py",
    "tests/tools/comparator_properties.cpp",
    "tests/tools/frontier_probe.cpp",
}

REQUIRED_BINARY: Set[str] = {
    "CHANGELOG.md",
    "CITATION.cff",
    "LICENSE",
    "README.md",
    "THIRD_PARTY_NOTICES.md",
    "VERSION",
    "acor.py",
    "benchmarks/v0.2/high_k_sparse_matrix.tsv",
    "benchmarks/v0.2/v010_v020_a800.tsv",
    "bin/acor_runner",
    "bin/ed_pairs",
    "docs/ALGORITHM.md",
    "docs/INPUT_FORMAT.md",
    "docs/LIMITATIONS.md",
    "docs/OUTPUT_FORMAT.md",
    "docs/PERFORMANCE_V0.2.md",
    "docs/REAL10_REGRESSION.md",
    "docs/REPRODUCIBILITY.md",
    "provenance/ALGORITHM_LOCK.json",
    "provenance/FROZEN_REFERENCE.json",
    "provenance/PATCH_LEDGER.tsv",
    "provenance/SOURCE_PROVENANCE.tsv",
    "regression/EXPECTED_PREDICTION_HASHES.tsv",
    "regression/FINAL_STATUS.json",
    "regression/RANDOM_REPLAY_LEDGER.tsv",
    "regression/REAL10_DATA_MANIFEST.tsv",
    "regression/REAL10_REGRESSION_LEDGER.tsv",
    "regression/REAL10_RUN_LEDGER.tsv",
    "regression/THREAD_IDENTITY_AUDIT.tsv",
    "tests/fixtures/tiny_aim/reads.tsv",
    "tests/fixtures/tiny_aim/truth.tsv",
    "tests/fixtures/tiny_noaim/reads.tsv",
    "tests/fixtures/tiny_tie/frontier_reads.tsv",
    "tests/fixtures/tiny_tie/frontier_tiny_ordered.tsv",
}

FORBIDDEN_SOURCE_COMPONENTS = {
    "__pycache__",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
    "backup",
    "backups",
    "bin",
    "build",
    "cache",
    "data",
    "logs",
    "results",
}
FORBIDDEN_MODEL_SUFFIXES = {
    ".joblib",
    ".model",
    ".onnx",
    ".pickle",
    ".pkl",
    ".pt",
    ".pth",
}
SOURCE_CODE_SUFFIXES = {
    ".c",
    ".cc",
    ".cmake",
    ".cpp",
    ".cxx",
    ".h",
    ".hpp",
    ".py",
    ".sh",
    ".yaml",
    ".yml",
}
MAX_TEXT_SCAN = 16 * 1024 * 1024
MAX_ARCHIVE_MEMBER = 128 * 1024 * 1024
MAX_ARCHIVE_TOTAL = 256 * 1024 * 1024
MAX_COMPRESSED_ARCHIVE = 256 * 1024 * 1024
MAX_WORKTREE_FILE = 128 * 1024 * 1024
MAX_EXPANDED_ARCHIVE = MAX_ARCHIVE_TOTAL + 16 * 1024 * 1024

# These are independent release-controller anchors.  They are deliberately
# compiled into the verifier instead of being learned from editable package
# provenance.  A package author therefore cannot change an algorithm file and
# its in-tree ledger together and still pass this release gate.
APPROVED_ALGORITHM_SHA256: Mapping[str, str] = {
    "acor.py": "366d9e80ff625b95e85d929813eb01a49b772ef7783fb752968329ee3ac87e13",
    "src/acor_runner.cpp": "5daa71f1181395ff2856fcc0d721bccc938e658e3bb609b9faa68c1e030d398e",
    "src/acor_kernel.cpp": "878d11bb899e4a370394482432400292a56bbc170df9c1282312b48407eae15b",
    "src/acor_kernel.hpp": "7fc1fbed15f032c8c3526356cb3dddc960976250c68c071709c159493bc785c1",
    "src/ed_pairs.cpp": "a191b6743448a3902f9ad1cd951fab2464ff3ad2f41f1f105f18f0753acf85a1",
}

EXPECTED_ALGORITHM_LOCK: Mapping[str, object] = {
    "beam_width": 16,
    "canonical": "A/a->A,C/c->C,G/g->G,T/t->T,other->N",
    "config": "kmc_k9_b16_lognormal (default)",
    "execution_model": "CAUSAL_OFFLINE_STREAM_REPLAY",
    "execution_backend": "LEDGER_ADAPTIVE_POOL",
    "final_decode": "FinalOnlyExact",
    "k_default": 9,
    "k_supported": [5, 15],
    "length_prior": "lognormal",
    "mode": "NO_STOP",
    "sparse_ledger_threshold": 12,
    "seed_derivation": "sha256(acor|minimal|seed|dataset|cluster_id)[:16]",
    "support_bonus": 0.15,
    "version": "0.2.0",
}

EXPECTED_FROZEN_REFERENCE: Mapping[str, object] = {
    "acor_py_sha256": "3b1f00eeab20d7b3ea65bab6c3880229f67eee02df554a8077316203c3de9c17",
    "ed_pairs_binary_sha256": "a9ccb311ac47f34a3934bb868bc3a21af0f77de0cf1355d5b680cd9bef516874",
    "ed_pairs_source_sha256": "9c9e60bf16738b57458125c0ea16a41c71943bb9937e938d4d26c0cd30c33c3f",
    "fixed_regression_seed": "6c2f91a47bd803e52a9c14f8d730be65",
    "kernel_source_sha256": "debde33a5cab10fdde4ec63eb16d4bbe8399ca65ea5f2fce9b7003968eb2384e",
    "runner_binary_sha256": "23611c1d783167ec16913c8b151a1e4895d9bd7ce4f4a1843bc0a9947be260bf",
    "runner_source_sha256": "9acc4dc286859fa786202935c7d5affe53092e132af03dc70e53f1b9a055b04e",
}

EXPECTED_LICENSE = """MIT License

Copyright (c) 2026 Jarvis5272

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
""".encode("utf-8")

EXPECTED_CITATION = """cff-version: 1.2.0
message: "If you use ACOR Reconstruction Backbone, please cite this software."
title: "ACOR Reconstruction Backbone"
version: 0.2.0
date-released: 2026-09-22
authors:
  - family-names: "Jarvis5272"
license: MIT
repository-code: "https://github.com/Jarvis5272/ACOR-Real-Time"
""".encode("utf-8")

EXPECTED_DATASETS: Set[str] = {
    "oligo0",
    "ncomms19_365dishes",
    "binned_nanopore",
    "microsoft_cnr",
    "trellisbma",
    "dnaformer_nanopore",
    "dnaformer_illumina",
    "zenodo10943282",
    "err1816980_fountain",
    "err1816980_official",
}
EXPECTED_THREADS = {1, 4, 16, 20, 32}
FIXED_REGRESSION_SEED = "6c2f91a47bd803e52a9c14f8d730be65"
FROZEN_RUNNER_SHA256 = "23611c1d783167ec16913c8b151a1e4895d9bd7ce4f4a1843bc0a9947be260bf"
REGRESSION_RELEASE_RUNNER_SHA256 = "66ce485be663e65967a2e17dee1a6b23b84789e405ffef8bb97fe80ea74a1b2e"

SECRET_PATTERNS: Sequence[Tuple[str, re.Pattern[str]]] = (
    ("private key", re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----")),
    ("GitHub token", re.compile(r"\b(?:ghp|gho|ghu|ghs|ghr)_[A-Za-z0-9]{20,}\b")),
    ("GitHub fine-grained token", re.compile(r"\bgithub_pat_[A-Za-z0-9_]{20,}\b")),
    ("AWS access key", re.compile(r"\b(?:AKIA|ASIA)[A-Z0-9]{16}\b")),
    ("OpenAI-style secret", re.compile(r"\bsk-[A-Za-z0-9_-]{20,}\b")),
)

PRIVATE_PATH_PATTERNS: Sequence[Tuple[str, re.Pattern[str]]] = (
    ("Unix account path", re.compile(r"/(?:home|Users|root)/[A-Za-z0-9._-]+(?:/|\b)")),
    ("Windows absolute path", re.compile(r"(?i)\b[A-Z]:[\\/][^\s'\"<>]+")),
    (
        "project server alias",
        re.compile(
            r"(?i)\b(?:"
            + re.escape("my" + "server")
            + "|"
            + re.escape("116-" + "PR4768GW")
            + "|"
            + re.escape("han" + "linxuan")
            + r")\b"
        ),
    ),
    ("private IPv4 address", re.compile(r"\b(?:10\.(?:\d{1,3}\.){2}\d{1,3}|192\.168\.(?:\d{1,3}\.)\d{1,3}|172\.(?:1[6-9]|2\d|3[01])\.(?:\d{1,3}\.)\d{1,3})\b")),
)

MODEL_CODE_PATTERNS: Sequence[Tuple[str, re.Pattern[str]]] = (
    ("joblib dependency", re.compile(r"(?m)^\s*(?:from\s+joblib|import\s+joblib)\b")),
    ("scikit-learn dependency", re.compile(r"(?m)^\s*(?:from\s+sklearn|import\s+sklearn)\b")),
    ("learned prediction call", re.compile(r"\bpredict_proba\s*\(")),
    (
        "historical early-stop model",
        re.compile(
            r"\b(?:"
            + re.escape("MODEL_" + "v4")
            + "|"
            + re.escape("S6_" + "ML")
            + "|"
            + re.escape("HistGradient" + "Boosting")
            + r")\b"
        ),
    ),
)


class VerificationError(RuntimeError):
    pass


def fail(message: str) -> None:
    raise VerificationError(message)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def safe_relpath(value: str) -> str:
    if not value or value != value.strip():
        fail(f"empty or whitespace-padded member path: {value!r}")
    if "\\" in value or "\x00" in value or "//" in value:
        fail(f"non-canonical member path: {value!r}")
    if any(ord(character) < 32 or ord(character) == 127 for character in value):
        fail(f"control character in member path: {value!r}")
    if unicodedata.normalize("NFC", value) != value:
        fail(f"non-NFC member path: {value!r}")
    path = PurePosixPath(value)
    if path.is_absolute() or value.startswith("/"):
        fail(f"absolute member path: {value!r}")
    if any(part in {"", ".", ".."} for part in path.parts):
        fail(f"path traversal or non-canonical component: {value!r}")
    if re.match(r"^[A-Za-z]:", value):
        fail(f"drive-qualified member path: {value!r}")
    canonical = path.as_posix()
    if canonical != value:
        fail(f"member path is not canonical: {value!r}")
    return canonical


def _read_regular_file_secure(path: Path, limit: int, label: str) -> bytes:
    """Read a bounded regular file without following a Unix parent symlink.

    On POSIX, every path component is opened relative to the already-opened
    parent directory (`openat` via ``dir_fd``).  A concurrent directory/symlink
    swap therefore cannot redirect the final open.  Windows lacks this Python
    interface, so its fallback retains the lstat/fstat identity checks used by
    local developer tests; release creation itself is Linux-only.
    """

    absolute = Path(os.path.abspath(path))
    if os.name == "posix" and os.open in os.supports_dir_fd:
        directory_flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
        directory_flags |= getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_CLOEXEC", 0)
        file_flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_CLOEXEC", 0)
        descriptors: List[int] = []
        try:
            parent_fd = os.open(absolute.anchor, directory_flags)
            descriptors.append(parent_fd)
            for component in absolute.parts[1:-1]:
                next_fd = os.open(component, directory_flags, dir_fd=parent_fd)
                descriptors.append(next_fd)
                parent_fd = next_fd
            descriptor = os.open(absolute.parts[-1], file_flags, dir_fd=parent_fd)
            descriptors.append(descriptor)
        except OSError as exc:
            for opened_fd in reversed(descriptors):
                os.close(opened_fd)
            fail(f"cannot securely open {label}: {exc}")
    else:
        reject_symlink_chain(absolute, label)
        try:
            before = os.lstat(absolute)
        except OSError as exc:
            fail(f"cannot lstat {label}: {exc}")
        if stat.S_ISLNK(before.st_mode):
            fail(f"symlink forbidden for {label}: {absolute}")
        file_flags = os.O_RDONLY | getattr(os, "O_BINARY", 0) | getattr(os, "O_NOFOLLOW", 0)
        try:
            descriptor = os.open(absolute, file_flags)
        except OSError as exc:
            fail(f"cannot securely open {label}: {exc}")
        descriptors = [descriptor]

    try:
        opened = os.fstat(descriptor)
        if not stat.S_ISREG(opened.st_mode):
            fail(f"{label} is not a regular file")
        if opened.st_size > limit:
            fail(f"{label} exceeds safety limit ({opened.st_size} > {limit})")
        chunks: List[bytes] = []
        remaining = limit + 1
        while remaining:
            chunk = os.read(descriptor, min(1024 * 1024, remaining))
            if not chunk:
                break
            chunks.append(chunk)
            remaining -= len(chunk)
        data = b"".join(chunks)
        if len(data) > limit:
            fail(f"{label} exceeds safety limit while reading")
        after = os.fstat(descriptor)
        if (
            opened.st_dev,
            opened.st_ino,
            opened.st_size,
            getattr(opened, "st_mtime_ns", int(opened.st_mtime * 1_000_000_000)),
        ) != (
            after.st_dev,
            after.st_ino,
            after.st_size,
            getattr(after, "st_mtime_ns", int(after.st_mtime * 1_000_000_000)),
        ):
            fail(f"{label} changed while being read")
        return data
    finally:
        for opened_fd in reversed(descriptors):
            os.close(opened_fd)


def read_root_member(root: Path, rel: str) -> bytes:
    safe_relpath(rel)
    return _read_regular_file_secure(root / rel, MAX_WORKTREE_FILE, f"worktree member {rel}")


def reject_symlink_chain(path: Path, label: str) -> None:
    absolute = path.absolute()
    current = Path(absolute.anchor)
    for component in absolute.parts[1:]:
        current = current / component
        try:
            metadata = os.lstat(current)
        except FileNotFoundError:
            # The caller's strict resolve/open reports the missing path with the
            # same fail-closed top-level error handling.
            return
        except OSError as exc:
            fail(f"cannot lstat {label} path component {current}: {exc}")
        if stat.S_ISLNK(metadata.st_mode):
            fail(f"symlink component forbidden in {label} path: {current}")


def read_manifest(root: Path, rel: str) -> Tuple[List[str], bytes]:
    raw = read_root_member(root, rel)
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        fail(f"manifest is not UTF-8: {rel}: {exc}")
    if "\r" in text:
        fail(f"manifest must use LF line endings: {rel}")
    if "# UNFINALIZED_ALLOWLIST" in text:
        fail(f"manifest remains unfinalized: {rel}")
    members: List[str] = []
    seen: Set[str] = set()
    seen_casefold: Dict[str, str] = {}
    for line_no, line in enumerate(text.splitlines(), 1):
        if not line:
            continue
        if line.startswith("#"):
            continue
        if line != line.strip():
            fail(f"whitespace-padded manifest entry {rel}:{line_no}")
        value = safe_relpath(line)
        if value in seen:
            fail(f"duplicate manifest entry {rel}:{line_no}: {value}")
        folded = value.casefold()
        if folded in seen_casefold and seen_casefold[folded] != value:
            fail(
                f"case-fold-colliding manifest entries: "
                f"{seen_casefold[folded]!r} and {value!r}"
            )
        seen_casefold[folded] = value
        seen.add(value)
        members.append(value)
    if members != sorted(members):
        fail(f"manifest entries are not bytewise sorted: {rel}")
    return members, raw


def parent_directories(files: Iterable[str]) -> Set[str]:
    result: Set[str] = set()
    for value in files:
        path = PurePosixPath(value)
        for parent in path.parents:
            if parent == PurePosixPath("."):
                break
            result.add(parent.as_posix())
    return result


def expected_mode(rel: str, kind: str) -> int:
    if kind == "binary" and rel in {"acor.py", "bin/acor_runner", "bin/ed_pairs"}:
        return 0o755
    if kind == "source" and (
        rel == "acor.py"
        or (rel.startswith("scripts/") and (rel.endswith(".sh") or rel == "scripts/verify_release.py"))
        or (rel.startswith("regression/") and rel.endswith(".py"))
    ):
        return 0o755
    return 0o644


def scan_text(data: bytes, rel: str, source_code: bool) -> None:
    if len(data) > MAX_TEXT_SCAN:
        fail(f"text file exceeds scan limit: {rel} ({len(data)} bytes)")
    if b"\x00" in data:
        fail(f"unexpected binary content in text member: {rel}")
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        fail(f"non-UTF-8 non-binary file: {rel}")
    for label, pattern in SECRET_PATTERNS:
        if pattern.search(text):
            fail(f"{label} detected in {rel}")
    for label, pattern in PRIVATE_PATH_PATTERNS:
        if pattern.search(text):
            fail(f"{label} detected in {rel}")
    if source_code:
        for label, pattern in MODEL_CODE_PATTERNS:
            if pattern.search(text):
                fail(f"{label} detected in source code: {rel}")


def scan_source_path(rel: str, data: bytes) -> None:
    path = PurePosixPath(rel)
    lower_parts = {part.lower() for part in path.parts}
    forbidden = sorted(lower_parts & FORBIDDEN_SOURCE_COMPONENTS)
    if forbidden:
        fail(f"forbidden source path component {forbidden[0]!r}: {rel}")
    if path.suffix.lower() in FORBIDDEN_MODEL_SUFFIXES:
        fail(f"model or pickle asset forbidden in source release: {rel}")
    if data.startswith(b"\x7fELF"):
        fail(f"ELF binary forbidden in source release: {rel}")
    if path.suffix.lower() == ".gz":
        if rel != "tests/golden/tiny_fixed_order.jsonl.gz" or not data.startswith(b"\x1f\x8b"):
            fail(f"unregistered or malformed compressed source fixture: {rel}")
        try:
            with gzip.GzipFile(fileobj=io.BytesIO(data), mode="rb") as stream:
                payload = stream.read(MAX_TEXT_SCAN + 1)
        except (EOFError, OSError) as exc:
            fail(f"malformed compressed source fixture {rel}: {exc}")
        if len(payload) > MAX_TEXT_SCAN:
            fail(f"decompressed source fixture exceeds scan limit: {rel}")
        scan_text(payload, rel + "::<gzip>", True)
        return
    code_like = path.suffix.lower() in SOURCE_CODE_SUFFIXES or path.name == "CMakeLists.txt"
    scan_text(data, rel, code_like)


def scan_elf_payload(data: bytes, rel: str) -> None:
    """Reject private build paths and excluded learned-model residue in ELF."""
    if len(data) < 64:
        fail(f"implausibly short ELF executable: {rel}")
    if data[0:4] != b"\x7fELF" or data[4] != 2 or data[5] != 1 or data[6] != 1:
        fail(f"ELF must be ELF64 little-endian version 1: {rel}")
    elf_type = int.from_bytes(data[16:18], "little")
    machine = int.from_bytes(data[18:20], "little")
    elf_version = int.from_bytes(data[20:24], "little")
    if elf_type not in {2, 3} or machine != 62 or elf_version != 1:
        fail(
            f"ELF is not an x86_64 ET_EXEC/ET_DYN executable: {rel} "
            f"type={elf_type} machine={machine} version={elf_version}"
        )
    forbidden_fragments = {
        b"/home/": "Unix account path",
        (b"C:" + b"\\" + b"Users" + b"\\"): "Windows account path",
        ("my" + "server").encode("ascii"): "project server alias",
        ("116-" + "PR4768GW").encode("ascii"): "project server hostname",
        ("han" + "linxuan").encode("ascii"): "private account name",
        ("MODEL_" + "v4").encode("ascii"): "historical early-stop model",
        ("S6_" + "ML").encode("ascii"): "historical early-stop model",
        b"predict_proba": "learned prediction call",
    }
    lowered = data.lower()
    for fragment, label in forbidden_fragments.items():
        probe = fragment.lower()
        if probe in lowered:
            fail(f"{label} detected in ELF payload: {rel}")


def _decode_utf8(payloads: Mapping[str, bytes], rel: str) -> str:
    try:
        return payloads[rel].decode("utf-8")
    except (KeyError, UnicodeDecodeError) as exc:
        fail(f"cannot decode required semantic payload {rel}: {exc}")


def _reject_duplicate_json_keys(pairs: Sequence[Tuple[str, object]]) -> Dict[str, object]:
    result: Dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            fail(f"duplicate JSON key in semantic release contract: {key}")
        result[key] = value
    return result


def _load_semantic_json(payloads: Mapping[str, bytes], rel: str) -> Dict[str, object]:
    text = _decode_utf8(payloads, rel)
    try:
        value = json.loads(
            text,
            object_pairs_hook=_reject_duplicate_json_keys,
            parse_constant=lambda token: fail(f"non-finite JSON constant in {rel}: {token}"),
        )
    except json.JSONDecodeError as exc:
        fail(f"invalid JSON in {rel}: {exc}")
    if not isinstance(value, dict):
        fail(f"semantic JSON contract must have object root: {rel}")
    return value


def _canonical_json(value: object) -> str:
    try:
        return json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        fail(f"non-canonical semantic JSON value: {exc}")


def _parse_exact_tsv(
    payloads: Mapping[str, bytes], rel: str, expected_fields: Sequence[str]
) -> List[Dict[str, str]]:
    text = _decode_utf8(payloads, rel)
    if "\r" in text:
        fail(f"semantic TSV must use LF line endings: {rel}")
    try:
        reader = csv.DictReader(io.StringIO(text), delimiter="\t")
        if reader.fieldnames != list(expected_fields):
            fail(f"unexpected TSV schema in {rel}")
        rows = list(reader)
    except csv.Error as exc:
        fail(f"invalid TSV in {rel}: {exc}")
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        fail(f"ragged TSV row in {rel}")
    return rows


REAL10_FIELDS = """dataset	threads	seed	old_order_sha256	new_order_sha256	reads_sha256	truth_sha256	truth_original_sha256	truth_effective_sha256	truth_only_excluded	old_runner_sha256	new_runner_sha256	old_runner_workers	new_runner_workers	old_runner_worker_cpus	new_runner_worker_cpus	old_raw_pred_sha256	new_raw_pred_sha256	old_canonical_pred_sha256	new_canonical_pred_sha256	missing	extra	duplicate	mismatch_clusters	mismatch_ids	old_exact	new_exact	old_accuracy	new_accuracy	old_mean_ed	new_mean_ed	old_engine_time_s	new_engine_time_s	old_round_time_s	new_round_time_s	protocol_old	protocol_new	schema_only_raw_difference	status	old_exact_count	new_exact_count	old_total_ed	new_total_ed	old_total_truth_bases	new_total_truth_bases	old_metric_counts_source	new_metric_counts_source""".split("\t")
RANDOM_FIELDS = """dataset	threads	seed	random_order_sha256	replay_order_sha256	random_order_inode	replay_order_inode	order_regenerated	random_canonical_pred_sha256	replay_canonical_pred_sha256	random_exact	replay_exact	status""".split("\t")


def verify_regression_ledgers(payloads: Mapping[str, bytes]) -> None:
    fixed_rows = _parse_exact_tsv(
        payloads, "regression/REAL10_REGRESSION_LEDGER.tsv", REAL10_FIELDS
    )
    if len(fixed_rows) != 50:
        fail(f"REAL10 regression ledger must contain exactly 50 fixed pairs, got {len(fixed_rows)}")
    expected_pairs = {(dataset, threads) for dataset in EXPECTED_DATASETS for threads in EXPECTED_THREADS}
    observed_pairs: Set[Tuple[str, int]] = set()
    sha_pattern = re.compile(r"^[0-9a-f]{64}$")
    fixed_sha_fields = (
        "old_order_sha256",
        "new_order_sha256",
        "reads_sha256",
        "truth_sha256",
        "truth_original_sha256",
        "truth_effective_sha256",
        "old_runner_sha256",
        "new_runner_sha256",
        "old_raw_pred_sha256",
        "new_raw_pred_sha256",
        "old_canonical_pred_sha256",
        "new_canonical_pred_sha256",
    )
    equality_pairs = (
        ("old_order_sha256", "new_order_sha256"),
        ("old_canonical_pred_sha256", "new_canonical_pred_sha256"),
        ("old_exact", "new_exact"),
        ("old_accuracy", "new_accuracy"),
        ("old_mean_ed", "new_mean_ed"),
        ("old_exact_count", "new_exact_count"),
        ("old_total_ed", "new_total_ed"),
        ("old_total_truth_bases", "new_total_truth_bases"),
    )
    for row in fixed_rows:
        try:
            threads = int(row["threads"])
        except ValueError:
            fail("non-integer threads in REAL10 regression ledger")
        pair = (row["dataset"], threads)
        if pair in observed_pairs:
            fail(f"duplicate fixed regression pair: {pair}")
        observed_pairs.add(pair)
        if row["seed"] != FIXED_REGRESSION_SEED:
            fail(f"unexpected fixed regression seed: {pair}")
        if row["status"] != "PASS":
            fail(f"non-PASS fixed regression row: {pair}")
        for field in ("missing", "extra", "duplicate", "mismatch_clusters"):
            if row[field] != "0":
                fail(f"non-zero {field} in fixed regression row: {pair}")
        if row["mismatch_ids"] != "NONE":
            fail(f"mismatch ids present in fixed regression row: {pair}")
        if row["old_runner_sha256"] != FROZEN_RUNNER_SHA256:
            fail(f"unexpected old runner SHA in fixed regression row: {pair}")
        if row["new_runner_sha256"] != REGRESSION_RELEASE_RUNNER_SHA256:
            fail(f"unexpected release runner SHA in fixed regression row: {pair}")
        if any(not sha_pattern.fullmatch(row[field]) for field in fixed_sha_fields):
            fail(f"invalid SHA-256 field in fixed regression row: {pair}")
        if row["old_runner_workers"] != row["new_runner_workers"] or row["old_runner_workers"] != str(threads):
            fail(f"worker-count mismatch in fixed regression row: {pair}")
        expected_cpus = ",".join(str(value) for value in range(threads))
        if row["old_runner_worker_cpus"] != row["new_runner_worker_cpus"] or row["old_runner_worker_cpus"] != expected_cpus:
            fail(f"worker CPU mapping mismatch in fixed regression row: {pair}")
        if row["schema_only_raw_difference"] != "1":
            fail(f"raw prediction difference is not schema-only: {pair}")
        if row["old_raw_pred_sha256"] == row["new_raw_pred_sha256"]:
            fail(f"schema-only raw prediction evidence is not actually distinct: {pair}")
        if row["protocol_old"] != "LEGACY_HARDCODED_EXCLUDED" or row["protocol_new"] != "NOT_APPLICABLE":
            fail(f"unexpected protocol ledger value: {pair}")
        if (
            row["old_metric_counts_source"]
            != "PAIRED_NEW_EVALUATOR_AFTER_CANONICAL_IDENTITY"
            or row["new_metric_counts_source"] != "NATIVE_EVALUATOR"
        ):
            fail(f"unexpected metric-count provenance in fixed regression row: {pair}")
        expected_truth_only = "680" if row["dataset"] == "zenodo10943282" else "0"
        if row["truth_only_excluded"] != expected_truth_only:
            fail(f"truth-only exclusion mismatch in fixed regression row: {pair}")
        for old_field, new_field in equality_pairs:
            if row[old_field] != row[new_field]:
                fail(f"fixed regression zero-delta violation {old_field}/{new_field}: {pair}")
        for field in ("old_exact", "new_exact", "old_accuracy", "new_accuracy"):
            try:
                value = float(row[field])
            except ValueError:
                fail(f"invalid quality value {field} in fixed regression row: {pair}")
            if not math.isfinite(value) or not 0.0 <= value <= 1.0:
                fail(f"out-of-range quality value {field} in fixed regression row: {pair}")
        for field in ("old_mean_ed", "new_mean_ed"):
            try:
                value = float(row[field])
            except ValueError:
                fail(f"invalid ED value {field} in fixed regression row: {pair}")
            if not math.isfinite(value) or value < 0.0:
                fail(f"out-of-range ED value {field} in fixed regression row: {pair}")
        for field in ("old_engine_time_s", "new_engine_time_s", "old_round_time_s", "new_round_time_s"):
            try:
                value = float(row[field])
            except ValueError:
                fail(f"invalid timing value {field} in fixed regression row: {pair}")
            if not math.isfinite(value) or value <= 0.0:
                fail(f"non-positive timing value {field} in fixed regression row: {pair}")
        for field in ("old_exact_count", "new_exact_count", "old_total_ed", "new_total_ed"):
            try:
                value = int(row[field])
            except ValueError:
                fail(f"invalid count {field} in fixed regression row: {pair}")
            if value < 0:
                fail(f"negative count {field} in fixed regression row: {pair}")
        for field in ("old_total_truth_bases", "new_total_truth_bases"):
            try:
                value = int(row[field])
            except ValueError:
                fail(f"invalid truth-base count {field} in fixed regression row: {pair}")
            if value <= 0:
                fail(f"non-positive truth-base count {field} in fixed regression row: {pair}")
    if observed_pairs != expected_pairs:
        fail("REAL10 fixed regression does not cover the exact 10-source x 5-thread matrix")

    random_rows = _parse_exact_tsv(
        payloads, "regression/RANDOM_REPLAY_LEDGER.tsv", RANDOM_FIELDS
    )
    if len(random_rows) != 10:
        fail(f"random replay ledger must contain exactly 10 pairs, got {len(random_rows)}")
    observed_datasets: Set[str] = set()
    observed_seeds: Set[str] = set()
    seed_pattern = re.compile(r"^[0-9a-f]{32}$")
    for row in random_rows:
        dataset = row["dataset"]
        if dataset in observed_datasets:
            fail(f"duplicate random replay dataset: {dataset}")
        observed_datasets.add(dataset)
        if row["threads"] != "20" or row["status"] != "PASS" or row["order_regenerated"] != "1":
            fail(f"invalid random replay status/threads/regeneration: {dataset}")
        seed = row["seed"]
        if not seed_pattern.fullmatch(seed) or seed in observed_seeds:
            fail(f"invalid or duplicate random replay seed: {dataset}")
        observed_seeds.add(seed)
        if row["random_order_sha256"] != row["replay_order_sha256"] or not sha_pattern.fullmatch(row["random_order_sha256"]):
            fail(f"random order replay mismatch: {dataset}")
        if row["random_canonical_pred_sha256"] != row["replay_canonical_pred_sha256"] or not sha_pattern.fullmatch(row["random_canonical_pred_sha256"]):
            fail(f"random prediction replay mismatch: {dataset}")
        if row["random_exact"] != row["replay_exact"]:
            fail(f"random replay exact mismatch: {dataset}")
        try:
            exact_value = float(row["random_exact"])
        except ValueError:
            fail(f"invalid random replay exact: {dataset}")
        if not math.isfinite(exact_value) or not 0.0 <= exact_value <= 1.0:
            fail(f"out-of-range random replay exact: {dataset}")
        try:
            random_inode = int(row["random_order_inode"])
            replay_inode = int(row["replay_order_inode"])
            if random_inode <= 0 or replay_inode <= 0:
                fail(f"non-positive random replay inode: {dataset}")
            if random_inode == replay_inode:
                fail(f"random replay order was not independently regenerated: {dataset}")
        except ValueError:
            fail(f"invalid random replay inode: {dataset}")
    if observed_datasets != EXPECTED_DATASETS:
        fail("random replay ledger does not cover the exact 10-source set")


def verify_semantic_contract(payloads: Mapping[str, bytes]) -> None:
    lock = _load_semantic_json(payloads, "provenance/ALGORITHM_LOCK.json")
    frozen = _load_semantic_json(payloads, "provenance/FROZEN_REFERENCE.json")
    final_status = _load_semantic_json(payloads, "regression/FINAL_STATUS.json")
    readme_text = _decode_utf8(payloads, "README.md")

    if _canonical_json(lock) != _canonical_json(EXPECTED_ALGORITHM_LOCK):
        fail("ALGORITHM_LOCK.json differs from the exact approved NO_STOP replay contract")
    if _canonical_json(frozen) != _canonical_json(EXPECTED_FROZEN_REFERENCE):
        fail("FROZEN_REFERENCE.json differs from the exact approved reference object")
    for token in ("NO_STOP", "CAUSAL_OFFLINE_STREAM_REPLAY"):
        if token not in readme_text:
            fail(f"README.md misses public scope token: {token}")

    for rel, approved_sha in APPROVED_ALGORITHM_SHA256.items():
        if rel in payloads and sha256_bytes(payloads[rel]) != approved_sha:
            fail(f"algorithm payload differs from independent frozen SHA-256: {rel}")

    if payloads.get("LICENSE") != EXPECTED_LICENSE:
        fail("LICENSE is not the exact MIT License grant for Copyright (c) 2026 Jarvis5272")
    if payloads.get("CITATION.cff") != EXPECTED_CITATION:
        fail("CITATION.cff is not the exact MIT/Jarvis5272 release citation")

    provenance_fields = ("component", "frozen_source_sha256", "release_source", "relation")
    provenance_rows = _parse_exact_tsv(
        payloads, "provenance/SOURCE_PROVENANCE.tsv", provenance_fields
    )
    expected_provenance = [
        {
            "component": "python_cli",
            "frozen_source_sha256": "3b1f00eeab20d7b3ea65bab6c3880229f67eee02df554a8077316203c3de9c17",
            "release_source": "acor.py",
            "relation": "v0.1-compatible defaults plus public k selection and backend audit",
        },
        {
            "component": "runner",
            "frozen_source_sha256": "9acc4dc286859fa786202935c7d5affe53092e132af03dc70e53f1b9a055b04e",
            "release_source": "src/acor_runner.cpp",
            "relation": "v0.1 algorithm plus persistent adaptive NUMA-aware scheduling",
        },
        {
            "component": "kernel",
            "frozen_source_sha256": "debde33a5cab10fdde4ec63eb16d4bbe8399ca65ea5f2fce9b7003968eb2384e",
            "release_source": "src/acor_kernel.cpp",
            "relation": "v0.1 algorithm plus exact sparse high-k evidence storage and merge",
        },
        {
            "component": "evaluator",
            "frozen_source_sha256": "9c9e60bf16738b57458125c0ea16a41c71943bb9937e938d4d26c0cd30c33c3f",
            "release_source": "src/ed_pairs.cpp",
            "relation": "byte-equivalent two-row implementation",
        },
    ]
    if provenance_rows != expected_provenance:
        fail("SOURCE_PROVENANCE.tsv differs from the exact approved four-row ledger")

    expected_status_keys = {
        "actual_completed_task_markers",
        "completed",
        "fixed_pairs",
        "ledger_rebuild_wall_s",
        "random_pairs",
        "retained_failed_attempt_markers",
        "status",
        "task_outer_wall_sum_s",
        "total",
        "truth_only_clusters_excluded",
    }
    if set(final_status) != expected_status_keys:
        fail("FINAL_STATUS.json schema differs from the exact approved field set")
    exact_status_values: Mapping[str, object] = {
        "actual_completed_task_markers": 120,
        "completed": 120,
        "fixed_pairs": 50,
        "random_pairs": 10,
        "retained_failed_attempt_markers": 1,
        "status": "PASS",
        "total": 120,
        "truth_only_clusters_excluded": 680,
    }
    for key, expected in exact_status_values.items():
        actual = final_status.get(key)
        if type(actual) is not type(expected) or actual != expected:
            fail(f"FINAL_STATUS.json has wrong {key}: {actual!r}")
    for key in ("ledger_rebuild_wall_s", "task_outer_wall_sum_s"):
        value = final_status.get(key)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            fail(f"FINAL_STATUS.json has invalid positive timing field: {key}")

    verify_regression_ledgers(payloads)


def check_root(root_arg: str) -> Path:
    root = Path(root_arg)
    reject_symlink_chain(root, "repository root")
    root = root.resolve(strict=True)
    if not root.is_dir():
        fail(f"repository root is not a directory: {root}")
    return root


def walk_worktree(root: Path) -> Tuple[Set[str], Set[str]]:
    files: Set[str] = set()
    directories: Set[str] = set()

    def visit(directory: Path, rel_dir: PurePosixPath) -> None:
        try:
            entries = sorted(os.scandir(directory), key=lambda item: item.name)
        except OSError as exc:
            fail(f"cannot scan {directory}: {exc}")
        for entry in entries:
            rel = (rel_dir / entry.name).as_posix()
            if rel == ".git":
                if entry.is_symlink():
                    fail(".git may not be a symlink")
                continue
            if entry.is_symlink():
                fail(f"symlink forbidden in publication worktree: {rel}")
            if entry.is_dir(follow_symlinks=False):
                directories.add(rel)
                visit(Path(entry.path), rel_dir / entry.name)
            elif entry.is_file(follow_symlinks=False):
                files.add(rel)
            else:
                fail(f"non-regular worktree object: {rel}")

    visit(root, PurePosixPath("."))
    return files, directories


def verify_worktree(root: Path, source_allowed: Set[str]) -> None:
    files, directories = walk_worktree(root)
    release_allowed = {
        f"release/{TOP}-source.tar.gz",
        f"release/{TOP}-linux-x86_64.tar.gz",
        "release/SHA256SUMS.txt",
    }
    actual_source_files = {value for value in files if not value.startswith("release/")}
    extra_source = sorted(actual_source_files - source_allowed)
    missing_source = sorted(source_allowed - actual_source_files)
    if extra_source:
        fail(f"unregistered worktree file: {extra_source[0]}")
    if missing_source:
        fail(f"registered worktree file missing: {missing_source[0]}")
    extra_release = sorted({value for value in files if value.startswith("release/")} - release_allowed)
    if extra_release:
        fail(f"unregistered generated release file: {extra_release[0]}")

    allowed_dirs = parent_directories(source_allowed) | {"release"}
    extra_dirs = sorted(directories - allowed_dirs)
    if extra_dirs:
        fail(f"unregistered worktree directory: {extra_dirs[0]}")

    payloads: Dict[str, bytes] = {}
    for rel in sorted(source_allowed):
        payloads[rel] = read_root_member(root, rel)
        scan_source_path(rel, payloads[rel])
    verify_semantic_contract(payloads)


def archive_member_rel(name: str) -> Tuple[str, bool]:
    safe_relpath(name)
    if name == TOP:
        return "", True
    prefix = TOP + "/"
    if not name.startswith(prefix):
        fail(f"archive has wrong top-level directory: {name}")
    rel = name[len(prefix):]
    if not rel:
        fail(f"non-canonical empty member below top: {name}")
    return safe_relpath(rel), False


def _ustar_cstring(field: bytes, label: str) -> str:
    head, separator, tail = field.partition(b"\x00")
    if separator and any(tail):
        fail(f"non-zero bytes after NUL in ustar {label}")
    try:
        return head.decode("utf-8")
    except UnicodeDecodeError as exc:
        fail(f"non-UTF-8 ustar {label}: {exc}")


def _ustar_octal(field: bytes, label: str) -> int:
    value = field.strip(b" \x00")
    if not value:
        return 0
    if any(byte < ord("0") or byte > ord("7") for byte in value):
        fail(f"non-octal ustar {label}")
    return int(value, 8)


def validate_single_gzip_canonical_ustar(archive: Path) -> Tuple[bytes, List[str]]:
    # Size/type/symlink checks happen before the first byte is read.  The
    # secure reader is bounded, so neither a giant input nor a path swap can
    # force an unbounded read before validation.
    compressed = _read_regular_file_secure(
        archive, MAX_COMPRESSED_ARCHIVE, f"compressed release archive {archive}"
    )
    if len(compressed) < 10 or compressed[:3] != b"\x1f\x8b\x08":
        fail(f"archive is not a gzip-compressed tar: {archive}")
    gzip_flags = compressed[3]
    gzip_mtime = int.from_bytes(compressed[4:8], "little")
    if gzip_flags != 0 or gzip_mtime != 0:
        fail(
            f"non-deterministic gzip header for {archive}: "
            f"flags={gzip_flags} mtime={gzip_mtime}"
        )
    inflater = zlib.decompressobj(wbits=16 + zlib.MAX_WBITS)
    payload_chunks: List[bytes] = []
    expanded = 0
    try:
        for offset in range(0, len(compressed), 1024 * 1024):
            if inflater.eof:
                fail(f"trailing or concatenated gzip payload forbidden: {archive}")
            chunk = compressed[offset : offset + 1024 * 1024]
            remaining = MAX_EXPANDED_ARCHIVE - expanded
            part = inflater.decompress(chunk, remaining + 1)
            if len(part) > remaining or inflater.unconsumed_tail:
                fail(f"gzip-expanded archive exceeds safety limit: {archive}")
            payload_chunks.append(part)
            expanded += len(part)
            if inflater.unused_data:
                fail(f"trailing or concatenated gzip payload forbidden: {archive}")
        remaining = MAX_EXPANDED_ARCHIVE - expanded
        tail = inflater.flush(remaining + 1)
        if len(tail) > remaining:
            fail(f"gzip-expanded archive exceeds safety limit: {archive}")
        payload_chunks.append(tail)
        expanded += len(tail)
    except zlib.error as exc:
        fail(f"invalid gzip stream {archive}: {exc}")
    payload = b"".join(payload_chunks)
    if not inflater.eof:
        fail(f"truncated gzip stream: {archive}")
    if inflater.unused_data:
        fail(f"trailing or concatenated gzip payload forbidden: {archive}")
    if inflater.unconsumed_tail:
        fail(f"unconsumed gzip input forbidden: {archive}")
    if len(payload) % 512 != 0:
        fail(f"ustar payload is not 512-byte aligned: {archive}")

    names: List[str] = []
    offset = 0
    zero_blocks = 0
    while offset + 512 <= len(payload):
        header = payload[offset : offset + 512]
        if header == b"\x00" * 512:
            zero_blocks += 1
            offset += 512
            if zero_blocks >= 2:
                if any(payload[offset:]):
                    fail(f"non-zero data hidden after canonical ustar EOF: {archive}")
                return payload, names
            continue
        if zero_blocks:
            fail(f"non-zero ustar header after an EOF zero block: {archive}")
        if header[257:263] != b"ustar\x00" or header[263:265] != b"00":
            fail(f"archive is not canonical POSIX ustar: {archive}")
        if any(header[500:512]):
            fail(f"non-zero reserved ustar header bytes: {archive}")
        stored_checksum = _ustar_octal(header[148:156], "checksum")
        checksum_header = header[:148] + b" " * 8 + header[156:]
        if sum(checksum_header) != stored_checksum:
            fail(f"ustar header checksum mismatch at byte {offset}: {archive}")
        typeflag = header[156:157]
        if typeflag not in {b"0", b"\x00", b"5"}:
            fail(f"extended/link/special ustar type forbidden: {typeflag!r}")
        if _ustar_cstring(header[157:257], "linkname"):
            fail(f"non-empty ustar linkname forbidden: {archive}")
        _ustar_octal(header[100:108], "mode")
        if _ustar_octal(header[108:116], "uid") != 0 or _ustar_octal(header[116:124], "gid") != 0:
            fail(f"non-zero ustar uid/gid: {archive}")
        if _ustar_octal(header[136:148], "mtime") != RELEASE_MTIME:
            fail(f"non-canonical ustar mtime: {archive}")
        uname = _ustar_cstring(header[265:297], "uname")
        gname = _ustar_cstring(header[297:329], "gname")
        if (uname, gname) not in {("", ""), ("root", "root")}:
            fail(f"non-canonical ustar uname/gname: {(uname, gname)!r}")
        if _ustar_octal(header[329:337], "devmajor") or _ustar_octal(header[337:345], "devminor"):
            fail(f"non-zero ustar device fields forbidden: {archive}")
        name = _ustar_cstring(header[0:100], "name")
        prefix = _ustar_cstring(header[345:500], "prefix")
        full_name = f"{prefix}/{name}" if prefix else name
        safe_relpath(full_name.rstrip("/") if typeflag == b"5" else full_name)
        names.append(full_name.rstrip("/") if typeflag == b"5" else full_name)
        size = _ustar_octal(header[124:136], "size")
        if typeflag == b"5" and size != 0:
            fail(f"ustar directory has non-zero size: {full_name}")
        data_start = offset + 512
        padded_size = ((size + 511) // 512) * 512
        data_end = data_start + size
        padded_end = data_start + padded_size
        if padded_end > len(payload):
            fail(f"truncated ustar member: {full_name}")
        if any(payload[data_end:padded_end]):
            fail(f"non-zero ustar member padding: {full_name}")
        offset = padded_end
    fail(f"canonical two-block ustar EOF missing: {archive}")


def read_archive_regulars(
    archive: Path,
    kind: str,
    allowed: Set[str],
) -> Tuple[Dict[str, bytes], List[str]]:
    if archive.is_symlink() or not archive.is_file():
        fail(f"archive missing or symlinked: {archive}")
    tar_payload, raw_ustar_names = validate_single_gzip_canonical_ustar(archive)
    expected_dirs = parent_directories(allowed)
    expected_names = {TOP} | {f"{TOP}/{rel}" for rel in allowed | expected_dirs}
    seen_names: Set[str] = set()
    seen_casefold: Dict[str, str] = {}
    regulars: Dict[str, bytes] = {}
    member_names: List[str] = []
    total_size = 0
    try:
        handle = tarfile.open(fileobj=io.BytesIO(tar_payload), mode="r:")
    except (tarfile.TarError, OSError) as exc:
        fail(f"cannot open archive {archive}: {exc}")
    with handle:
        for member in handle.getmembers():
            name = member.name.rstrip("/") if member.isdir() else member.name
            if name in seen_names:
                fail(f"duplicate tar member: {name}")
            folded = name.casefold()
            if folded in seen_casefold and seen_casefold[folded] != name:
                fail(
                    f"case-fold-colliding tar members: "
                    f"{seen_casefold[folded]!r} and {name!r}"
                )
            seen_casefold[folded] = name
            seen_names.add(name)
            member_names.append(name)
            rel, is_top = archive_member_rel(name)
            if name not in expected_names:
                fail(f"unregistered {kind} archive member: {name}")
            if member.issym() or member.islnk():
                fail(f"link forbidden in archive: {name}")
            if member.isdev() or member.isfifo() or not (member.isdir() or member.isfile()):
                fail(f"special member forbidden in archive: {name}")
            if member.uid != 0 or member.gid != 0:
                fail(f"non-zero archive owner/group: {name}")
            if member.mtime != RELEASE_MTIME:
                fail(f"non-deterministic mtime for {name}: {member.mtime}")
            expected_member_mode = 0o755 if member.isdir() else expected_mode(rel, kind)
            if stat.S_IMODE(member.mode) != expected_member_mode:
                fail(
                    f"wrong mode for {name}: {oct(stat.S_IMODE(member.mode))}; "
                    f"expected {oct(expected_member_mode)}"
                )
            if member.isdir():
                if not is_top and rel not in expected_dirs:
                    fail(f"unregistered archive directory: {name}")
                continue
            if member.size > MAX_ARCHIVE_MEMBER:
                fail(f"archive member exceeds size limit: {name}")
            total_size += member.size
            if total_size > MAX_ARCHIVE_TOTAL:
                fail(f"archive exceeds total uncompressed size limit: {archive}")
            stream = handle.extractfile(member)
            if stream is None:
                fail(f"cannot read regular archive member: {name}")
            data = stream.read(MAX_ARCHIVE_MEMBER + 1)
            if len(data) != member.size:
                fail(f"short or oversized archive member read: {name}")
            regulars[rel] = data

    if member_names != sorted(member_names):
        fail(f"tar members are not sorted: {archive}")
    if member_names != raw_ustar_names:
        fail(f"tar semantic member list differs from canonical ustar headers: {archive}")
    if seen_names != expected_names:
        missing = sorted(expected_names - seen_names)
        fail(f"archive member set incomplete: {missing[0] if missing else 'unknown'}")
    if set(regulars) != allowed:
        fail(f"archive regular-file set does not equal {kind} allowlist")
    return regulars, member_names


def verify_source_archive(
    root: Path,
    archive: Path,
    source_allowed: Set[str],
    source_manifest_raw: bytes,
) -> None:
    try:
        regulars, _ = read_archive_regulars(archive, "source", source_allowed)
    except VerificationError:
        raise
    except (EOFError, OSError, tarfile.TarError) as exc:
        fail(f"malformed source archive {archive}: {type(exc).__name__}: {exc}")
    if regulars[SOURCE_MANIFEST] != source_manifest_raw:
        fail("source archive allowlist is not identical to its worktree trust anchor")
    for rel, data in sorted(regulars.items()):
        scan_source_path(rel, data)
        worktree_data = read_root_member(root, rel)
        if sha256_bytes(data) != sha256_bytes(worktree_data):
            fail(f"source archive/worktree SHA mismatch: {rel}")
    verify_semantic_contract(regulars)


def verify_binary_archive(
    root: Path,
    archive: Path,
    binary_allowed: Set[str],
) -> None:
    try:
        regulars, _ = read_archive_regulars(archive, "binary", binary_allowed)
    except VerificationError:
        raise
    except (EOFError, OSError, tarfile.TarError) as exc:
        fail(f"malformed binary archive {archive}: {type(exc).__name__}: {exc}")
    expected_elf = {"bin/acor_runner", "bin/ed_pairs"}
    actual_elf = {rel for rel, data in regulars.items() if data.startswith(b"\x7fELF")}
    if actual_elf != expected_elf:
        fail(
            "binary archive ELF set mismatch: "
            f"expected {sorted(expected_elf)}, got {sorted(actual_elf)}"
        )
    for rel, data in sorted(regulars.items()):
        if rel in expected_elf:
            scan_elf_payload(data, rel)
            continue
        path = PurePosixPath(rel)
        if path.suffix.lower() in FORBIDDEN_MODEL_SUFFIXES:
            fail(f"model or pickle asset forbidden in binary release: {rel}")
        scan_text(data, rel, rel == "acor.py")
        if rel not in expected_elf:
            worktree_data = read_root_member(root, rel)
            if sha256_bytes(data) != sha256_bytes(worktree_data):
                fail(f"binary archive/worktree payload mismatch: {rel}")
    verify_semantic_contract(regulars)


def run_internal_self_test() -> None:
    cases = 0

    def must_fail(label: str, operation) -> None:
        nonlocal cases
        cases += 1
        try:
            operation()
        except VerificationError:
            return
        fail(f"internal verifier self-test unexpectedly accepted: {label}")

    must_fail("path traversal", lambda: safe_relpath("../escape"))
    must_fail(
        "model call in shell",
        lambda: scan_source_path(
            "scripts/evil.sh",
            b"python -c 'model.predict_" + b"proba(x)'\n",
        ),
    )
    must_fail(
        "model call hidden in gzip fixture",
        lambda: scan_source_path(
            "tests/golden/tiny_fixed_order.jsonl.gz",
            gzip.compress(b"predict_" + b"proba(x)\n", mtime=0),
        ),
    )
    fake_elf = bytearray(128)
    fake_elf[0:7] = b"\x7fELF\x02\x01\x01"
    fake_elf[16:18] = (1).to_bytes(2, "little")  # ET_REL is forbidden.
    fake_elf[18:20] = (183).to_bytes(2, "little")  # AArch64 is forbidden.
    fake_elf[20:24] = (1).to_bytes(4, "little")
    must_fail("non-x86_64 non-executable ELF", lambda: scan_elf_payload(bytes(fake_elf), "bin/acor_runner"))

    with tempfile.TemporaryDirectory(prefix="acor-verifier-self-test-") as tmp:
        base = Path(tmp)
        target = base / "target"
        target.write_text("target\n", encoding="utf-8")
        link = base / "link"
        try:
            os.symlink(target, link)
        except OSError:
            pass
        else:
            must_fail("symlink path component", lambda: reject_symlink_chain(link, "self-test"))
        canonical_empty_tar = b"\x00" * 1024
        concat = base / "concat.tar.gz"
        concat.write_bytes(
            gzip.compress(canonical_empty_tar, mtime=0)
            + gzip.compress(b"hidden\n", mtime=0)
        )
        must_fail("concatenated gzip member", lambda: validate_single_gzip_canonical_ustar(concat))

        trailer = base / "trailer.tar.gz"
        hidden = b"hidden-after-tar-eof"
        trailer.write_bytes(
            gzip.compress(
                canonical_empty_tar + hidden + b"\x00" * ((-len(hidden)) % 512),
                mtime=0,
            )
        )
        must_fail("non-zero payload after ustar EOF", lambda: validate_single_gzip_canonical_ustar(trailer))

        truncated = base / "truncated.tar.gz"
        valid_gzip = gzip.compress(canonical_empty_tar, mtime=0)
        truncated.write_bytes(valid_gzip[:-8])
        must_fail("truncated gzip", lambda: validate_single_gzip_canonical_ustar(truncated))

    print(f"VERIFIER_SELF_TEST=PASS cases={cases}")


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--worktree", action="store_true", help="audit publication worktree")
    parser.add_argument("--source-tar", type=Path, help="audit deterministic source archive")
    parser.add_argument("--binary-tar", type=Path, help="audit deterministic Linux archive")
    parser.add_argument("--self-test", action="store_true", help="run built-in attack fixtures")
    args = parser.parse_args(argv)
    if not (args.worktree or args.source_tar or args.binary_tar or args.self_test):
        root = Path(args.root).resolve()
        args.worktree = True
        args.source_tar = root / "release" / f"{TOP}-source.tar.gz"
        args.binary_tar = root / "release" / f"{TOP}-linux-x86_64.tar.gz"
    return args


def main(argv: Sequence[str] = ()) -> int:
    args = parse_args(argv or sys.argv[1:])
    if args.self_test:
        run_internal_self_test()
        if not (args.worktree or args.source_tar or args.binary_tar):
            return 0
    root = check_root(args.root)
    source_members, source_manifest_raw = read_manifest(root, SOURCE_MANIFEST)
    binary_members, _ = read_manifest(root, BINARY_MANIFEST)
    source_allowed = set(source_members)
    binary_allowed = set(binary_members)

    version_bytes = read_root_member(root, "VERSION")
    if version_bytes not in {VERSION.encode("ascii"), (VERSION + "\n").encode("ascii")}:
        fail(f"VERSION payload is not exactly {VERSION}")

    missing_source_required = sorted(REQUIRED_SOURCE - source_allowed)
    if missing_source_required:
        fail(f"source allowlist misses required member: {missing_source_required[0]}")
    extra_source_allowed = sorted(source_allowed - REQUIRED_SOURCE)
    if extra_source_allowed:
        fail(f"source allowlist contains unapproved member: {extra_source_allowed[0]}")
    missing_binary_required = sorted(REQUIRED_BINARY - binary_allowed)
    if missing_binary_required:
        fail(f"binary allowlist misses required member: {missing_binary_required[0]}")
    if binary_allowed != REQUIRED_BINARY:
        extra = sorted(binary_allowed - REQUIRED_BINARY)
        fail(f"binary allowlist contains unapproved member: {extra[0]}")

    if args.worktree:
        verify_worktree(root, source_allowed)
        print("WORKTREE_AUDIT=PASS")
    if args.source_tar:
        reject_symlink_chain(args.source_tar, "source archive")
        source_archive = Path(os.path.abspath(args.source_tar))
        verify_source_archive(root, source_archive, source_allowed, source_manifest_raw)
        source_archive_bytes = _read_regular_file_secure(
            source_archive, MAX_COMPRESSED_ARCHIVE, "source archive for SHA-256"
        )
        print(f"SOURCE_ARCHIVE_AUDIT=PASS sha256={sha256_bytes(source_archive_bytes)}")
    if args.binary_tar:
        reject_symlink_chain(args.binary_tar, "binary archive")
        binary_archive = Path(os.path.abspath(args.binary_tar))
        verify_binary_archive(root, binary_archive, binary_allowed)
        binary_archive_bytes = _read_regular_file_secure(
            binary_archive, MAX_COMPRESSED_ARCHIVE, "binary archive for SHA-256"
        )
        print(f"BINARY_ARCHIVE_AUDIT=PASS sha256={sha256_bytes(binary_archive_bytes)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except VerificationError as exc:
        print(f"RELEASE_AUDIT=FAIL reason={exc}", file=sys.stderr)
        raise SystemExit(2)
    except (EOFError, OSError, tarfile.TarError, UnicodeError, zlib.error) as exc:
        print(
            f"RELEASE_AUDIT=FAIL reason={type(exc).__name__}: {exc}",
            file=sys.stderr,
        )
        raise SystemExit(2)
