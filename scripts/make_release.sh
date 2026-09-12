#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

readonly RELEASE_VERSION="0.1.0"
readonly RELEASE_TOP="acor-reconstruction-backbone-${RELEASE_VERSION}"
readonly RELEASE_MTIME="1704067200"
readonly SMOKE_SEED="0123456789abcdef0123456789abcdef"

die() {
  printf 'release error: %s\n' "$*" >&2
  exit 1
}

[[ ! -L "${BASH_SOURCE[0]}" ]] || die "make_release.sh may not be invoked through a symlink"
script_dir="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(CDPATH= cd -- "${script_dir}/.." && pwd -P)"
source_manifest="${script_dir}/release_source_members.txt"
binary_manifest="${script_dir}/release_binary_members.txt"
verify_script="${script_dir}/verify_release.py"
output_dir="${repo_root}/release"

if [[ $# -gt 0 ]]; then
  if [[ $# -ne 2 || "$1" != "--output-dir" ]]; then
    die "usage: scripts/make_release.sh [--output-dir DIR]"
  fi
  output_dir="$2"
fi

for tool in cmake cmp ctest env flock gzip make python3 sha256sum tar uname; do
  command -v "${tool}" >/dev/null 2>&1 || die "required tool not found: ${tool}"
done

[[ "$(uname -s)" == "Linux" && "$(uname -m)" == "x86_64" ]] ||
  die "the linux-x86_64 binary release must be built on Linux x86_64"

[[ -f "${repo_root}/VERSION" ]] || die "VERSION is missing"
[[ "$(tr -d '\r\n' < "${repo_root}/VERSION")" == "${RELEASE_VERSION}" ]] ||
  die "VERSION does not equal ${RELEASE_VERSION}"
[[ -f "${source_manifest}" && -f "${binary_manifest}" ]] ||
  die "release member manifests are missing"
[[ -f "${verify_script}" && ! -L "${verify_script}" ]] ||
  die "verify_release.py is missing or symlinked"

# A draft allowlist carrying this marker is deliberately non-releasable.
if grep -q '^# UNFINALIZED_ALLOWLIST$' "${source_manifest}" "${binary_manifest}"; then
  die "release member allowlists are still marked UNFINALIZED_ALLOWLIST"
fi

[[ ! -L "${output_dir}" ]] || die "output directory may not be a symlink"
mkdir -p -- "${output_dir}"
output_dir="$(CDPATH= cd -- "${output_dir}" && pwd -P)"
if [[ "${output_dir}" == "${repo_root}" ||
      ( "${output_dir}" == "${repo_root}"/* && "${output_dir}" != "${repo_root}/release" ) ]]; then
  die "an in-repository output directory must be exactly ${repo_root}/release"
fi

# The output directory, rather than the checkout, is the shared publication
# resource.  Lock its actual directory inode (not a predictable /tmp pathname),
# so two different checkouts targeting one destination serialize safely and a
# pre-planted lock-file symlink cannot redirect/truncate another file.
exec 9<"${output_dir}" || die "cannot open output directory for locking"
flock -n 9 || die "another release build holds the output-directory lock"

tmp_root="$(mktemp -d "${TMPDIR:-/tmp}/acor-release.XXXXXXXX")"
publish_stage=""
cleanup() {
  chmod -R u+w "${tmp_root}" 2>/dev/null || true
  rm -rf -- "${tmp_root}"
  if [[ -n "${publish_stage}" && -d "${publish_stage}" ]]; then
    chmod -R u+w "${publish_stage}" 2>/dev/null || true
    rm -rf -- "${publish_stage}"
  fi
}
trap cleanup EXIT HUP INT TERM

source_stage_parent="${tmp_root}/source-stage"
source_stage="${source_stage_parent}/${RELEASE_TOP}"
source_unpack_parent="${tmp_root}/source-unpack"
build_dir="${tmp_root}/portable-build"
install_parent="${tmp_root}/binary-stage"
install_dir="${install_parent}/${RELEASE_TOP}"
mkdir -p -- "${source_stage}" "${source_unpack_parent}" "${build_dir}" "${install_dir}"

read_manifest() {
  local manifest="$1"
  sed -e 's/\r$//' -e '/^[[:space:]]*#/d' -e '/^[[:space:]]*$/d' "${manifest}"
}

source_mode() {
  case "$1" in
    acor.py|scripts/*.sh|scripts/verify_release.py|regression/*.py)
      printf '0755\n'
      ;;
    *)
      printf '0644\n'
      ;;
  esac
}

copy_registered_file() {
  local from_root="$1"
  local to_root="$2"
  local rel="$3"
  local mode="$4"
  local src="${from_root}/${rel}"
  local dst="${to_root}/${rel}"

  [[ "${rel}" != /* && "${rel}" != *'\\'* && "${rel}" != *'..'* ]] ||
    die "unsafe member path: ${rel}"
  [[ -f "${src}" && ! -L "${src}" ]] ||
    die "registered file missing or not a regular file: ${rel}"
  install -D -m "${mode}" -- "${src}" "${dst}"
}

snapshot_registered_files() {
  local manifest="$1"
  local root="$2"
  local out="$3"
  : > "${out}"
  while IFS= read -r rel; do
    [[ -n "${rel}" ]] || continue
    [[ -f "${root}/${rel}" && ! -L "${root}/${rel}" ]] ||
      die "cannot snapshot registered file: ${rel}"
    sha256sum -- "${root}/${rel}" | sed "s#  ${root}/#  #" >> "${out}"
  done < <(read_manifest "${manifest}")
  LC_ALL=C sort -o "${out}" "${out}"
}

normalise_stage() {
  local root="$1"
  find "${root}" -type d -exec chmod 0755 {} +
  find "${root}" -exec touch -h -d "@${RELEASE_MTIME}" {} +
}

make_deterministic_tar() {
  local parent="$1"
  local archive="$2"
  local tmp_archive="${archive}.tmp"
  rm -f -- "${tmp_archive}"
  (
    cd -- "${parent}"
    LC_ALL=C tar \
      --sort=name \
      --format=ustar \
      --mtime="@${RELEASE_MTIME}" \
      --owner=0 --group=0 --numeric-owner \
      --mode='u+rwX,go+rX,go-w' \
      -cf - "${RELEASE_TOP}"
  ) | gzip -n -9 > "${tmp_archive}"
  mv -f -- "${tmp_archive}" "${archive}"
}

snapshot_before="${tmp_root}/source-before.sha256"
snapshot_after="${tmp_root}/source-after.sha256"

# Strictly audit manifests and the publication worktree before reading any
# manifest-selected path.
python3 -I "${verify_script}" --root "${repo_root}" --worktree
snapshot_registered_files "${source_manifest}" "${repo_root}" "${snapshot_before}"

while IFS= read -r rel; do
  [[ -n "${rel}" ]] || continue
  copy_registered_file "${repo_root}" "${source_stage}" "${rel}" "$(source_mode "${rel}")"
done < <(read_manifest "${source_manifest}")
normalise_stage "${source_stage_parent}"

source_a="${tmp_root}/source-a.tar.gz"
source_b="${tmp_root}/source-b.tar.gz"
make_deterministic_tar "${source_stage_parent}" "${source_a}"
make_deterministic_tar "${source_stage_parent}" "${source_b}"
source_a_sha256="$(sha256sum "${source_a}" | cut -d' ' -f1)"
source_b_sha256="$(sha256sum "${source_b}" | cut -d' ' -f1)"
[[ "${source_a_sha256}" == "${source_b_sha256}" ]] ||
  die "two source archive passes produced different SHA-256 values"
python3 -I "${verify_script}" --root "${repo_root}" --source-tar "${source_a}"

tar -xzf "${source_a}" -C "${source_unpack_parent}"
# Release inputs must not inherit machine-specific or LTO flags.  CMake also
# rejects explicit cache injection, while this clean environment closes the
# compiler/linker environment channel.
env -u CC -u CXX -u CFLAGS -u CXXFLAGS -u CPPFLAGS -u LDFLAGS \
    -u CMAKE_TOOLCHAIN_FILE -u CMAKE_GENERATOR -u CMAKE_GENERATOR_PLATFORM \
    -u CMAKE_GENERATOR_TOOLSET \
    -u CMAKE_CXX_FLAGS -u CMAKE_EXE_LINKER_FLAGS \
    -u CMAKE_SHARED_LINKER_FLAGS -u CMAKE_MODULE_LINKER_FLAGS \
  cmake -S "${source_unpack_parent}/${RELEASE_TOP}" -B "${build_dir}" \
  -G "Unix Makefiles" \
  -DCMAKE_BUILD_TYPE=Release \
  -DACOR_BUILD_PROFILE=portable-release \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DBUILD_TESTING=ON
python3 -I -c 'import json,pathlib,sys; rows=json.loads(pathlib.Path(sys.argv[1]).read_text()); cmds=[row["command"] for row in rows]; forbidden=("-march","-mtune","-flto"); ok=bool(cmds) and all("-O3" in command and "-DNDEBUG" in command for command in cmds) and not any(flag in command for command in cmds for flag in forbidden); print("PORTABLE_RELEASE_FLAGS=PASS") if ok else sys.exit("invalid portable-release compile commands")' \
  "${build_dir}/compile_commands.json" ||
  die "portable clean-room compile flags are not exactly architecture-neutral -O3/-DNDEBUG"
cmake --build "${build_dir}" --parallel
ctest --test-dir "${build_dir}" --output-on-failure
cmake --install "${build_dir}" --prefix "${install_dir}"

# The install tree must contain exactly the binary allowlist.  Re-copying
# documentation/provenance is explicit and does not auto-discover files.
while IFS= read -r rel; do
  [[ -n "${rel}" ]] || continue
  if [[ ! -e "${install_dir}/${rel}" ]]; then
    case "${rel}" in
      bin/acor_runner|bin/ed_pairs)
        die "installed executable missing: ${rel}"
        ;;
      *)
        copy_registered_file "${source_stage}" "${install_dir}" "${rel}" "$(source_mode "${rel}")"
        ;;
    esac
  fi
done < <(read_manifest "${binary_manifest}")
normalise_stage "${install_parent}"

binary_a="${tmp_root}/binary-a.tar.gz"
binary_b="${tmp_root}/binary-b.tar.gz"
make_deterministic_tar "${install_parent}" "${binary_a}"
make_deterministic_tar "${install_parent}" "${binary_b}"
binary_a_sha256="$(sha256sum "${binary_a}" | cut -d' ' -f1)"
binary_b_sha256="$(sha256sum "${binary_b}" | cut -d' ' -f1)"
[[ "${binary_a_sha256}" == "${binary_b_sha256}" ]] ||
  die "two binary archive passes produced different SHA-256 values"
python3 -I "${verify_script}" --root "${repo_root}" --binary-tar "${binary_a}"

# Exercise the exact binary-archive layout, runner, evaluator, CLI, and bundled
# tiny fixture before any publication file is committed.
binary_smoke_parent="${tmp_root}/binary-smoke-unpack"
binary_smoke_out="${tmp_root}/binary-smoke-results"
mkdir -p -- "${binary_smoke_parent}"
tar -xzf "${binary_a}" -C "${binary_smoke_parent}"
python3 "${binary_smoke_parent}/${RELEASE_TOP}/acor.py" run \
  --data "${binary_smoke_parent}/${RELEASE_TOP}/tests/fixtures/tiny_aim" \
  --aim --threads 1 --seeds "${SMOKE_SEED}" \
  --results "${binary_smoke_out}"
binary_smoke_pred="$(find "${binary_smoke_out}" -type f -name pred.tsv -print -quit)"
binary_smoke_summary="$(find "${binary_smoke_out}" -type f -name SUMMARY.tsv -print -quit)"
test -s "${binary_smoke_pred}" ||
  die "binary archive smoke did not produce pred.tsv"
test -s "${binary_smoke_summary}" ||
  die "binary archive smoke did not produce SUMMARY.tsv"
cmp -s "${binary_smoke_pred}" "${source_stage}/tests/golden/tiny_fixed_pred.tsv" ||
  die "binary archive smoke prediction differs from the registered golden output"
python3 -I -c 'import csv,sys; rows=list(csv.DictReader(open(sys.argv[1],encoding="utf-8"),delimiter="\t")); metrics=[float(rows[-1][key]) for key in ("exact","accuracy","mean_ed")] if rows else []; pred_lines=sum(1 for _ in open(sys.argv[2],encoding="utf-8")); ok=bool(rows) and len(metrics)==3 and pred_lines>1; print("BINARY_SMOKE_SCHEMA=PASS") if ok else sys.exit("binary smoke schema/metrics invalid")' \
  "${binary_smoke_summary}" "${binary_smoke_pred}" ||
  die "binary archive smoke metrics/prediction schema validation failed"

# Fail closed if a registered input changed while configure/build/test ran.
snapshot_registered_files "${source_manifest}" "${repo_root}" "${snapshot_after}"
cmp -s "${snapshot_before}" "${snapshot_after}" ||
  die "registered source files changed during release creation"
# Re-walk the full tree as the last pre-publication read.  This catches a new
# unregistered file/directory created after the initial audit, not only edits
# to manifest-selected files.
python3 -I "${verify_script}" --root "${repo_root}" --worktree

source_final="${output_dir}/${RELEASE_TOP}-source.tar.gz"
binary_final="${output_dir}/${RELEASE_TOP}-linux-x86_64.tar.gz"
publish_stage="$(mktemp -d "${output_dir}/.release-stage.XXXXXXXX")"
source_staged="${publish_stage}/$(basename -- "${source_final}")"
binary_staged="${publish_stage}/$(basename -- "${binary_final}")"
install -m 0644 -- "${source_a}" "${source_staged}"
install -m 0644 -- "${binary_a}" "${binary_staged}"
(
  cd -- "${publish_stage}"
  sha256sum -- "$(basename -- "${source_staged}")" \
                "$(basename -- "${binary_staged}")" > SHA256SUMS.txt
)

python3 -I "${verify_script}" --root "${repo_root}" \
  --source-tar "${source_staged}" --binary-tar "${binary_staged}"

# SHA256SUMS is the commit marker: remove the old marker first, atomically
# replace each named archive, then install the new marker last.  A crash can
# therefore yield only a fail-closed missing/mismatching marker, never a false
# successful release.
rm -f -- "${output_dir}/SHA256SUMS.txt"
mv -f -- "${source_staged}" "${source_final}"
mv -f -- "${binary_staged}" "${binary_final}"
mv -f -- "${publish_stage}/SHA256SUMS.txt" "${output_dir}/SHA256SUMS.txt"
rmdir -- "${publish_stage}"
publish_stage=""

python3 -I "${verify_script}" --root "${repo_root}" \
  --source-tar "${source_final}" --binary-tar "${binary_final}"
(
  cd -- "${output_dir}"
  sha256sum --check SHA256SUMS.txt
)

printf 'source  %s\n' "$(sha256sum "${source_final}")"
printf 'binary  %s\n' "$(sha256sum "${binary_final}")"
printf 'release verification: PASS\n'
