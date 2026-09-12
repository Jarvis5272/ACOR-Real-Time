#!/usr/bin/env bash
set -euo pipefail

profile="${1:-portable-release}"
build_dir="build/${profile}"
cmake -S . -B "${build_dir}" -DCMAKE_BUILD_TYPE=Release -DACOR_BUILD_PROFILE="${profile}"
cmake --build "${build_dir}" --parallel
cmake --install "${build_dir}" --prefix .
