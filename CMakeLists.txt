cmake_minimum_required(VERSION 3.16)

project(acor-reconstruction-backbone VERSION 0.2.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(CTest)

find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

set(ACOR_BUILD_PROFILE "portable-release" CACHE STRING
    "Release compiler profile: portable-release or frozen-compatible")
set_property(CACHE ACOR_BUILD_PROFILE PROPERTY STRINGS
             portable-release frozen-compatible)

if(CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_CONFIGURATION_TYPES "Release" CACHE STRING
      "ACOR release package supports the Release configuration" FORCE)
else()
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "ACOR release package build type" FORCE)
  elseif(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
    message(FATAL_ERROR
      "ACOR_BUILD_PROFILE is a release-package contract; use -DCMAKE_BUILD_TYPE=Release")
  endif()
endif()

if(NOT ACOR_BUILD_PROFILE STREQUAL "portable-release" AND
   NOT ACOR_BUILD_PROFILE STREQUAL "frozen-compatible")
  message(FATAL_ERROR
          "ACOR_BUILD_PROFILE must be portable-release or frozen-compatible")
endif()

# These two profiles select compiler engineering only.  They do not change the
# ACOR algorithm, its constants, or its input/output semantics.
#
# portable-release matches the effective -O3 -DNDEBUG compilation flags
# exercised by the 120-run release regression.
# frozen-compatible documents the historical reference compilation flags.  It
# is intentionally described as a reference profile; rebuilding new source is
# not claimed to reproduce a historical executable byte-for-byte.
if(CMAKE_INTERPROCEDURAL_OPTIMIZATION OR
   CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE)
  message(FATAL_ERROR
    "ACOR release profiles forbid interprocedural/LTO optimisation")
endif()
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF CACHE BOOL
    "LTO is disabled for ACOR release profiles" FORCE)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE OFF CACHE BOOL
    "Release LTO is disabled for ACOR release profiles" FORCE)

if(MSVC)
  string(CONCAT _acor_inherited_flags
    "${CMAKE_CXX_FLAGS} "
    "${CMAKE_CXX_FLAGS_RELEASE} "
    "${CMAKE_EXE_LINKER_FLAGS} "
    "${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
  string(TOLOWER "${_acor_inherited_flags}" _acor_inherited_flags_lower)
  if(_acor_inherited_flags_lower MATCHES "(^|[ ,;])/arch:" OR
     _acor_inherited_flags_lower MATCHES "(^|[ ,;])/gl($|[ ,;])")
    message(FATAL_ERROR
      "ACOR release profiles forbid architecture-specific /arch and /GL flags")
  endif()
  if(ACOR_BUILD_PROFILE STREQUAL "portable-release")
    set(_acor_release_flags "/O2 /DNDEBUG")
  else()
    # MSVC has no literal GCC -O2 equivalent.  /O2 is the closest supported
    # reference setting and remains architecture-neutral.
    set(_acor_release_flags "/O2 /DNDEBUG")
  endif()
  set(CMAKE_CXX_FLAGS_RELEASE "${_acor_release_flags}" CACHE STRING
      "ACOR release flags selected by ACOR_BUILD_PROFILE" FORCE)
else()
  # Both public profiles are architecture-neutral and non-LTO.  Check all
  # inherited compile/link variables before replacing the Release defaults so
  # an environment CXXFLAGS/LDFLAGS injection cannot be hidden by the cache.
  string(CONCAT _acor_inherited_flags
    "${CMAKE_CXX_FLAGS} "
    "${CMAKE_CXX_FLAGS_RELEASE} "
    "${CMAKE_EXE_LINKER_FLAGS} "
    "${CMAKE_EXE_LINKER_FLAGS_RELEASE} "
    "${CMAKE_SHARED_LINKER_FLAGS} "
    "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} "
    "${CMAKE_MODULE_LINKER_FLAGS} "
    "${CMAKE_MODULE_LINKER_FLAGS_RELEASE}")
  if(_acor_inherited_flags MATCHES "(^|[ ,;])-march(=|[ ,;])" OR
     _acor_inherited_flags MATCHES "(^|[ ,;])-mtune(=|[ ,;])" OR
     _acor_inherited_flags MATCHES "(^|[ ,;])-flto($|[ =,;])")
    message(FATAL_ERROR
      "ACOR release profiles forbid -march, -mtune, and -flto; clear inherited flags")
  endif()
  if(ACOR_BUILD_PROFILE STREQUAL "portable-release")
    set(_acor_release_flags "-O3 -DNDEBUG")
  else()
    set(_acor_release_flags "-O2 -DNDEBUG")
  endif()

  # Replace CMake's platform Release defaults instead of stacking -O2 after
  # an implicit -O3.  This makes the selected profile unambiguous in verbose
  # compiler logs and prevents accidental -march/-flto additions.
  set(CMAKE_CXX_FLAGS_RELEASE "${_acor_release_flags}" CACHE STRING
      "ACOR release flags selected by ACOR_BUILD_PROFILE" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(acor_runner
  src/acor_runner.cpp
  src/acor_kernel.cpp
)
target_include_directories(acor_runner PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(acor_runner PRIVATE Threads::Threads ZLIB::ZLIB)

add_executable(ed_pairs
  src/ed_pairs.cpp
)

add_executable(test_ledger_reduction
  tests/test_ledger_reduction.cpp
  src/acor_kernel.cpp
)
target_include_directories(test_ledger_reduction PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(test_ledger_reduction PRIVATE Threads::Threads ZLIB::ZLIB)

if(MSVC)
  target_compile_options(acor_runner PRIVATE /W4 /permissive-)
  target_compile_options(ed_pairs PRIVATE /W4 /permissive-)
else()
  target_compile_options(acor_runner PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(ed_pairs PRIVATE -Wall -Wextra -Wpedantic)
endif()

install(TARGETS acor_runner ed_pairs
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(PROGRAMS acor.py DESTINATION .)
install(FILES README.md LICENSE CITATION.cff VERSION CHANGELOG.md
              THIRD_PARTY_NOTICES.md
        DESTINATION .)
install(DIRECTORY provenance/ DESTINATION provenance)

if(BUILD_TESTING)
  find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter)
  set(ACOR_TEST_PACKAGE_ROOT "${CMAKE_CURRENT_BINARY_DIR}/test-package")
  add_custom_target(acor_test_package ALL
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${ACOR_TEST_PACKAGE_ROOT}/bin"
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/acor.py" "${ACOR_TEST_PACKAGE_ROOT}/acor.py"
    COMMAND "${CMAKE_COMMAND}" -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/src" "${ACOR_TEST_PACKAGE_ROOT}/src"
    COMMAND "${CMAKE_COMMAND}" -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/tests" "${ACOR_TEST_PACKAGE_ROOT}/tests"
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different
            "$<TARGET_FILE:acor_runner>" "${ACOR_TEST_PACKAGE_ROOT}/bin/acor_runner"
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different
            "$<TARGET_FILE:ed_pairs>" "${ACOR_TEST_PACKAGE_ROOT}/bin/ed_pairs"
    DEPENDS acor_runner ed_pairs
    VERBATIM
  )
  add_test(
    NAME acor_python_unittest
    COMMAND "${CMAKE_COMMAND}" -E env
            "ACOR_RUNNER=$<TARGET_FILE:acor_runner>"
            "ACOR_RUNNER_BIN=$<TARGET_FILE:acor_runner>"
            "ACOR_ED_PAIRS=$<TARGET_FILE:ed_pairs>"
            "ACOR_ED_PAIRS_BIN=$<TARGET_FILE:ed_pairs>"
            "ACOR_PACKAGE_BIN_DIR=$<TARGET_FILE_DIR:acor_runner>"
            "ACOR_PACKAGE_ROOT=${ACOR_TEST_PACKAGE_ROOT}"
            "${Python3_EXECUTABLE}" -m unittest discover
            -s "${ACOR_TEST_PACKAGE_ROOT}/tests" -p "test_*.py" -v
  )
  set_tests_properties(acor_python_unittest PROPERTIES
    WORKING_DIRECTORY "${ACOR_TEST_PACKAGE_ROOT}")
  add_test(NAME acor_ledger_reduction COMMAND test_ledger_reduction)
endif()

message(STATUS "ACOR build profile: ${ACOR_BUILD_PROFILE}")
message(STATUS "Release flags: ${_acor_release_flags}")
if(ACOR_BUILD_PROFILE STREQUAL "frozen-compatible")
  message(STATUS
    "frozen-compatible supplies reference flags only; no historical binary SHA claim is made")
endif()
