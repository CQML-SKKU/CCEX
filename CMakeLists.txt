# CCEX: the CMake build. The same binaries as do_compile.sh (main.out, and
# main_gpu.out where there is CUDA), found the standard way:
#
#     cmake -B build                     # finds the MPI wrapper and, if there is one, nvcc
#     cmake --build build -j 8           # -> build/bin/main.out
#
# Options (-D on the configure line):
#     CCEX_CUDA=AUTO|ON|OFF   the CUDA backend binary; AUTO builds it when a CUDA compiler of 13.2 or newer is found
#                             (any version when CMAKE_CUDA_ARCHITECTURES names only capabilities below 120),
#                             ON stops on an older one
#     CCEX_CUDA_SPLIT_COMPILE=N  nvcc -split-compile threads for the CUDA files (5, do_compile.sh's CCEX_JOBS default;
#                             0 every CPU, 1 off): the optimizer takes a file's kernels at once
#     CCEX_NATIVE=ON|OFF      -march=native on the CPU objects, -mcpu=native on Apple Silicon (ON)
#     CCEX_SANITIZE=ON|OFF    AddressSanitizer and UBSan on the CPU, adaptive and SQLite objects (OFF)
#     CCEX_BLAS=AUTO|MKL|ACCELERATE|OFF  large cCCE/gCCE clusters by Intel MKL (zgemm in the cCCE ensemble Hahn echo's
#                             recursive path from bath dimension 21, zheevd from 63: the bath dimension for cCCE, qubit
#                             x bath for gCCE; include/common/blas_option.h), as do_compile.sh's CCEX_BLAS. AUTO (default): on when
#                             oneAPI's MKLConfig.cmake is found (MKLROOT from setvars.sh or the intel module), on macOS
#                             Apple's Accelerate; MKL, ACCELERATE (macOS): stop without it; OFF: never. lp64, sequential
#     CCEX_BIN_DIR=<dir>      where the executables go (default <build>/bin); -DCCEX_BIN_DIR=$PWD/bin
#                             puts them where do_compile.sh does and the checks look by default
#     CMAKE_CUDA_ARCHITECTURES=<cc>   the card's compute capability, e.g. 86 or 120 (default: native)
#     MPI_CXX_COMPILER=<wrapper>      when the wrapper is not the mpicxx on PATH, e.g. mpiicpc
#     CMAKE_CXX_COMPILER=icpx         Intel's compiler for the C++ (idm4-gpu: 15-25% faster than gcc 11 on
#                                     the CPU path); with CUDA add CMAKE_CUDA_HOST_COMPILER=g++, because
#                                     nvcc 12.8 accepts icpx as host compiler only from oneAPI 2024 on
#
# The binaries stay in the build folder unless CCEX_BIN_DIR says otherwise. Until
# 2026-10-01 they went to bin/ of the source tree, which do_compile.sh writes too, so a
# CMake build and a do_compile.sh build (or two CMake builds with other flags) replaced
# each other's main.out. The checks take the binary to test from CCEX_BIN, e.g.
# CCEX_BIN=$PWD/build/bin/main.out, and the unit-test drivers of check_rotation and
# check_pcce take the build folder from CCEX_BUILD_DIR.
# Eigen, uthash and SQLite are in external/. MKL when it is found (CCEX_BLAS, see do_compile.sh).

cmake_minimum_required(VERSION 3.18)
project(CCEX VERSION 2.0.1 LANGUAGES C CXX)

set(CCEX_CUDA AUTO CACHE STRING "CUDA backend binary: AUTO, ON or OFF")
set(CCEX_CUDA_SPLIT_COMPILE 5 CACHE STRING "nvcc -split-compile threads for the CUDA files (0: every CPU, 1: off)")
option(CCEX_NATIVE "-march=native on the CPU objects (-mcpu=native on Apple Silicon)" ON)
option(CCEX_SANITIZE "AddressSanitizer and UBSan on the CPU, adaptive and SQLite objects" OFF)
set(CCEX_BLAS AUTO CACHE STRING "BLAS/LAPACK for the large clusters: AUTO (when found), MKL (required), ACCELERATE (macOS) or OFF")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)                     # gnu++17, as the presets
set(CCEX_BIN_DIR ${CMAKE_BINARY_DIR}/bin CACHE PATH "where the executables go")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CCEX_BIN_DIR})
# Optimisation is set per target below (the refinement policies keep their own flags),
# so a build type adds nothing on top.
set(CMAKE_CXX_FLAGS_RELEASE "" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE "" CACHE STRING "" FORCE)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

find_package(MPI REQUIRED COMPONENTS CXX)
find_package(Threads REQUIRED)

# ---- CCEX_BLAS (A46): with MKL every target compiles with CCEX_BLAS_DIRECT and links MKL (lp64, sequential: one thread
#      per MPI rank), set before the targets below so that the binaries and the unit tests get the same
set(_ccex_blas_unit_libs)
set(_ccex_mkl OFF)
set(_ccex_accelerate OFF)
string(TOUPPER "${CCEX_BLAS}" _ccex_blas)
# macOS: Apple's Accelerate in MKL's place (do_compile.sh's mac preset), on every Mac with the system: the same Fortran
# zgemm_ and zheevd_ (src/common/blas_option.cpp). AUTO looked for MKL there too until 2026-10-19 and found Intel's macOS
# MKL 2023 under /opt/intel, x86_64 only, which the arm64 link then refused.
if(APPLE AND (_ccex_blas STREQUAL "AUTO" OR _ccex_blas STREQUAL "ACCELERATE"))
  find_library(CCEX_ACCELERATE_FRAMEWORK Accelerate REQUIRED)
  set(_ccex_accelerate ON)
elseif(_ccex_blas STREQUAL "ACCELERATE")
  message(FATAL_ERROR "CCEX: CCEX_BLAS=ACCELERATE: Apple's Accelerate is on macOS only")
elseif(APPLE AND _ccex_blas STREQUAL "MKL")
  message(FATAL_ERROR "CCEX: CCEX_BLAS=MKL: no MKL on macOS; ACCELERATE (the default here) or OFF")
elseif(_ccex_blas STREQUAL "AUTO" OR _ccex_blas STREQUAL "MKL")
  set(MKL_INTERFACE lp64)
  set(MKL_INTERFACE_FULL intel_lp64)            # the name older MKLConfig.cmake (2023) reads
  set(MKL_THREADING sequential)
  set(MKL_LINK dynamic)
  if(_ccex_blas STREQUAL "MKL")
    find_package(MKL CONFIG REQUIRED)
  else()
    find_package(MKL CONFIG QUIET)
  endif()
  set(_ccex_mkl ${MKL_FOUND})
elseif(NOT _ccex_blas STREQUAL "OFF")
  message(FATAL_ERROR "CCEX: CCEX_BLAS=${CCEX_BLAS}: expected AUTO, MKL, ACCELERATE or OFF")
endif()
if(_ccex_accelerate)
  add_compile_definitions(CCEX_BLAS_DIRECT CCEX_BLAS_ACCELERATE)
  link_libraries(${CCEX_ACCELERATE_FRAMEWORK})
  set(_ccex_blas_unit_libs -framework Accelerate)
  message(STATUS "CCEX: Accelerate on (${CCEX_ACCELERATE_FRAMEWORK})")
elseif(_ccex_mkl)
  add_compile_definitions(CCEX_BLAS_DIRECT)
  link_libraries(MKL::MKL)
  set(_mkl_lib ${MKL_ROOT}/lib/intel64)
  if(NOT IS_DIRECTORY ${_mkl_lib})
    set(_mkl_lib ${MKL_ROOT}/lib)
  endif()
  set(_ccex_blas_unit_libs -L${_mkl_lib} -Wl,-rpath,${_mkl_lib} -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core)
  message(STATUS "CCEX: MKL on (${_mkl_lib})")
elseif(_ccex_blas STREQUAL "OFF")
  message(STATUS "CCEX: MKL off (CCEX_BLAS=OFF)")
else()
  message(STATUS "CCEX: MKL off (MKLConfig.cmake not found; source setvars.sh or set MKL_DIR)")
endif()

# ---- CUDA: AUTO means "when nvcc is found". 12.8 miscompiles the pulse kernel applyPulse for compute capability
# 12.0 (A97, CHANGELOG.md), so CUDA 13.2 or newer is taken unless CMAKE_CUDA_ARCHITECTURES names only capabilities
# below 12.0 (e.g. 90 for GH200, whose driver for CUDA 12.x runs 13.2 only in compat mode, without MPS); make gpu-test
# shows an older nvcc compiles applyPulse right there. AUTO leaves a refused nvcc out with a warning; ON stops on it.
set(CCEX_CUDA_MIN 13.2)
set(_ccex_cuda_old_ok OFF)   # every named architecture below 120 ("native", "all" or none named: not known here)
foreach(_a ${CMAKE_CUDA_ARCHITECTURES})
  string(REGEX MATCH "^[0-9]+" _n "${_a}")
  if(_n AND _n LESS 120)
    set(_ccex_cuda_old_ok ON)
  else()
    set(_ccex_cuda_old_ok OFF)
    break()
  endif()
endforeach()
include(CheckLanguage)
if(CCEX_CUDA STREQUAL "AUTO")
  check_language(CUDA)
  if(CMAKE_CUDA_COMPILER)
    execute_process(COMMAND ${CMAKE_CUDA_COMPILER} --version OUTPUT_VARIABLE _ccex_nvcc ERROR_QUIET)
    if(_ccex_nvcc MATCHES "release ([0-9]+\\.[0-9]+)")
      set(_ccex_nvcc ${CMAKE_MATCH_1})
    else()
      set(_ccex_nvcc 0)
    endif()
    if(_ccex_nvcc VERSION_LESS CCEX_CUDA_MIN AND NOT _ccex_cuda_old_ok)
      message(WARNING "CCEX: ${CMAKE_CUDA_COMPILER} is CUDA ${_ccex_nvcc}; the CUDA backend needs ${CCEX_CUDA_MIN} or newer unless CMAKE_CUDA_ARCHITECTURES names only capabilities below 120, so main_gpu.out is not built (-DCMAKE_CUDA_COMPILER=<newer nvcc>, -DCMAKE_CUDA_ARCHITECTURES=90, or -DCCEX_CUDA=OFF)")
      set(_ccex_cuda OFF)
    else()
      set(_ccex_cuda ON)
    endif()
  else()
    set(_ccex_cuda OFF)
  endif()
elseif(CCEX_CUDA)
  set(_ccex_cuda ON)
else()
  set(_ccex_cuda OFF)
endif()
if(_ccex_cuda)
  if(NOT CMAKE_CUDA_ARCHITECTURES)
    # "native" is a CMake 3.24 value; older CMake (3.18 is enough for the CPU build) needs the number
    if(CMAKE_VERSION VERSION_LESS 3.24)
      message(FATAL_ERROR "CCEX: CMake ${CMAKE_VERSION} cannot detect the card; give -DCMAKE_CUDA_ARCHITECTURES=<cc> (e.g. 86, 120) or -DCCEX_CUDA=OFF")
    endif()
    set(CMAKE_CUDA_ARCHITECTURES native)
  endif()
  enable_language(CUDA)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS CCEX_CUDA_MIN)
    if(NOT _ccex_cuda_old_ok)
      message(FATAL_ERROR "CCEX: ${CMAKE_CUDA_COMPILER} is CUDA ${CMAKE_CUDA_COMPILER_VERSION}; the CUDA backend needs ${CCEX_CUDA_MIN} or newer (12.8 miscompiles applyPulse for compute capability 12.0, see CHANGELOG.md) unless CMAKE_CUDA_ARCHITECTURES names only capabilities below 120. Give -DCMAKE_CUDA_COMPILER=<CUDA ${CCEX_CUDA_MIN} or newer>/bin/nvcc, -DCMAKE_CUDA_ARCHITECTURES=<cc below 120>, or -DCCEX_CUDA=OFF")
    endif()
    find_package(CUDAToolkit REQUIRED)
  else()
    find_package(CUDAToolkit ${CCEX_CUDA_MIN} REQUIRED)
  endif()
  set(CMAKE_CUDA_STANDARD 17)
endif()
message(STATUS "CCEX: MPI ${MPI_CXX_COMPILER}; CUDA backend ${_ccex_cuda}")

# ---- sources: the core in five directories (src/README.md), the bundle clustering in
# clustering/bundle/, the compute backends in calculation/gpu/
set(CCEX_CORE_DIRS input physical clustering clustering/bundle calculation calculation/gpu common auto_convergence)
set(CCEX_CORE_SRCS ${CMAKE_SOURCE_DIR}/src/main.cpp)
foreach(d ${CCEX_CORE_DIRS})
  file(GLOB _srcs CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/${d}/*.cpp)
  list(APPEND CCEX_CORE_SRCS ${_srcs})
endforeach()
file(GLOB CCEX_ADAPTIVE_SRCS CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/clustering/algo_adaptive/*.cpp)
# A30: the batched small-cluster kernel writes out every fused multiply-add it needs; the compiler adds none
set_source_files_properties(${CMAKE_SOURCE_DIR}/src/calculation/cce_batch_kernel.cpp PROPERTIES COMPILE_OPTIONS -ffp-contract=off)
# A42: the small eigensolver's lanes must compute what one matrix alone computes: nothing fused
set_source_files_properties(${CMAKE_SOURCE_DIR}/src/calculation/small_eigensolver.cpp PROPERTIES COMPILE_OPTIONS -ffp-contract=off)

set(CCEX_INCLUDE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/external/eigen ${CMAKE_SOURCE_DIR}/external/uthash/include)
# the warnings the code trips on every compiler; -fpermissive is for the C idioms in src/
set(CCEX_CXX_WARN -fpermissive -Wno-deprecated-declarations -Wno-writable-strings -Wno-unused-result -Wno-format-security)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  list(APPEND CCEX_CXX_WARN -Wno-format-overflow)   # g++'s; clang does not know it and says so once per object
endif()
# -DCCEX_STRICT_WARNINGS=ON (the CI does it, D34): the suppressions above are dropped, the warnings come on, and the
# kinds that are defects -- a printf/scanf format that does not match, a format that is not a literal, an
# unsequenced write, an if that does not guard what its indentation says -- stop the build. The C idioms the
# suppressions are for (string literals as char*, unused results) stay warnings.
option(CCEX_STRICT_WARNINGS "warnings on, and format/sequence-point/misleading-indentation ones as errors" OFF)
if(CCEX_STRICT_WARNINGS)
  set(CCEX_CXX_WARN -fpermissive -Wall -Wformat=2 -Wno-write-strings
      -Werror=format -Werror=format-security -Werror=sequence-point -Werror=misleading-indentation
      -Wno-error=format-truncation)    # snprintf that may cut a message short: safe, so a warning only
endif()
set(CCEX_CPU_OPT -O3 -g)
if(CCEX_NATIVE)
  if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
    list(APPEND CCEX_CPU_OPT -mcpu=native)    # Apple Silicon: the core's instruction set and scheduling (as the mac preset)
  else()
    list(APPEND CCEX_CPU_OPT -march=native)
  endif()
endif()
# icpx defaults to -fp-model=fast; the idm4 preset keeps precise, so the CMake build does too
if(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM|Intel")
  list(APPEND CCEX_CPU_OPT -fp-model=precise)
  set(CCEX_ADAPTIVE_FP -fp-model=precise)
endif()
# the adaptive and SQLite objects get them too (ccex_adaptive below), so that a sanitized main.out
# checks the refinement code it runs in process (D20, 2026-10-02)
set(CCEX_SANITIZE_OPT)
if(CCEX_SANITIZE)
  set(CCEX_SANITIZE_OPT -fsanitize=address,undefined -fno-omit-frame-pointer)
  list(APPEND CCEX_CPU_OPT ${CCEX_SANITIZE_OPT})
  add_link_options(-fsanitize=address,undefined)
endif()

# ---- the refinement policies and SQLite: compiled once, linked into both binaries.
# adaptive_main.cpp comes with its own main() compiled out; what is wanted is
# adaptive_main, so that the binary runs the policies itself. Their flags decide float
# results (the port was verified against CPython under them) and are not shared.
add_library(ccex_adaptive OBJECT ${CCEX_ADAPTIVE_SRCS} ${CMAKE_SOURCE_DIR}/external/sqlite/sqlite3.c)
target_include_directories(ccex_adaptive PRIVATE ${CMAKE_SOURCE_DIR}/include/clustering/algo_adaptive ${CMAKE_SOURCE_DIR}/external/sqlite)
target_compile_definitions(ccex_adaptive PRIVATE CCEX_ADAPTIVE_NO_MAIN)
target_compile_options(ccex_adaptive PRIVATE
  $<$<COMPILE_LANGUAGE:CXX>:-O2 -g -ffp-contract=off -fmath-errno ${CCEX_ADAPTIVE_FP} ${CCEX_SANITIZE_OPT}>
  $<$<COMPILE_LANGUAGE:C>:-O2 -w ${CCEX_SANITIZE_OPT}>)
# the strict warnings reach the policies too (review 20261012-194049-af8e42 R2); warnings only, so the flags above stay
if(CCEX_STRICT_WARNINGS)
  foreach(_w ${CCEX_CXX_WARN})
    target_compile_options(ccex_adaptive PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${_w}>)
  endforeach()
endif()

# ---- the CPU binary: the core but main.cpp as an object library (ccex_core), which the unit tests link
#      as well (cpu-test, and the drivers through unit-build.env below; D26)
set(_ccex_lib_srcs ${CCEX_CORE_SRCS})
list(REMOVE_ITEM _ccex_lib_srcs ${CMAKE_SOURCE_DIR}/src/main.cpp)
add_library(ccex_core OBJECT ${_ccex_lib_srcs})
target_include_directories(ccex_core PRIVATE ${CCEX_INCLUDE})
target_compile_options(ccex_core PRIVATE ${CCEX_CPU_OPT} ${CCEX_CXX_WARN})
target_link_libraries(ccex_core PRIVATE MPI::MPI_CXX)
# A68: the banner's Build line (src/common/utilities.cpp): preset "cmake", and the commit from ccex_build_commit.h, which
# ccex_build_commit writes at every build and changes only when the commit (or "-dirty") does, as do_compile.sh's Makefile:
# git rev-parse --short HEAD, "-dirty" when the tree differs from it, "unknown" outside git
set(_ccex_generated ${CMAKE_BINARY_DIR}/generated)
file(WRITE ${_ccex_generated}/build_commit.cmake [=[
execute_process(COMMAND git rev-parse --short HEAD WORKING_DIRECTORY ${SRC} RESULT_VARIABLE r
                OUTPUT_VARIABLE c OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT r EQUAL 0 OR c STREQUAL "")
  set(c unknown)
else()
  execute_process(COMMAND git diff --quiet HEAD -- WORKING_DIRECTORY ${SRC} RESULT_VARIABLE d OUTPUT_QUIET ERROR_QUIET)
  if(NOT d EQUAL 0)
    set(c "${c}-dirty")
  endif()
endif()
file(CONFIGURE OUTPUT ${OUT} CONTENT "#define CCEX_BUILD_COMMIT \"${c}\"\n")
]=])
set(_ccex_commit_cmd ${CMAKE_COMMAND} -DSRC=${CMAKE_SOURCE_DIR} -DOUT=${_ccex_generated}/ccex_build_commit.h
    -P ${_ccex_generated}/build_commit.cmake)
execute_process(COMMAND ${_ccex_commit_cmd})
add_custom_target(ccex_build_commit COMMAND ${_ccex_commit_cmd} BYPRODUCTS ${_ccex_generated}/ccex_build_commit.h)
add_dependencies(ccex_core ccex_build_commit)
set_source_files_properties(${CMAKE_SOURCE_DIR}/src/common/utilities.cpp PROPERTIES
  INCLUDE_DIRECTORIES ${_ccex_generated} COMPILE_DEFINITIONS "CCEX_BUILD_PRESET=\"cmake\"")
add_executable(ccex ${CMAKE_SOURCE_DIR}/src/main.cpp $<TARGET_OBJECTS:ccex_core> $<TARGET_OBJECTS:ccex_adaptive>)
set_target_properties(ccex PROPERTIES OUTPUT_NAME main.out)
target_include_directories(ccex PRIVATE ${CCEX_INCLUDE})
target_compile_options(ccex PRIVATE ${CCEX_CPU_OPT} ${CCEX_CXX_WARN})
target_link_libraries(ccex PRIVATE MPI::MPI_CXX Threads::Threads ${CMAKE_DL_LIBS})

# How to compile a unit test against this build, for the drivers outside CMake (test/unit_build_env.sh,
# D26): the compiler and the flags, includes and libraries the objects of ccex_core were built with (MPI as
# FindMPI found it), and those objects with the refinement policies'. Written at generation, any generator.
set(_unit_flags ${CMAKE_CXX17_EXTENSION_COMPILE_OPTION} ${CCEX_CPU_OPT} ${CCEX_CXX_WARN} ${MPI_CXX_COMPILE_OPTIONS})
foreach(d ${MPI_CXX_COMPILE_DEFINITIONS})
  list(APPEND _unit_flags -D${d})
endforeach()
set(_unit_includes)
foreach(i ${CCEX_INCLUDE} ${MPI_CXX_INCLUDE_DIRS})
  list(APPEND _unit_includes -I${i})
endforeach()
set(_unit_libs ${MPI_CXX_LINK_FLAGS} ${MPI_CXX_LIBRARIES} ${_ccex_blas_unit_libs} -lpthread -ldl -lm)
if(_ccex_mkl)
  list(APPEND _unit_flags -DCCEX_BLAS_DIRECT)
elseif(_ccex_accelerate)
  list(APPEND _unit_flags -DCCEX_BLAS_DIRECT -DCCEX_BLAS_ACCELERATE)
endif()
if(CCEX_SANITIZE)
  list(APPEND _unit_libs -fsanitize=address,undefined)
endif()
# and the flags of the configuration the objects were built in (CMAKE_CXX_FLAGS_DEBUG etc.), which CMake adds to
# every compile line: -D_GLIBCXX_DEBUG or an alignment flag there changes the layout the objects use (review
# 78d583 R4). A multi-config generator gets one file per configuration, unit-build-<config>.env (R3).
set(_unit_config_flags)
foreach(c ${CMAKE_CONFIGURATION_TYPES} ${CMAKE_BUILD_TYPE})
  string(TOUPPER "${c}" _c)
  list(APPEND _unit_config_flags "$<$<CONFIG:${c}>:${CMAKE_CXX_FLAGS_${_c}}>")
endforeach()
string(JOIN " " _unit_flags_s ${CMAKE_CXX_FLAGS} ${_unit_config_flags} ${_unit_flags})
string(JOIN " " _unit_includes_s ${_unit_includes})
string(JOIN " " _unit_libs_s ${_unit_libs})
# FindMPI quotes the rpath of its link flags; the drivers split on spaces and would keep the quotes
string(REPLACE "\"" "" _unit_libs_s "${_unit_libs_s}")
if(CMAKE_CONFIGURATION_TYPES)
  set(_unit_env_file ${CMAKE_BINARY_DIR}/unit-build-$<CONFIG>.env)
else()
  set(_unit_env_file ${CMAKE_BINARY_DIR}/unit-build.env)
endif()
file(GENERATE OUTPUT ${_unit_env_file} CONTENT
"CCEX_UNIT_CXX='${CMAKE_CXX_COMPILER}'
CCEX_UNIT_ENV=''
CCEX_UNIT_CXXFLAGS='${_unit_flags_s}'
CCEX_UNIT_INCLUDES='${_unit_includes_s}'
CCEX_UNIT_LIBS='${_unit_libs_s}'
CCEX_UNIT_OBJS='$<JOIN:$<TARGET_OBJECTS:ccex_core>, > $<JOIN:$<TARGET_OBJECTS:ccex_adaptive>, >'
")

# ---- the CUDA binary: the same sources with -DCCEX_USE_CUDA, plus the kernels
if(_ccex_cuda)
  add_library(ccex_gpu_backend OBJECT ${CMAKE_SOURCE_DIR}/src/calculation/gpu/gpu_backend.cu
              ${CMAKE_SOURCE_DIR}/src/calculation/gpu/small_cluster_gpu.cu
              ${CMAKE_SOURCE_DIR}/src/calculation/gpu/warp_cluster_gpu.cu)
  target_include_directories(ccex_gpu_backend PRIVATE ${CMAKE_SOURCE_DIR}/include)
  target_compile_options(ccex_gpu_backend PRIVATE -O3 -lineinfo -split-compile=${CCEX_CUDA_SPLIT_COMPILE})
  # A68: the architecture in the Build line; "native" leaves it to the list nvcc compiled for (gpu_backend.cu)
  if(NOT CMAKE_CUDA_ARCHITECTURES STREQUAL "native")
    string(REPLACE ";" "," _ccex_arch "${CMAKE_CUDA_ARCHITECTURES}")
    target_compile_definitions(ccex_gpu_backend PRIVATE CCEX_CUDA_ARCH=${_ccex_arch})
  endif()

  add_executable(ccex_gpu ${CCEX_CORE_SRCS} $<TARGET_OBJECTS:ccex_adaptive> $<TARGET_OBJECTS:ccex_gpu_backend>)
  add_dependencies(ccex_gpu ccex_build_commit)
  set_target_properties(ccex_gpu PROPERTIES OUTPUT_NAME main_gpu.out)
  target_include_directories(ccex_gpu PRIVATE ${CCEX_INCLUDE})
  target_compile_definitions(ccex_gpu PRIVATE CCEX_USE_CUDA)
  target_compile_options(ccex_gpu PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${CCEX_CPU_OPT} ${CCEX_CXX_WARN}>)
  target_link_libraries(ccex_gpu PRIVATE MPI::MPI_CXX Threads::Threads ${CMAKE_DL_LIBS}
                        CUDA::cusolver CUDA::cublas CUDA::cudart)

  # the GPU unit tests need the card: cmake --build build --target gpu-test, inside a job
  set(_gpu_tests test_gcce_cuda test_cce_cuda test_gcce_accuracy test_gcce_refinement test_gcce_batch test_review_failure test_review_memory
                 test_small_cluster_cuda test_warp_cluster_cuda)
  file(GLOB _small_cluster_dumps ${CMAKE_SOURCE_DIR}/test/gpu_verification/check_gpu/data/small_cluster/*.bin)
  file(GLOB _warp_cluster_dumps ${CMAKE_SOURCE_DIR}/test/gpu_verification/check_gpu/data/warp_cluster/*.bin)
  foreach(t ${_gpu_tests})
    add_executable(${t} EXCLUDE_FROM_ALL ${CMAKE_SOURCE_DIR}/test/gpu_verification/check_gpu/${t}.cpp $<TARGET_OBJECTS:ccex_gpu_backend>)
    set_target_properties(${t} PROPERTIES OUTPUT_NAME ${t}.out)
    target_include_directories(${t} PRIVATE ${CCEX_INCLUDE})
    target_compile_definitions(${t} PRIVATE CCEX_USE_CUDA)
    target_compile_options(${t} PRIVATE ${CCEX_CPU_OPT} ${CCEX_CXX_WARN})
    target_link_libraries(${t} PRIVATE MPI::MPI_CXX Threads::Threads ${CMAKE_DL_LIBS} CUDA::cusolver CUDA::cublas CUDA::cudart)
  endforeach()
  add_custom_target(gpu-test
    COMMAND $<TARGET_FILE:test_scaled_product>
    COMMAND $<TARGET_FILE:test_hamiltonian_check>
    COMMAND $<TARGET_FILE:test_gcce_cuda>
    COMMAND $<TARGET_FILE:test_cce_cuda>
    COMMAND $<TARGET_FILE:test_gcce_accuracy>
    COMMAND $<TARGET_FILE:test_gcce_accuracy> ${CMAKE_SOURCE_DIR}/test/gpu_verification/check_gpu/data/gcce_scale_separation_three_spins.txt
    COMMAND $<TARGET_FILE:test_gcce_refinement>
    COMMAND $<TARGET_FILE:test_gcce_batch>
    COMMAND $<TARGET_FILE:test_review_failure>
    COMMAND $<TARGET_FILE:test_review_memory>
    COMMAND $<TARGET_FILE:test_small_cluster_cuda> ${_small_cluster_dumps}
    COMMAND $<TARGET_FILE:test_warp_cluster_cuda> ${_warp_cluster_dumps}
    DEPENDS test_scaled_product test_hamiltonian_check ${_gpu_tests}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif()

# ---- the numeric unit tests (CPU) of the scaled product and of the Hamiltonian check the CUDA
#      cCCE path makes: cmake --build build --target numeric-test
foreach(t test_scaled_product test_hamiltonian_check)
  add_executable(${t} EXCLUDE_FROM_ALL ${CMAKE_SOURCE_DIR}/test/gpu_verification/check_gpu/${t}.cpp)
  set_target_properties(${t} PROPERTIES OUTPUT_NAME ${t}.out)
  target_include_directories(${t} PRIVATE ${CCEX_INCLUDE})
  target_compile_options(${t} PRIVATE ${CCEX_CPU_OPT} ${CCEX_CXX_WARN})
  target_link_libraries(${t} PRIVATE MPI::MPI_CXX)
endforeach()
# the packed kernel's arithmetic on the host (test_small_cluster_cuda without the card): dumps, dump reader, Jacobi scaling
add_executable(test_small_cluster_host EXCLUDE_FROM_ALL ${CMAKE_SOURCE_DIR}/test/gpu_verification/check_gpu/test_small_cluster_cuda.cpp)
set_target_properties(test_small_cluster_host PROPERTIES OUTPUT_NAME test_small_cluster_host.out)
target_include_directories(test_small_cluster_host PRIVATE ${CCEX_INCLUDE})
target_compile_definitions(test_small_cluster_host PRIVATE SMALL_CLUSTER_HOST_ONLY)
target_compile_options(test_small_cluster_host PRIVATE ${CCEX_CPU_OPT} ${CCEX_CXX_WARN})
file(GLOB _small_cluster_host_dumps ${CMAKE_SOURCE_DIR}/test/gpu_verification/check_gpu/data/small_cluster/*.bin)
add_custom_target(numeric-test
  COMMAND $<TARGET_FILE:test_scaled_product>
  COMMAND $<TARGET_FILE:test_hamiltonian_check>
  COMMAND $<TARGET_FILE:test_small_cluster_host> ${_small_cluster_host_dumps}
  DEPENDS test_scaled_product test_hamiltonian_check test_small_cluster_host)

# ---- the CPU regression tests of check_gpu (D17): cmake --build build --target cpu-test.
#      The flags of the CPU binary and -UNDEBUG for their asserts; test_host_lifetimes links the
#      objects of main.out but main.cpp (ccex_core). Seed sharing and the empty-rank
#      ownership run on 1 and 4 ranks (MPIEXEC_EXECUTABLE is the launcher FindMPI found).
set(_cpu_tests test_cluster_log test_mpi_seed test_review_controls test_host_lifetimes)
foreach(t ${_cpu_tests})
  if(t STREQUAL "test_host_lifetimes")
    add_executable(${t} EXCLUDE_FROM_ALL ${CMAKE_SOURCE_DIR}/test/gpu_verification/check_gpu/${t}.cpp
                   $<TARGET_OBJECTS:ccex_core> $<TARGET_OBJECTS:ccex_adaptive>)
  else()
    add_executable(${t} EXCLUDE_FROM_ALL ${CMAKE_SOURCE_DIR}/test/gpu_verification/check_gpu/${t}.cpp)
  endif()
  set_target_properties(${t} PROPERTIES OUTPUT_NAME ${t}.out)
  target_include_directories(${t} PRIVATE ${CCEX_INCLUDE})
  target_compile_options(${t} PRIVATE ${CCEX_CPU_OPT} ${CCEX_CXX_WARN} -UNDEBUG)
  target_link_libraries(${t} PRIVATE MPI::MPI_CXX Threads::Threads ${CMAKE_DL_LIBS})
endforeach()
# the bundle partition, connection and high-order expansion against their own oracles (D31)
set(_bundle_tests test_adaptive_partition test_adaptive_connections test_adaptive_high_order)
foreach(t ${_bundle_tests})
  add_executable(${t} EXCLUDE_FROM_ALL ${CMAKE_SOURCE_DIR}/test/cluster_verification/check_pcce/${t}.cpp
                 ${CMAKE_SOURCE_DIR}/src/clustering/bundle/bundle_partition.cpp ${CMAKE_SOURCE_DIR}/src/clustering/bundle/bundle_expansion.cpp
                 ${CMAKE_SOURCE_DIR}/src/clustering/bundle/bundle_powers.cpp)
  set_target_properties(${t} PROPERTIES OUTPUT_NAME ${t}.out)
  target_include_directories(${t} PRIVATE ${CCEX_INCLUDE})
  target_compile_options(${t} PRIVATE ${CCEX_CPU_OPT} ${CCEX_CXX_WARN} -UNDEBUG)
  target_link_libraries(${t} PRIVATE MPI::MPI_CXX)
endforeach()
# the pCCE partition cache file codec (D32): the test includes pcce_partition_cache.cpp, so it links the other
# objects of ccex_core and writes (and removes) a file in the build folder
add_executable(test_partition_cache_format EXCLUDE_FROM_ALL
  ${CMAKE_SOURCE_DIR}/test/cluster_verification/check_pcce/test_partition_cache_format.cpp
  $<FILTER:$<TARGET_OBJECTS:ccex_core>,EXCLUDE,pcce_partition_cache>)
set_target_properties(test_partition_cache_format PROPERTIES OUTPUT_NAME test_partition_cache_format.out)
target_include_directories(test_partition_cache_format PRIVATE ${CCEX_INCLUDE})
target_compile_options(test_partition_cache_format PRIVATE ${CCEX_CPU_OPT} ${CCEX_CXX_WARN} -UNDEBUG)
target_link_libraries(test_partition_cache_format PRIVATE MPI::MPI_CXX Threads::Threads ${CMAKE_DL_LIBS})
# the fixed-size single-state evolution against the dynamic one, bit for bit (P1-04 of review 20261012-233058-452278):
# the test includes simulator_cce.cpp, so it links the other objects of ccex_core
add_executable(test_fixed_size EXCLUDE_FROM_ALL
  ${CMAKE_SOURCE_DIR}/test/code_verification/check_numerics/test_fixed_size.cpp
  $<FILTER:$<TARGET_OBJECTS:ccex_core>,EXCLUDE,simulator_cce\\.>)
set_target_properties(test_fixed_size PROPERTIES OUTPUT_NAME test_fixed_size.out)
target_include_directories(test_fixed_size PRIVATE ${CCEX_INCLUDE})
target_compile_options(test_fixed_size PRIVATE ${CCEX_CPU_OPT} ${CCEX_CXX_WARN} -UNDEBUG)
target_link_libraries(test_fixed_size PRIVATE MPI::MPI_CXX Threads::Threads ${CMAKE_DL_LIBS})
# A30: the batched small-cluster cCCE kernel against the scalar path, bit for bit, on dumped clusters
add_executable(test_cce_batch EXCLUDE_FROM_ALL
  ${CMAKE_SOURCE_DIR}/test/code_verification/check_numerics/test_cce_batch.cpp $<TARGET_OBJECTS:ccex_core>)
set_target_properties(test_cce_batch PROPERTIES OUTPUT_NAME test_cce_batch.out)
target_include_directories(test_cce_batch PRIVATE ${CCEX_INCLUDE})
target_compile_options(test_cce_batch PRIVATE ${CCEX_CPU_OPT} ${CCEX_CXX_WARN} -UNDEBUG)
target_link_libraries(test_cce_batch PRIVATE MPI::MPI_CXX Threads::Threads ${CMAKE_DL_LIBS})
# A42: the small cCCE eigensolver (fixed-sweep Jacobi) against Eigen and lane against lane
add_executable(test_small_eigensolver EXCLUDE_FROM_ALL
  ${CMAKE_SOURCE_DIR}/test/code_verification/check_numerics/test_small_eigensolver.cpp $<TARGET_OBJECTS:ccex_core>)
set_target_properties(test_small_eigensolver PROPERTIES OUTPUT_NAME test_small_eigensolver.out)
target_include_directories(test_small_eigensolver PRIVATE ${CCEX_INCLUDE})
target_compile_options(test_small_eigensolver PRIVATE ${CCEX_CPU_OPT} ${CCEX_CXX_WARN} -UNDEBUG)
target_link_libraries(test_small_eigensolver PRIVATE MPI::MPI_CXX Threads::Threads ${CMAKE_DL_LIBS})
# A62: the column-block one-state cCCE evolution against the one-time-point path (1e-12 relative)
add_executable(test_cce_columns EXCLUDE_FROM_ALL
  ${CMAKE_SOURCE_DIR}/test/code_verification/check_numerics/test_cce_columns.cpp $<TARGET_OBJECTS:ccex_core>)
set_target_properties(test_cce_columns PROPERTIES OUTPUT_NAME test_cce_columns.out)
target_include_directories(test_cce_columns PRIVATE ${CCEX_INCLUDE})
target_compile_options(test_cce_columns PRIVATE ${CCEX_CPU_OPT} ${CCEX_CXX_WARN} -UNDEBUG)
target_link_libraries(test_cce_columns PRIVATE MPI::MPI_CXX Threads::Threads ${CMAKE_DL_LIBS})
file(GLOB _small_cluster_dumps ${CMAKE_SOURCE_DIR}/test/gpu_verification/check_gpu/data/small_cluster/*.bin)
add_custom_target(cpu-test
  COMMAND $<TARGET_FILE:test_cluster_log>
  COMMAND $<TARGET_FILE:test_review_controls>
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_mpi_seed>
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_mpi_seed>
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_host_lifetimes>
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_host_lifetimes>
  COMMAND $<TARGET_FILE:test_adaptive_partition>
  COMMAND $<TARGET_FILE:test_adaptive_connections>
  COMMAND $<TARGET_FILE:test_adaptive_high_order>
  COMMAND $<TARGET_FILE:test_partition_cache_format> ${CMAKE_BINARY_DIR}/test_partition_cache.bin
  COMMAND $<TARGET_FILE:test_fixed_size>
  COMMAND $<TARGET_FILE:test_cce_batch> ${_small_cluster_dumps}
  COMMAND $<TARGET_FILE:test_small_eigensolver> ${_small_cluster_dumps}
  COMMAND $<TARGET_FILE:test_cce_columns>
  DEPENDS ${_cpu_tests} ${_bundle_tests} test_partition_cache_format test_fixed_size test_cce_batch test_small_eigensolver test_cce_columns)
