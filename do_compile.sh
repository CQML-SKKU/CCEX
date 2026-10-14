#!/usr/bin/env bash
# do_compile.sh <preset>  -- writes the Makefile for this machine, then runs make.
#
# The build rules are one template (the end of this file). A preset is only what
# differs from machine to machine: the compiler, its flags, what to load, and
# whether there is CUDA.
#
#   generic  any machine with an MPI C++ wrapper on PATH. Taken from the environment:
#            CXX (default mpicxx), CC (default cc), CXXFLAGS (default below),
#            CXXFLAGS_EXTRA, LIBRARY_EXTRA, and for the GPU binary CUDA_HOME
#            (or nvcc on PATH), CUDA_ARCH (compute capability without the dot,
#            e.g. 86) and NVCC. Examples:
#                bash do_compile.sh generic
#                CXX=mpic++ bash do_compile.sh generic
#                CUDA_HOME=/usr/local/cuda CUDA_ARCH=86 bash do_compile.sh generic
#   wsl-gnu  WSL/Ubuntu with Intel oneAPI: Intel MPI's mpicxx (g++ behind it); "wsl" is accepted
#            as the old name
#   wsl-icpx the same with icpx behind Intel MPI's mpicxx (needs the oneAPI C++ compiler:
#            apt install intel-oneapi-compiler-dpcpp-cpp); closer to idm's icpx build and faster
#   idm4     idm: module intel-23.2, mpiicpc with icpx
#   gpu      idm4-gpu: OpenMPI + CUDA 12.8, the C++ with icpx when oneAPI is there (else gcc); builds the GPU binary too
#   nurion   KISTI Nurion: intel/19.0.5 + impi/19.0.5
#   mac      macOS with an MPI wrapper (Homebrew open-mpi or mpich)
#
# Every preset builds bin/main.out (CPU). A preset with CUDA also builds
# bin/main_gpu.out, the same sources with -DCCEX_USE_CUDA and the CUDA kernels,
# for "backend": "cuda" / "auto" in the input. Eigen, uthash and sqlite are in
# external/; nothing else is needed. Intel MKL is used when it is found (CCEX_BLAS
# below; every preset decides after its modules/setvars.sh, CCEX_BLAS=off turns it off).
# The old presets' EIGEN_USE_MKL_ALL was defined after Eigen was included and had no
# effect, and no MKL symbol was linked (2026-09-24: ldd, and the 25-case numerics check
# gave the same numbers without it); it was removed then.
#
#   CCEX_JOBS=N       make -jN (default 5: idm's login node is shared)
#   CCEX_SANITIZE=1   -fsanitize=address,undefined on the CPU objects, the adaptive ones and SQLite
#   CCEX_BLAS=auto|mkl|off   large cCCE/gCCE clusters by Intel MKL (zgemm, zheevd): auto (default) when
#                     MKL is found, mkl insists (stops without it), off never; see below

if [ $# -ne 1 ]; then
    echo "usage: bash do_compile.sh <preset>" >&2
    echo "  preset: generic / wsl-gnu (old name wsl) / wsl-icpx / idm4 / gpu / nurion / mac" >&2
    exit 1
fi
PRESET="$1"

# what every preset may set; the template reads these
GPU=0
MAKE_EXPORTS=""                       # a line the Makefile needs, e.g. which compiler an MPI wrapper calls
CC="${CC:-cc}"
NVCC="${NVCC:-}"
CUDA_ARCH="${CUDA_ARCH:-}"
LIBRARY_EXTRA="${LIBRARY_EXTRA:-}"
NOTE=""
# the warnings the code trips on every compiler; -fpermissive is for the C idioms in src/
WARN="-fpermissive -Wno-deprecated-declarations -Wno-writable-strings -Wno-unused-result -Wno-format-overflow -Wno-format-security"

load_modules() {
    # the group's servers select compilers with environment modules; a machine without
    # the module command, or without these modules, keeps what is on PATH
    if ! type module >/dev/null 2>&1 && [ -r /etc/profile.d/modules.sh ]; then
        . /etc/profile.d/modules.sh
    fi
    if type module >/dev/null 2>&1; then
        module purge
        module load "$@"
    fi
}

case "$PRESET" in

generic)
    CXX="${CXX:-mpicxx}"
    CXXFLAGS="${CXXFLAGS:--std=gnu++17 -O3 -g -march=native $WARN}"
    if [ -z "${CUDA_HOME:-}" ] && command -v nvcc >/dev/null 2>&1; then
        CUDA_HOME="$(dirname "$(dirname "$(command -v nvcc)")")"
    fi
    if [ -n "${CUDA_HOME:-}" ]; then
        GPU=1
        if [ -z "$CUDA_ARCH" ]; then
            echo "[generic] CUDA_HOME is set: also building bin/main_gpu.out. Set CUDA_ARCH to the" >&2
            echo "[generic] card's compute capability without the dot (e.g. 86); assuming 80." >&2
            CUDA_ARCH=80
        fi
    fi
    NOTE="generic: CXX=$CXX"
    ;;

wsl-gnu|wsl)
    # Intel oneAPI: setvars.sh puts Intel MPI's mpicxx on PATH. hash -r, because a
    # /usr/bin/mpicxx (Ubuntu's MPICH) seen earlier in the shell would shadow it; the
    # absolute path is baked into the Makefile for the same reason. I_MPI_FABRICS=shm
    # keeps single-node MPI off the WSL network namespace (needed with
    # networkingMode=mirrored).
    . /opt/intel/oneapi/setvars.sh > /dev/null 2>&1
    hash -r
    export I_MPI_FABRICS=shm
    CXX="${I_MPI_ROOT}/bin/mpicxx"
    if [ ! -x "$CXX" ]; then
        echo "[wsl] Intel MPI's mpicxx not found at $CXX; is oneAPI installed and setvars.sh working?" >&2
        exit 1
    fi
    CXXFLAGS="-std=gnu++17 -O3 -g -march=native -mtune=native -D_FORTIFY_SOURCE=0 $WARN"
    NOTE="wsl-gnu: g++ behind Intel MPI $I_MPI_ROOT, I_MPI_FABRICS=shm"
    ;;

wsl-icpx)
    # As wsl, with icpx behind Intel MPI's mpicxx (I_MPI_CXX, baked into the Makefile with its
    # full path so a later make need not source setvars.sh). icpx's default is -fp-model=fast,
    # which folds isinf/isnan; precise is set here for the core and the template keeps it for
    # the adaptive objects (2026-09-30). The oneAPI version differs from idm's 2023.2, so the
    # numbers are close to idm's but not bit for bit: check_numerics' head reference was recorded
    # with the wsl-gnu (g++) build, and icpx 2026.1 differs from it by up to 3.5e-9 in noDiv in three
    # P1 cases (2026-10-07), so run the standing checks with the wsl-gnu preset. --force: a shell that sourced setvars.sh
    # before the compiler was installed would otherwise skip it and leave icpx off PATH.
    . /opt/intel/oneapi/setvars.sh --force > /dev/null 2>&1
    hash -r
    export I_MPI_FABRICS=shm
    CXX="${I_MPI_ROOT}/bin/mpicxx"
    if ! command -v icpx >/dev/null 2>&1; then
        echo "[wsl-icpx] icpx not found: apt install intel-oneapi-compiler-dpcpp-cpp, or use the wsl-gnu preset" >&2
        exit 1
    fi
    MAKE_EXPORTS="export I_MPI_CXX := $(command -v icpx)"
    CC="$(command -v icx)"
    CXXFLAGS="-std=gnu++17 -O3 -g -march=native -fp-model=precise $WARN"
    NOTE="wsl-icpx: icpx $(icpx --version | head -1 | grep -o '[0-9][0-9.]*' | head -1) behind Intel MPI $I_MPI_ROOT, I_MPI_FABRICS=shm"
    ;;

idm4)
    load_modules intel-23.2/icc-23.2
    CXX=mpiicpc
    # icpx, not icpc: mpiicpc defaults to icpc (Intel C++ Classic, deprecated in 2023),
    # which is 3.2x slower on this Eigen-template-heavy code (19.5 s against 6.1 s on
    # the P1 10 ppm gCCE regression case). icpx is in the same module.
    MAKE_EXPORTS="export I_MPI_CXX := icpx"
    CXXFLAGS="-std=c++17 -O3 -g -DNDEBUG -fp-model=precise -xCORE-AVX512 -Wno-deprecated-declarations -Wno-writable-strings"
    NOTE="idm4: mpiicpc with icpx"
    ;;

gpu)
    # idm4-gpu: OpenMPI and CUDA 12.8 from modules; the C++ with icpx when oneAPI is
    # installed (/opt/intel/oneapi; there is no module for it), which is 15-28% faster
    # than gcc 11 on the CPU path (2026-09-24, 2 ranks: p1_cce4 26.6 s -> 19.1 s), and
    # gcc otherwise. The .cu keeps nvcc's g++ host: nvcc 12.8 takes icpx as host compiler
    # only from oneAPI 2024 on. RTX PRO 6000 Blackwell is compute capability 12.0. The
    # login node is the compute node, so keep CCEX_JOBS small while jobs are running.
    load_modules gcc/gcc-11.5.0 gcc/openmpi-4.1.5 cuda/cuda-12.8
    CXX="$(command -v mpicxx)"          # OpenMPI's wrapper, found before setvars.sh puts Intel MPI's first
    CC=gcc
    CXX_NAME=gcc
    if [ -r /opt/intel/oneapi/setvars.sh ]; then
        . /opt/intel/oneapi/setvars.sh > /dev/null 2>&1
        if command -v icpx >/dev/null 2>&1; then
            MAKE_EXPORTS="export OMPI_CXX := $(command -v icpx)"   # the full path: a later make need not source setvars.sh
            CC="$(command -v icx)"
            CXX_NAME=icpx
        fi
    fi
    CUDA_HOME="${CUDA_HOME:-/usr/local/cuda-12.8}"
    CUDA_ARCH="${CUDA_ARCH:-120}"
    GPU=1
    CXXFLAGS="-std=gnu++17 -O3 -g -march=native $WARN"
    [ "$CXX_NAME" = icpx ] && CXXFLAGS="$CXXFLAGS -fp-model=precise"
    NOTE="gpu: $CXX_NAME + OpenMPI + CUDA $CUDA_HOME, arch $CUDA_ARCH"
    ;;

nurion)
    load_modules craype-network-opa craype-mic-knl intel/19.0.5 impi/19.0.5
    CXX=mpiicpc
    # C++17: the code has needed it since the bundle clustering (2026-09); the preset
    # said -std=c++11 until 2026-09-24, which no longer compiles.
    CXXFLAGS="-std=c++17 -O3 -g -diag-disable=2196 -diag-disable=10441 -Wno-deprecated-declarations"
    NOTE="nurion: mpiicpc, intel/19.0.5"
    ;;

mac)
    CXX="${CXX:-mpicxx}"
    CXXFLAGS="-std=gnu++17 -O2 -g $WARN"
    NOTE="mac: $CXX"
    ;;

*)
    # going on would run make with whatever Makefile the last preset left behind
    echo "unknown preset: $PRESET" >&2
    echo "  preset: generic / wsl-gnu (old name wsl) / wsl-icpx / idm4 / gpu / nurion / mac" >&2
    exit 1
    ;;
esac

# CCEX_BLAS (A46): the cCCE ensemble Hahn step's D^3 product (the recursive path) from bath dimension 21
# and the cluster eigensolve from dimension 63 (cCCE: the bath dimension; gCCE: the whole qubit x bath
# dimension) by MKL's zgemm/zheevd (include/common/blas_option.h).
# idm job 17337: h-BN r 0.72-0.88 where the time is in such clusters, P1 untouched; the curves move
# in the last digits (h-BN gCCE late times up to 9e-8). auto (the default; user decision 2026-10-13):
# on when $MKLROOT (set by setvars.sh or an intel module, so decided here, after the preset) holds
# the three libraries; mkl: the same, but stop when they are not there; off: never. MKL is linked
# sequential (one thread per MPI rank); the binary then needs MKL's libraries at run time (rpath).
BLAS_MODE="$(echo "${CCEX_BLAS:-auto}" | tr '[:upper:]' '[:lower:]')"
MKL_LIB=""
MKL_WHY=""
case "$BLAS_MODE" in
off) MKL_WHY="CCEX_BLAS=off" ;;
auto|mkl)
    if [ -z "${MKLROOT:-}" ]; then
        MKL_WHY="MKLROOT is not set"
    else
        for d in "$MKLROOT/lib/intel64" "$MKLROOT/lib"; do
            if [ -e "$d/libmkl_core.so" ] && [ -e "$d/libmkl_intel_lp64.so" ] && [ -e "$d/libmkl_sequential.so" ]; then MKL_LIB="$d"; break; fi
        done
        [ -n "$MKL_LIB" ] || MKL_WHY="no libmkl_core/intel_lp64/sequential.so under $MKLROOT/lib(/intel64)"
    fi
    if [ "$BLAS_MODE" = mkl ] && [ -z "$MKL_LIB" ]; then
        echo "[build] CCEX_BLAS=mkl: $MKL_WHY (source setvars.sh or load the intel module)" >&2; exit 1
    fi
    ;;
*) echo "[build] CCEX_BLAS=$CCEX_BLAS: expected auto, mkl or off" >&2; exit 1 ;;
esac
if [ -n "$MKL_LIB" ]; then
    LIBRARY_EXTRA="$LIBRARY_EXTRA -L$MKL_LIB -Wl,-rpath,$MKL_LIB -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
    CXXFLAGS="$CXXFLAGS -DCCEX_BLAS_DIRECT"
    NOTE="$NOTE, MKL on ($MKL_LIB)"
    echo "[build] MKL on ($MKL_LIB)"
else
    NOTE="$NOTE, MKL off ($MKL_WHY)"
    echo "[build] MKL off ($MKL_WHY)"
fi
CXXFLAGS="$CXXFLAGS ${CXXFLAGS_EXTRA:-}"
# The refinement policies get their own flags (below), not CXXFLAGS, and they carry
# infinities and NaNs as values (an unmeasured error is inf). A compiler whose default
# assumes neither folds isnan()/isinf() to false: icpx 2023.2 at -O2 does both, icpc
# 2023.2 isnan (measured 2026-09-30); on idm4-gpu a merge run crashed writing an
# infinity into summary.json. So the compiler the Makefile will call is asked, with a
# small program built with the policies' own flags; if it gets either wrong,
# -fp-model=precise is added and asked again, and the build stops if that does not help.
adaptive_fp_probe() {                   # $1: extra flags; prints "11" when isinf and isnan hold
    probe_dir=$(mktemp -d) || return 1
    cat > "$probe_dir/fp.cpp" <<'PROBE'
#include <cmath>
#include <cstdio>
#include <cstdlib>
int main(int argc, char** argv) {
    volatile double big = std::atof(argc > 1 ? argv[1] : "1e308");
    double x = big*10.0, y = x-x;
    std::printf("%d%d", (int)std::isinf(x), (int)std::isnan(y));
    return 0;
}
PROBE
    probe_env=""
    case "$MAKE_EXPORTS" in export\ *) probe_env="$(echo "${MAKE_EXPORTS#export }" | sed 's/ *:= */=/')" ;; esac
    if env $probe_env $CXX -std=gnu++17 -O2 -ffp-contract=off -fmath-errno $1 \
           "$probe_dir/fp.cpp" -o "$probe_dir/fp" > /dev/null 2>&1; then
        "$probe_dir/fp" 2> /dev/null
    else
        echo "nobuild"
    fi
    rm -f "$probe_dir/fp.cpp" "$probe_dir/fp"; rmdir "$probe_dir"
}
# "export I_MPI_CXX := icpx" -> "I_MPI_CXX=icpx": what the unit tests compiled outside make need too (unit-env)
UNIT_ENV=""
case "$MAKE_EXPORTS" in export\ *) UNIT_ENV="$(echo "${MAKE_EXPORTS#export }" | sed 's/ *:= */=/')" ;; esac
ADAPTIVE_FP=""
fp_answer="$(adaptive_fp_probe "")"
# a probe that does not even build says nothing about isinf/isnan: the wrapper or its libraries are broken
if [ "$fp_answer" = "nobuild" ]; then
    echo "[build] $CXX cannot build a small C++ program (the MPI wrapper, its compiler or its libraries);" >&2
    echo "        try: $CXX -std=gnu++17 hello.cpp" >&2
    exit 1
fi
if [ "$fp_answer" != "11" ]; then
    ADAPTIVE_FP="-fp-model=precise"
    fp_answer="$(adaptive_fp_probe "$ADAPTIVE_FP")"
    if [ "$fp_answer" != "11" ]; then
        echo "[build] the compiler behind $CXX drops isinf/isnan even with -fp-model=precise" >&2
        echo "        (probe answered '$fp_answer', expected '11'); src/clustering/algo_adaptive needs both" >&2
        exit 1
    fi
    echo "[build] src/clustering/algo_adaptive: $ADAPTIVE_FP (the compiler's default drops isinf/isnan)"
fi
# the sanitizer flags go on the adaptive and SQLite objects too (their own rules below), so that a
# sanitized main.out checks the refinement code it runs in process (D20, 2026-10-02)
SANFLAGS=""
if [ "${CCEX_SANITIZE:-0}" = "1" ]; then
    SANFLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"
    CXXFLAGS="$CXXFLAGS $SANFLAGS"
fi
if [ "$GPU" = "1" ]; then
    NVCC="${NVCC:-$CUDA_HOME/bin/nvcc}"
    TARGETS='$(TARGET) $(TARGET_GPU)'
else
    TARGETS='$(TARGET)'
fi

# ----------------------------------------------------------------------------
# The template. $VAR is filled in now; \$(VAR) stays for make.
# ----------------------------------------------------------------------------
cat << EOF > ./Makefile
# Generated by "bash do_compile.sh $PRESET" on $(date -u +"%Y-%m-%dT%H:%M:%SZ"); do not edit, rerun the script.
# $NOTE
CXX := $CXX
CC  := $CC
$MAKE_EXPORTS
CXXFLAGS := $CXXFLAGS
NVCC      := $NVCC
CUDA_HOME := ${CUDA_HOME:-}
CUDA_ARCH := $CUDA_ARCH

SRC_DIR := ./src
OBJ_DIR := ./obj
OBJ_GPU_DIR := ./obj_gpu
BIN_DIR := ./bin
TARGET     := \$(BIN_DIR)/main.out
TARGET_GPU := \$(BIN_DIR)/main_gpu.out

# The core in five directories (src/README.md); the bundle clustering is
# clustering/bundle/, the compute backends are calculation/gpu/.
CORE_DIRS := \$(SRC_DIR)/input \$(SRC_DIR)/physical \$(SRC_DIR)/clustering \$(SRC_DIR)/clustering/bundle \$(SRC_DIR)/calculation \$(SRC_DIR)/calculation/gpu \$(SRC_DIR)/common \$(SRC_DIR)/auto_convergence
NAMES := \$(notdir \$(basename \$(wildcard \$(SRC_DIR)/*.cpp) \$(foreach d,\$(CORE_DIRS),\$(wildcard \$(d)/*.cpp))))
OBJS     := \$(addprefix \$(OBJ_DIR)/,\$(addsuffix .o,\$(NAMES)))
OBJS_GPU := \$(addprefix \$(OBJ_GPU_DIR)/,\$(addsuffix .o,\$(NAMES))) \$(OBJ_GPU_DIR)/gpu_backend.o \$(OBJ_GPU_DIR)/small_cluster_gpu.o \$(OBJ_GPU_DIR)/warp_cluster_gpu.o

# The refinement policies (src/clustering/algo_adaptive) and sqlite do not depend on the backend:
# compiled once, linked into both binaries. adaptive_main.cpp comes with its own main()
# compiled out; what is wanted is adaptive_main, so that this binary runs the policies
# itself. Their flags decide float results (the port was verified against CPython under
# them) and are not shared with src/.
ADAPTIVE_DIR  := ./src/clustering/algo_adaptive
ADAPTIVE_OBJS := \$(patsubst \$(ADAPTIVE_DIR)/%.cpp,\$(OBJ_DIR)/adaptive_%.o,\$(wildcard \$(ADAPTIVE_DIR)/*.cpp)) \$(OBJ_DIR)/sqlite3.o
ADAPTIVE_FPFLAGS := -ffp-contract=off -fmath-errno $ADAPTIVE_FP
SANFLAGS := $SANFLAGS
ADAPTIVE_INCLUDE := -I./include/clustering/algo_adaptive -I./external/sqlite

INCLUDE := -I./external/eigen -I./external/uthash/include/ -I./include/
LIBRARY := -lpthread -ldl -lm $LIBRARY_EXTRA
CUDA_LIBDIR := \$(firstword \$(wildcard \$(CUDA_HOME)/lib64 \$(CUDA_HOME)/lib))
LIBRARY_GPU := -L\$(CUDA_LIBDIR) -Wl,-rpath,\$(CUDA_LIBDIR) -lcusolver -lcublas -lcudart

# Header dependency tracking: -MMD writes obj/<name>.d with every header the unit
# included, -MP adds an empty target for each so that deleting a header does not stop
# the build. Without it make compares the .o only against its .cpp.
DEPFLAGS := -MMD -MP

all: $TARGETS

\$(TARGET): \$(OBJS) \$(ADAPTIVE_OBJS) | \$(BIN_DIR)
	\$(CXX) \$^ -o \$@ \$(CXXFLAGS) \$(LIBRARY)

\$(TARGET_GPU): \$(OBJS_GPU) \$(ADAPTIVE_OBJS) | \$(BIN_DIR)
	\$(CXX) \$^ -o \$@ \$(CXXFLAGS) \$(LIBRARY) \$(LIBRARY_GPU)

# The same sources twice: \$(1) is where the objects go, \$(2) what that build adds.
define COMPILE_RULES
\$(1)/%.o: \$(SRC_DIR)/%.cpp | \$(1)
	\$\$(CXX) -c \$\$< -o \$\$@ \$\$(CXXFLAGS) \$(2) \$\$(DEPFLAGS) \$\$(INCLUDE)
\$(1)/%.o: \$(SRC_DIR)/input/%.cpp | \$(1)
	\$\$(CXX) -c \$\$< -o \$\$@ \$\$(CXXFLAGS) \$(2) \$\$(DEPFLAGS) \$\$(INCLUDE)
\$(1)/%.o: \$(SRC_DIR)/physical/%.cpp | \$(1)
	\$\$(CXX) -c \$\$< -o \$\$@ \$\$(CXXFLAGS) \$(2) \$\$(DEPFLAGS) \$\$(INCLUDE)
\$(1)/%.o: \$(SRC_DIR)/clustering/%.cpp | \$(1)
	\$\$(CXX) -c \$\$< -o \$\$@ \$\$(CXXFLAGS) \$(2) \$\$(DEPFLAGS) \$\$(INCLUDE)
\$(1)/%.o: \$(SRC_DIR)/clustering/bundle/%.cpp | \$(1)
	\$\$(CXX) -c \$\$< -o \$\$@ \$\$(CXXFLAGS) \$(2) \$\$(DEPFLAGS) \$\$(INCLUDE)
\$(1)/%.o: \$(SRC_DIR)/calculation/%.cpp | \$(1)
	\$\$(CXX) -c \$\$< -o \$\$@ \$\$(CXXFLAGS) \$(2) \$\$(DEPFLAGS) \$\$(INCLUDE)
\$(1)/%.o: \$(SRC_DIR)/calculation/gpu/%.cpp | \$(1)
	\$\$(CXX) -c \$\$< -o \$\$@ \$\$(CXXFLAGS) \$(2) \$\$(DEPFLAGS) \$\$(INCLUDE)
\$(1)/%.o: \$(SRC_DIR)/common/%.cpp | \$(1)
	\$\$(CXX) -c \$\$< -o \$\$@ \$\$(CXXFLAGS) \$(2) \$\$(DEPFLAGS) \$\$(INCLUDE)
\$(1)/%.o: \$(SRC_DIR)/auto_convergence/%.cpp | \$(1)
	\$\$(CXX) -c \$\$< -o \$\$@ \$\$(CXXFLAGS) \$(2) \$\$(DEPFLAGS) \$\$(INCLUDE)
\$(1):
	mkdir -p \$(1)
endef
\$(eval \$(call COMPILE_RULES,\$(OBJ_DIR),))
\$(eval \$(call COMPILE_RULES,\$(OBJ_GPU_DIR),-DCCEX_USE_CUDA))
# A30: the batched small-cluster kernel writes out every fused multiply-add it needs; the compiler adds none
\$(OBJ_DIR)/cce_batch_kernel.o \$(OBJ_GPU_DIR)/cce_batch_kernel.o: CXXFLAGS += -ffp-contract=off
# A42: the small eigensolver's lanes must compute what one matrix alone computes: nothing fused
\$(OBJ_DIR)/small_eigensolver.o \$(OBJ_GPU_DIR)/small_eigensolver.o: CXXFLAGS += -ffp-contract=off

\$(OBJ_GPU_DIR)/gpu_backend.o: \$(SRC_DIR)/calculation/gpu/gpu_backend.cu \$(wildcard ./include/calculation/gpu/*.h) | \$(OBJ_GPU_DIR)
	\$(NVCC) -std=c++17 -O3 -lineinfo -I./include \\
	    -gencode arch=compute_\$(CUDA_ARCH),code=sm_\$(CUDA_ARCH) \\
	    -gencode arch=compute_\$(CUDA_ARCH),code=compute_\$(CUDA_ARCH) -c \$< -o \$@
# A26: the batched small-cluster cCCE kernel (one thread per cluster)
\$(OBJ_GPU_DIR)/small_cluster_gpu.o: \$(SRC_DIR)/calculation/gpu/small_cluster_gpu.cu \$(wildcard ./include/calculation/gpu/*.h) | \$(OBJ_GPU_DIR)
	\$(NVCC) -std=c++17 -O3 -lineinfo -I./include \\
	    -gencode arch=compute_\$(CUDA_ARCH),code=sm_\$(CUDA_ARCH) \\
	    -gencode arch=compute_\$(CUDA_ARCH),code=compute_\$(CUDA_ARCH) -c \$< -o \$@

# A40: the warp kernel (one cluster per warp; about 3.5 minutes of nvcc, the longest object of the build)
\$(OBJ_GPU_DIR)/warp_cluster_gpu.o: \$(SRC_DIR)/calculation/gpu/warp_cluster_gpu.cu \$(wildcard ./include/calculation/gpu/*.h) | \$(OBJ_GPU_DIR)
	\$(NVCC) -std=c++17 -O3 -lineinfo -I./include \\
	    -gencode arch=compute_\$(CUDA_ARCH),code=sm_\$(CUDA_ARCH) \\
	    -gencode arch=compute_\$(CUDA_ARCH),code=compute_\$(CUDA_ARCH) -c \$< -o \$@

\$(OBJ_DIR)/adaptive_%.o: \$(ADAPTIVE_DIR)/%.cpp | \$(OBJ_DIR)
	\$(CXX) -c \$< -o \$@ -std=gnu++17 -O2 -g -DCCEX_ADAPTIVE_NO_MAIN \$(ADAPTIVE_FPFLAGS) \$(SANFLAGS) \$(DEPFLAGS) \$(ADAPTIVE_INCLUDE)

\$(OBJ_DIR)/sqlite3.o: ./external/sqlite/sqlite3.c | \$(OBJ_DIR)
	\$(CC) -c \$< -o \$@ -O2 -w \$(SANFLAGS) -I./external/sqlite

\$(BIN_DIR):
	mkdir -p \$(BIN_DIR)

# The GPU unit tests need the card: run "make gpu-test" inside a job. test_review_memory
# includes cuda_runtime.h itself, hence the CUDA include folder.
GPU_TESTS := test_gcce_cuda test_cce_cuda test_gcce_accuracy test_gcce_refinement test_gcce_batch test_review_failure test_review_memory test_small_cluster_cuda test_warp_cluster_cuda
\$(BIN_DIR)/test_%.out: test/gpu_verification/check_gpu/test_%.cpp \$(OBJ_GPU_DIR)/gpu_backend.o | \$(BIN_DIR)
	\$(CXX) \$(INCLUDE) -I\$(CUDA_HOME)/include -DCCEX_USE_CUDA \$(CXXFLAGS) \$< \$(OBJ_GPU_DIR)/gpu_backend.o \$(LIBRARY) \$(LIBRARY_GPU) -o \$@
\$(BIN_DIR)/test_small_cluster_cuda.out: test/gpu_verification/check_gpu/test_small_cluster_cuda.cpp \$(OBJ_GPU_DIR)/small_cluster_gpu.o \$(OBJ_GPU_DIR)/warp_cluster_gpu.o | \$(BIN_DIR)
	\$(CXX) \$(INCLUDE) -I\$(CUDA_HOME)/include \$(CXXFLAGS) \$< \$(OBJ_GPU_DIR)/small_cluster_gpu.o \$(OBJ_GPU_DIR)/warp_cluster_gpu.o \$(LIBRARY) \$(LIBRARY_GPU) -o \$@
\$(BIN_DIR)/test_warp_cluster_cuda.out: test/gpu_verification/check_gpu/test_warp_cluster_cuda.cpp \$(OBJ_GPU_DIR)/small_cluster_gpu.o \$(OBJ_GPU_DIR)/warp_cluster_gpu.o | \$(BIN_DIR)
	\$(CXX) \$(INCLUDE) -I\$(CUDA_HOME)/include \$(CXXFLAGS) \$< \$(OBJ_GPU_DIR)/small_cluster_gpu.o \$(OBJ_GPU_DIR)/warp_cluster_gpu.o \$(LIBRARY) \$(LIBRARY_GPU) -o \$@
\$(BIN_DIR)/test_scaled_product.out: test/gpu_verification/check_gpu/test_scaled_product.cpp include/calculation/scaled_product.h | \$(BIN_DIR)
	\$(CXX) \$(INCLUDE) \$(CXXFLAGS) \$< -o \$@
\$(BIN_DIR)/test_hamiltonian_check.out: test/gpu_verification/check_gpu/test_hamiltonian_check.cpp include/calculation/hamiltonian_check.h | \$(BIN_DIR)
	\$(CXX) \$(INCLUDE) \$(CXXFLAGS) \$< -o \$@
# the packed kernel's arithmetic on the host (test_small_cluster_cuda without the card): dumps, dump reader, Jacobi scaling
\$(BIN_DIR)/test_small_cluster_host.out: test/gpu_verification/check_gpu/test_small_cluster_cuda.cpp include/calculation/gpu/small_cluster_math.h | \$(BIN_DIR)
	\$(CXX) \$(INCLUDE) \$(CXXFLAGS) -DSMALL_CLUSTER_HOST_ONLY \$< -o \$@
# The CPU regression tests of check_gpu (D17): the build's own flags, -UNDEBUG for their asserts;
# test_host_lifetimes links the objects of bin/main.out but main.o. MPIRUN=... picks the launcher.
CPU_TESTS := test_cluster_log test_mpi_seed test_review_controls test_host_lifetimes
# The bundle partition, connection and high-order expansion against their own oracles (D31): the two bundle
# sources compiled into each test with the build's flags, -UNDEBUG for their asserts.
BUNDLE_TESTS := test_adaptive_partition test_adaptive_connections test_adaptive_high_order
# The pCCE partition cache file codec (D32): the test includes pcce_partition_cache.cpp itself, so it links the other
# objects of main.out; it is given a file to write under obj/ and removes it.
CACHE_TEST_OBJS := \$(filter-out \$(OBJ_DIR)/main.o \$(OBJ_DIR)/pcce_partition_cache.o,\$(OBJS))
\$(BIN_DIR)/test_partition_cache_format.out: test/cluster_verification/check_pcce/test_partition_cache_format.cpp \$(CACHE_TEST_OBJS) | \$(BIN_DIR) \$(OBJ_DIR)
	\$(CXX) \$(INCLUDE) \$(CXXFLAGS) -UNDEBUG -MMD -MP -MF \$(OBJ_DIR)/test_partition_cache_format.test.d -MT \$@ \$< \$(CACHE_TEST_OBJS) \$(LIBRARY) -o \$@
# The fixed-size single-state evolution against the dynamic one, bit for bit (P1-04 of review 20261012-233058-452278):
# the test includes simulator_cce.cpp itself, so it links the other objects of main.out.
FIXED_TEST_OBJS := \$(filter-out \$(OBJ_DIR)/main.o \$(OBJ_DIR)/simulator_cce.o,\$(OBJS))
\$(BIN_DIR)/test_fixed_size.out: test/code_verification/check_numerics/test_fixed_size.cpp \$(FIXED_TEST_OBJS) | \$(BIN_DIR) \$(OBJ_DIR)
	\$(CXX) \$(INCLUDE) \$(CXXFLAGS) -UNDEBUG -MMD -MP -MF \$(OBJ_DIR)/test_fixed_size.test.d -MT \$@ \$< \$(FIXED_TEST_OBJS) \$(LIBRARY) -o \$@
# A30: the batched small-cluster cCCE kernel against the scalar path, bit for bit, on the dumped clusters of
# check_gpu/data/small_cluster; it links the objects of main.out but main.o.
\$(BIN_DIR)/test_cce_batch.out: test/code_verification/check_numerics/test_cce_batch.cpp \$(filter-out \$(OBJ_DIR)/main.o,\$(OBJS)) | \$(BIN_DIR) \$(OBJ_DIR)
	\$(CXX) \$(INCLUDE) \$(CXXFLAGS) -UNDEBUG -MMD -MP -MF \$(OBJ_DIR)/test_cce_batch.test.d -MT \$@ \$< \$(filter-out \$(OBJ_DIR)/main.o,\$(OBJS)) \$(LIBRARY) -o \$@
# A42: the small cCCE eigensolver (fixed-sweep Jacobi) against Eigen and lane against lane
\$(BIN_DIR)/test_small_eigensolver.out: test/code_verification/check_numerics/test_small_eigensolver.cpp \$(filter-out \$(OBJ_DIR)/main.o,\$(OBJS)) | \$(BIN_DIR) \$(OBJ_DIR)
	\$(CXX) \$(INCLUDE) \$(CXXFLAGS) -UNDEBUG -MMD -MP -MF \$(OBJ_DIR)/test_small_eigensolver.test.d -MT \$@ \$< \$(filter-out \$(OBJ_DIR)/main.o,\$(OBJS)) \$(LIBRARY) -o \$@
BUNDLE_SRCS := src/clustering/bundle/bundle_partition.cpp src/clustering/bundle/bundle_expansion.cpp src/clustering/bundle/bundle_powers.cpp
\$(addprefix \$(BIN_DIR)/,\$(addsuffix .out,\$(BUNDLE_TESTS))): \$(BIN_DIR)/%.out: test/cluster_verification/check_pcce/%.cpp \$(BUNDLE_SRCS) | \$(BIN_DIR) \$(OBJ_DIR)
	\$(CXX) \$(INCLUDE) \$(CXXFLAGS) -UNDEBUG -MMD -MP -MF \$(OBJ_DIR)/\$*.test.d -MT \$@ \$< \$(BUNDLE_SRCS) -o \$@
MPIRUN ?= mpirun
# What they test is mostly inline in headers, so each writes obj/<test>.test.d (included below with the
# other .d files) and is rebuilt when one of its headers changes.
\$(BIN_DIR)/test_cluster_log.out \$(BIN_DIR)/test_mpi_seed.out \$(BIN_DIR)/test_review_controls.out: \$(BIN_DIR)/%.out: test/gpu_verification/check_gpu/%.cpp | \$(BIN_DIR) \$(OBJ_DIR)
	\$(CXX) \$(INCLUDE) \$(CXXFLAGS) -UNDEBUG -MMD -MP -MF \$(OBJ_DIR)/\$*.test.d -MT \$@ \$< -o \$@
\$(BIN_DIR)/test_host_lifetimes.out: test/gpu_verification/check_gpu/test_host_lifetimes.cpp \$(filter-out \$(OBJ_DIR)/main.o,\$(OBJS)) \$(ADAPTIVE_OBJS) | \$(BIN_DIR) \$(OBJ_DIR)
	\$(CXX) \$(INCLUDE) \$(CXXFLAGS) -UNDEBUG -MMD -MP -MF \$(OBJ_DIR)/test_host_lifetimes.test.d -MT \$@ \$< \$(filter-out \$(OBJ_DIR)/main.o,\$(OBJS)) \$(ADAPTIVE_OBJS) \$(LIBRARY) -o \$@

# How to compile a unit test against this build, for the drivers outside make (test/unit_build_env.sh, D26):
# the compiler, the environment its wrapper needs, the flags, includes and libraries of bin/main.out, and its
# objects but main.o, all with absolute paths. main.out comes first, so that what is printed is up to date.
unit-env: \$(TARGET)
	@echo "CCEX_UNIT_CXX='\$(CXX)'"
	@echo "CCEX_UNIT_ENV='$UNIT_ENV'"
	@echo "CCEX_UNIT_CXXFLAGS='\$(CXXFLAGS)'"
	@echo "CCEX_UNIT_INCLUDES='\$(patsubst -I./%,-I\$(CURDIR)/%,\$(INCLUDE))'"
	@echo "CCEX_UNIT_LIBS='\$(LIBRARY)'"
	@echo "CCEX_UNIT_OBJS='\$(abspath \$(filter-out \$(OBJ_DIR)/main.o,\$(OBJS)) \$(ADAPTIVE_OBJS))'"

.PHONY: all clean info gpu-test numeric-test cpu-test unit-env
numeric-test: \$(BIN_DIR)/test_scaled_product.out \$(BIN_DIR)/test_hamiltonian_check.out \$(BIN_DIR)/test_small_cluster_host.out
	./\$(BIN_DIR)/test_scaled_product.out
	./\$(BIN_DIR)/test_hamiltonian_check.out
	./\$(BIN_DIR)/test_small_cluster_host.out test/gpu_verification/check_gpu/data/small_cluster/*.bin
cpu-test: \$(addprefix \$(BIN_DIR)/,\$(addsuffix .out,\$(CPU_TESTS) \$(BUNDLE_TESTS) test_partition_cache_format test_fixed_size test_cce_batch test_small_eigensolver))
	./\$(BIN_DIR)/test_cluster_log.out
	./\$(BIN_DIR)/test_review_controls.out
	\$(MPIRUN) -n 1 ./\$(BIN_DIR)/test_mpi_seed.out
	\$(MPIRUN) -n 4 ./\$(BIN_DIR)/test_mpi_seed.out
	\$(MPIRUN) -n 1 ./\$(BIN_DIR)/test_host_lifetimes.out
	\$(MPIRUN) -n 4 ./\$(BIN_DIR)/test_host_lifetimes.out
	./\$(BIN_DIR)/test_adaptive_partition.out
	./\$(BIN_DIR)/test_adaptive_connections.out
	./\$(BIN_DIR)/test_adaptive_high_order.out
	./\$(BIN_DIR)/test_partition_cache_format.out \$(OBJ_DIR)/test_partition_cache.bin
	./\$(BIN_DIR)/test_fixed_size.out
	./\$(BIN_DIR)/test_cce_batch.out test/gpu_verification/check_gpu/data/small_cluster/*.bin
	./\$(BIN_DIR)/test_small_eigensolver.out test/gpu_verification/check_gpu/data/small_cluster/*.bin
gpu-test: numeric-test \$(addprefix \$(BIN_DIR)/,\$(addsuffix .out,\$(GPU_TESTS)))
	./\$(BIN_DIR)/test_gcce_cuda.out
	./\$(BIN_DIR)/test_cce_cuda.out
	./\$(BIN_DIR)/test_gcce_accuracy.out
	./\$(BIN_DIR)/test_gcce_accuracy.out test/gpu_verification/check_gpu/data/gcce_scale_separation_three_spins.txt
	./\$(BIN_DIR)/test_gcce_refinement.out
	./\$(BIN_DIR)/test_gcce_batch.out
	./\$(BIN_DIR)/test_review_failure.out
	./\$(BIN_DIR)/test_review_memory.out
	./\$(BIN_DIR)/test_small_cluster_cuda.out test/gpu_verification/check_gpu/data/small_cluster/*.bin
	./\$(BIN_DIR)/test_warp_cluster_cuda.out test/gpu_verification/check_gpu/data/warp_cluster/*.bin

info:
	@echo "CXX      = \$(CXX)"
	@echo "CXXFLAGS = \$(CXXFLAGS)"
	@echo "targets  = $TARGETS"
	@\$(CXX) -show 2>/dev/null | head -1 || true

clean:
	rm -rf \$(OBJ_DIR) \$(OBJ_GPU_DIR) \$(BIN_DIR)

-include \$(wildcard \$(OBJ_DIR)/*.d \$(OBJ_GPU_DIR)/*.d)
EOF

# A changed Makefile has to invalidate every object: make compares each .o against its
# .cpp and headers, never against the flags it was built with, so switching preset or
# compiler would otherwise relink the old objects and hand back a binary that looks
# new and is not.
# The stamp is the Makefile without its first line, which holds the generation time (with it, every
# rerun with the same preset cleaned and rebuilt everything; D23, 2026-10-02), followed by the list of
# sources: the Makefile names them by wildcard, so a deleted source would leave its object and the
# old binary in place, and make would not relink.
build_fingerprint() {
    tail -n +2 Makefile
    find src external/sqlite \( -name '*.cpp' -o -name '*.c' -o -name '*.cu' \) | LC_ALL=C sort
}
if ! cmp -s <(build_fingerprint) .ccex_makefile.stamp 2>/dev/null; then
    echo "[build] build settings or sources changed -- cleaning objects first"
    make clean >/dev/null 2>&1 || true
fi
build_fingerprint > .ccex_makefile.stamp

echo "[build] $NOTE"
make -j"${CCEX_JOBS:-5}"
