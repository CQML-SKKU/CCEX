#!/bin/sh
# advanced/04_GPU/pCCE2_Hahn_backends: the same run with execution.backend cpu, auto and cuda
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
# cpu with bin/main.out; auto and cuda with bin/main_gpu.out when it exists (bash do_compile.sh gpu on idm4-gpu makes
# both; CCEX_GPU_BIN names another GPU binary). Without it auto runs with main.out (CPU fallback) and cuda is skipped.
CCEX_GPU=${CCEX_GPU_BIN:-$ROOT/bin/main_gpu.out}
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/cpu.log
cp ccein.jsonc results/cpu.jsonc
BACKENDS="cpu auto"; [ -x $CCEX_GPU ] && BACKENDS="cpu auto cuda"
for b in $BACKENDS; do
  [ $b = cpu ] && continue
  BIN=$CCEX; [ -x $CCEX_GPU ] && BIN=$CCEX_GPU
  python3 $ROOT/example/tools/InputKeys/set_input_keys.py ccein.jsonc results/$b.json execution.backend=$b
  mpirun -n $NP $BIN -f results/$b.json -o results > results/$b.log
done
grep -h '\[profile\]' results/*.log | head -20 || true      # timings per backend and dimension
set -- "cpu=results/ccein_wiDiv"; for b in $BACKENDS; do [ $b = cpu ] || set -- "$@" "$b=results/${b}_wiDiv"; done
python3 $PLOT results/coherence.png --unit us "$@"
