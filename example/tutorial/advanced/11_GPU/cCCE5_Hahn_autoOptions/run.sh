#!/bin/sh
# advanced/11_GPU/cCCE5_Hahn_autoOptions: backend auto without gpu_options, the values bin/main_gpu.out chooses ([gpu-auto])
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
CCEX_GPU=${CCEX_GPU_BIN:-$ROOT/bin/main_gpu.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-32}
mkdir -p results
python3 $ROOT/example/tools/InputKeys/set_input_keys.py ccein.jsonc results/cpu.json execution.backend=cpu
mpirun -n $NP $CCEX -f results/cpu.json -o results > results/cpu.log
if [ -x $CCEX_GPU ]; then
  mpirun -n $NP $CCEX_GPU -f ccein.jsonc -o results > results/ccein.log
  grep -h '\[gpu-auto\]' results/ccein.log || echo "no [gpu-auto] line: no GPU visible, or every key given"
else
  echo "no $CCEX_GPU: auto runs on the CPU (bash do_compile.sh gpu makes it)"
  mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
fi
python3 $PLOT results/coherence.png --unit us "cpu=results/cpu_wiDiv" "auto=results/ccein_wiDiv"
