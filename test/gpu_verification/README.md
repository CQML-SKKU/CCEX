# test/gpu_verification/

Checks of the CUDA backend (`src/calculation/gpu/`, `backend` in the input). Three folders.

| folder | what | needs |
|---|---|---|
| `check_gpu/` | 14 `test_*.cpp` unit tests (which `make` target runs which: `check_gpu/README.md`) and nine `check_*.py` that run CPU and CUDA builds on the same inputs and compare them (pulses, gCCE4 accuracy and batching, pCCE, seeds and sampled states, options, host leaks). `examples/` and `data/` hold their inputs, the seven `validation-*.json` the recorded results. | a CUDA build (`bin/main_gpu.out`, idm4-gpu) for most. Four run on CPU builds: `check_review_controls.py --cpu`, `check_seed_states.py --cpu`, `check_host_leaks.py` without `--cuda` (`--binary` and `--baseline` are executables), `check_review_baseline.py` (`--baseline` an executable, `--results` an earlier results folder). |
| `gpu_benchmarks/` | timing, not correctness: CPU against CUDA, rank scaling, the PL6 pCCE case | idm4-gpu |
| `server/` | the three standing Slurm jobs for idm4-gpu (the `gpu` partition of idm's Slurm): build, CCE2 measurement, high orders, with the environment they use; two measurement jobs (`min_dim`, `gpu_sharing`). Last runs: `server/README.md`. | idm4-gpu |
