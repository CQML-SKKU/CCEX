# example/

This folder has three parts. `tutorial/` holds worked examples to run, `analysis/` holds the analyzer for the output,
and `tools/` holds the programs that prepare the input.

## Start here

```
bash do_compile.sh wsl-gnu                               # build bin/main.out (other presets: see the top README)
cd example/tutorial/basic/01_first_run && NP=4 bash run.sh    # about 1 s; writes results/ and plots results/coherence.png
```

Compare `results/coherence.png` with `expected.png` in the same folder (T2 1.084 ms). Then work through `basic/`
in order and pick from `advanced/`. [tutorial/README.md](tutorial/README.md) lists them. Each topic's README says what each calculation shows,
which keys matter, how long it runs and what to try next.

To start an input of your own, run `bin/main.out -e`. It writes `ccex_example.jsonc`, the first-run input with every
other key commented out under its explanation and default. It runs once its `bathfile` names a bath file.

To find the method and the tutorial that fit your system, see
[What should I use for my case?](../docs/user/Examples.md#WhatToUse).

## analysis/ and tools/

| Folder | Contents |
|---|---|
| `analysis/Analyzer/` | `CoherenceAnalyzer.v1/run.py` reads the `_noDiv` and `_wiDiv` files, corrects division errors, averages over configurations, fits and plots. `coherence.py` offers the same as a library. |
| `tools/AutoConvergence/` | The automatic convergence of the time window, rdip, rbath and the number of configurations. It runs when an input has an `auto_convergence` block, and its README lists the keys and outputs. |
| `tools/BathGenerator/` | `BathGenerator.py` writes a bath file from a POSCAR unit cell. |
| `tools/SpinTensorReader/` | `main.py` reads hyperfine and quadrupole tensors from VASP or QE output. |
| `tools/IsotopeTable/` | `make_isotope_table.py` writes the built-in isotope table (`src/input/isotope_table.cpp`). |
| `tools/StateFiles/` | `make_state_files.py` writes the bath state, axis and defect state files of an input before the run. |
| `tools/InputKeys/` | `convert_input_keys.py` rewrites inputs with older keys. `set_input_keys.py` writes a copy of an input with some keys changed. |

[Tools](../docs/user/Tools.md) and [Analysis](../docs/user/Analysis.md) describe them in detail.
