# Analyzer

Post-processing of the coherence files CCEX writes: error correction, ensemble
average, fitting and plots. Nothing here runs CCEX.

## What is where

| File | Holds |
|---|---|
| `CoherenceAnalyzer.v1/coherence.py` | the part worth reusing: read and write a coherence file, the error-correction rule, the ensemble average, the file-naming convention. Needs numpy only |
| `CoherenceAnalyzer.v1/DataSpace.py` | the `DataSpace` class the CLI works on: collects files by variable, fits (scipy), draws (matplotlib). Its reading, correction and averaging call `coherence.py` |
| `CoherenceAnalyzer.v1/CoherencePlot_v2.py` | the plot classes |
| `CoherenceAnalyzer.v1/run.py` | the command line (`python3 run.py -h`) |
| `do_*.sh` | one example call each: ensemble average, error correction, single and multi plot, intensity plot, fitting table |
| `data/` | the example results the `do_*.sh` read: P1 1 ppm, CCE2, pulse counts 1 to 128 |

## The file format

One line per time step: the time in ms, a tab, the complex coherence written as
`+re+imj`. Two real columns (`re  im`) are also read. CCEX writes `_noDiv` (the
product of cluster coherences before the sub-cluster division) and `_wiDiv`
(after it); the analyzer writes `_ErrCorr` and `..._EnsAvg_conv<N>`.

## The error-correction rule

CCE divides each cluster's coherence by that of its sub-clusters, and the division
blows up where a sub-cluster coherence passes through zero. So the divided curve
(wiDiv) is checked against the undivided one (noDiv), point by point, with a
tolerance (default 2.0):

1. A point where |wiDiv − noDiv| ≥ tolerance is a divergence.
2. Before the loop, a divergent last point is pulled to noDiv + tolerance − 1e-6,
   so that step 4 always finds a next good point; it then stays at that value.
3. A divergent first point becomes 0.
4. Any other divergent point becomes the mean of the previous corrected value and
   the next non-divergent wiDiv value.

This is the group's convention, not a physical law. Recorded results were made with it, so it must not change; the one
implementation is `coherence.error_correct`.

## Using it from Python

```python
import sys; sys.path.insert(0, 'example/analysis/Analyzer/CoherenceAnalyzer.v1')
from coherence import read_coherence, error_correct, ensemble_average, write_coherence, expand_names

t, L_wD = read_coherence('results/raw/CCE2_NV_conf1_wiDiv')      # t in ms, L complex
_, L_nD = read_coherence('results/raw/CCE2_NV_conf1_noDiv')
L_EC = error_correct(L_wD, L_nD, tolerance=2.0)
write_coherence('results/raw/CCE2_NV_conf1_ErrCorr', t, L_EC)

files = expand_names('results/raw/', 'CCE2_NV_confCONF_ErrCorr', [('CONF', range(1, 21))])
t, L_avg, n = ensemble_average(read_coherence(p) for p in files.values())
```

## Using the command line

Two typical calls (`test/code_verification/check_input_forms/errcorr.sh` uses the first):

```
run.py -d DIR -fi PATTERN_wiDiv -err_d DIR -err_fi PATTERN_noDiv -vn CONF -v 1 2 3 \
       -uc ms -err -err_tol 2.0 -ntn -fo DIR/PATTERN_ErrCorr        # error correction
run.py -d DIR -fi PATTERN_ErrCorr -vn CONF -v 1 2 3 -uc ms -ea -ntn -fo DIR/NAME_EnsAvg   # ensemble average
```

`-vn` names the variables, `-v` their values (one `-v` group per variable), and
each variable name in the pattern is replaced by each value. Plots need a display
or `MPLBACKEND=Agg`.

## File names in `data/`

`wD` is the divided result (the final one), `nD` the undivided one used to correct
it; `EC` is error-corrected; `conv<N>` an average over N samples.
