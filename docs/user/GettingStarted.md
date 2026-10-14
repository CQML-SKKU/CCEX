# Getting Started
<a id="GettingStarted"></a>
  
This is the getting started guide.

<a id="HowToRun"></a>
## How to run

To run the code using <tt> mpirun </tt>, use the following command :  

``` 
$ mpirun -n process_number code_path -f input_file 
```  

- process_number : The number of processes to run.  
- code_path : The path to your executable code.  
- input_file : The input file required for the code to run. <br>
  (input file should be <tt>.json</tt> file)  

**Example** 

``` 
$ mpirun -n 4 /usr/local/bin/mycode -f data/input.json 
```

**Possble** **command-line** **options**
- -h, --help : Print the options and exit  
- -f : [Input file](InputFile.md) option (Required) (See [Input file](InputFile.md))  
- -m : Method (See [method](InputFile_General.md#tag_method))  
- -q : Quantity (See [quantity](InputFile_General.md#tag_quantity))   
- -I : A bath file: the first `-I` replaces the input's `bathfile` list, each further `-I` adds one (See [bathfile](InputFile_General.md#tag_bathfile))  
- -s : File path to set the spin state of bath spins (See [bath_state_file](InputFile_General.md#tag_bath_state_file))   
- -a : File path to set principal axies of the defects in the bath (See [Defect.axis_file](InputFile_Defect.md#tag_df_axis_file))  
- -S : File path to set the spin state of on-site spins of the defects in the bath (See [Defect.state_file](InputFile_Defect.md#tag_df_state_file))  
- -N : The number of state configuration to do time ensemble average (See [nstate](InputFile_General.md#tag_nstate))  
- -B : Magnetic field in z direction (Unit : G) (See [bfield](InputFile_General.md#tag_bfield))  
- -o : Result file path (See [outfile](InputFile_Output.md#tag_outfile))  
- -v : Verbose output (longer log)  
- -P : Write the cluster plan to `<outfile>.clusters.json` and stop before the dynamics (cCCE, gCCE and cCCE with `cluster_mode` pcce and explicit connections)  

Three run modes take the input file and replace the ordinary calculation:
- --auto_conv -f FILE : The automatic convergence; the same as `-f FILE` with an `auto_convergence` block (See [Automatic convergence tags](InputFile_AutoConvergence.md)). It takes no other option.  
- --auto_status -f FILE : Where that automatic convergence is (running, interrupted, stopped or how it ended); nothing is computed.  
- --adaptive ... : The adaptive refinement; every word after `--adaptive` goes to it (e.g. `--adaptive --input in.json --output DIR ...`). An `adaptive` block in the input asks for the same; giving both is refused (See [the `adaptive` block](../methods/adaptive-input-block.md)).  

**About** **options**

The code read options twice :
- Read options from input file called as In-file options
- Read optoins from command line called as Command-line options <br>
> **Note:** Command-line options are "priority"

<a id="EnvironmentVariables"></a>
**Environment** **variables**

None is needed for an ordinary run. They switch an optimisation off for a comparison, point the automatic convergence
at its tool, or serve the checks; set them in front of the command (`CCEX_SMALL_EIG=0 mpirun -n 32 main.out -f ccein.json`)
or export them in a job script.

| variable | when it is read | effect |
|---|---|---|
| `CCEX_CPU_BATCH=0` | every run | the CPU batch of small cCCE clusters (bath dimension 2 and 4, A30) off: each cluster alone. The batch gives the same bits as that path; where it would not (another compiler or its options), it checks itself at start and stays off by itself, with a line in the log |
| `CCEX_SMALL_EIG=0` | every run | the fixed-sweep Jacobi eigensolver of the same small clusters (A42) off: Eigen as before 2026-10-13, the output exactly as then. Matrices outside a safe range or not converged already go to Eigen one by one |
| `CUDA_VISIBLE_DEVICES` | GPU build (`main_gpu.out`) | which card the run may use; empty or not matching any card hides them all, and `backend` `auto` then stays on the CPU |
| `CCEX_PYTHON` | `--auto_conv`, `--auto_status` | the Python that runs the automatic-convergence tool (default `python3`; it needs numpy, scipy, matplotlib) |
| `CCEX_AUTOCONV_DIR` | `--auto_conv`, `--auto_status` | the tool's folder when it is not found next to the binary (`example/tools/AutoConvergence`); a folder others can write to is refused |
| `CCEX_CHECK_DISORDER=1` | every run with `nstate` > 0 | the checks only: also builds the static mean field of the bath on one rank the old way and stops on any bit that differs from the split one (A13) |
| `CCEX_SMALL_CLUSTER_DUMP=FILE`, `CCEX_SMALL_CLUSTER_DUMP_N=N` | cCCE runs | the checks only: write the first N (400) small clusters with the CPU result to FILE (the CPU batch is off meanwhile) |
| `CCEX_WARP_CLUSTER_DUMP=FILE`, `CCEX_WARP_CLUSTER_DUMP_N=N` | gCCE runs | the checks only: the same for the clusters of the GPU warp kernel |
| `CCEX_TEST_GPU_OOM_BYTES=N` | GPU build | the checks only: every GPU allocation from N bytes fails as if the card were full, to exercise the fall-back to the CPU |

At build time `do_compile.sh` reads `CCEX_JOBS` (parallel compile jobs, default 5), `CCEX_BLAS` (`auto` by default: Intel MKL
for the large clusters when it is found, `off` to leave it out, `mkl` to require it), `CXX`, `CXXFLAGS_EXTRA`,
`LIBRARY_EXTRA`, `CUDA_HOME` and `CUDA_ARCH` for the GPU binary, and `CCEX_SANITIZE=1` (address and undefined-behaviour
checks); the header of `do_compile.sh` describes each.

<hr>
<a id="Workflow"></a>
## Workflow

The CCEX software follow the below work process. 

When you give the arguments in the command line, the code frist read the arguments. <br>
Note that if you don't give the input file accesible with -f option, then the code doesn't work. <br>

Code read the all tags in input file and update the ...  

<img src="CCEX_workflow.png" alt="" width="600">
