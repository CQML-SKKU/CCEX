# make_natab_baths.py N OUTDIR SIZE: the MPI half of make_natab_baths.sh (run it as mpirun -n NP python3 make_natab_baths.py).
# The configurations are shared over the ranks (rank r makes every k with (k - 1) % size == r), and each one is made by
# BathGenerator.py on MPI.COMM_SELF with the random seed k, so the files are the ones a single process makes.
import os, sys, random, runpy, traceback
import numpy as np
from mpi4py import MPI

n, out, size = int(sys.argv[1]), sys.argv[2], sys.argv[3]
root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
gen = os.path.join(root, "example/tutorial/basic/01_first_run/bath/genbath")
script = os.path.join(root, "example/tools/BathGenerator/BathGenerator.py")
comm = MPI.COMM_WORLD
failed = []
for k in range(1 + comm.rank, n + 1, comm.size):
    cfg = f"{out}/config_{k}"
    if os.path.isfile(cfg) and os.path.getsize(cfg) > 0:
        continue                                               # files already there are kept
    stdout, stderr = sys.stdout, sys.stderr
    with open(f"{out}/process_{k}", "w") as log:
        sys.stdout = sys.stderr = log
        np.random.seed(k); random.seed(k)
        sys.argv = [script, "-o", cfg, "-i", f"{gen}/POSCAR_pure", "-I", f"{gen}/POSCAR_NV", "-p", f"{gen}/prob_natab",
                    "--size", size, "--Defect", "V", "--Defect_center", "0.5", "0.5", "0.5"]
        try:
            runpy.run_path(script, run_name="__main__", init_globals={"bathgen_comm": MPI.COMM_SELF})
        except BaseException:                                  # sys.exit() of an input error, or a crash: the file decides
            traceback.print_exc()
        sys.stdout.flush()
    sys.stdout, sys.stderr = stdout, stderr
    if os.path.isfile(cfg) and os.path.getsize(cfg) > 0:
        os.replace(f"{cfg}_defect", f"{out}/config_defect")    # the same NV position for every k
    else:
        failed.append(k)
        print(f"bath {k} was not made: see {out}/process_{k}", flush=True)
nfailed = comm.allreduce(len(failed))
# numpy + mpi4py can abort or hang at interpreter exit ("double free"); finalize and leave without Python's teardown
sys.stdout.flush(); sys.stderr.flush()
MPI.Finalize()
os._exit(1 if nfailed else 0)
