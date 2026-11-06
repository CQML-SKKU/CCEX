# gCCE precision fixture

`gcce_scale_separation.txt` contains a 6×6 Hamiltonian, 6×6 initial density
matrix and 3×3 qubit pulse, in column-major complex-double text form. The first
line gives total dimension, qubit dimension and pulse count.

It was extracted from the repository's NV–13C regression input at B=300 G,
using loaded bath site 100 (zero-based), ensemble, and a Hahn echo. The qubit
splitting is approximately 1e7 in the simulator's frequency units while bath
energy differences are much smaller. In the full CCE4 calculation this site's
normalised coherence enters with power -43, amplifying small propagation errors.

The test compares against an extended-precision eigendecomposition of the same
stored Hamiltonian. The fixture is only test data. No simulator branch depends
on this material, site number, matrix size or field.

`gcce_scale_separation_three_spins.txt` stores the analogous 24×24 case for
loaded sites 100, 116 and 134, which retains eigenvector sensitivity after
eigenvalue-only correction. Both fixtures use the original 30-point time grid.

`small_cluster/` and `warp_cluster/` hold the CPU-dumped clusters that `test_small_cluster_cuda` and
`test_warp_cluster_cuda` compare the GPU kernels against. Each folder's `make_dumps.py` writes its files.
