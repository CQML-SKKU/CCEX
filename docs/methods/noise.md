# Gaussian phase noise (method noise)

`"method": "noise"` treats the bath as a source of classical noise on the qubit's transition frequency and gives the
coherence in the Gaussian approximation, the semi-classical CCE: a sum over the bath's modes, the second cumulant, or
the noise spectrum. For the correlation function it follows PyCCE's `CCENoise`. The keys are in
[InputFile_General](../user/InputFile_General.md#tag_noise_spectrum); the tutorial is
`example/tutorial/advanced/09_noise/`.

## The coherence

Under a sequence of ideal π pulses the qubit picks up the phase $\varphi = \int f(t)\,\delta\omega(t)\,dt$, where
$f(t) = \pm 1$ changes sign at each pulse. For Gaussian noise the coherence is

$$W(T) = e^{-\chi(T)}, \qquad \chi = \langle \varphi^2 \rangle / 2 .$$

This is exact only for Gaussian noise. It leaves out the qubit's back-action on the bath and the non-Gaussian part of
the noise, so it suits many weakly coupled noise sources. Strongly coupled spins need cCCE or gCCE.

## The noise from the bath

With `hamiltonian.noise.spectrum` `"bath"` (the default) every cCCE cluster contributes

$$\chi_C = \frac{1}{2D} \sum_{ab} |\langle a|\delta\omega|b\rangle|^2\, |Y(\omega_{ab})|^2, \qquad
Y(\omega) = \int f(t)\, e^{i\omega t}\, dt ,$$

over the eigenstates $|a\rangle$ of $H_\text{bath} + p \sum_i \mathbf{A}_{z,i}\cdot\mathbf{I}_i$ with
$p = (|m_\alpha| + |m_\beta|)/2$, at infinite temperature. $D$ is the cluster's dimension.

- The noise operator is $\delta\omega = (m_\alpha - m_\beta) \sum_i A_{zz,i} I_{z,i}$.
- `pseudo_secular` true adds $A_{zx,i} I_{x,i} + A_{zy,i} I_{y,i}$.

The $\chi_C$ are summed over the CCE terms. This equals the cCCE product of the $e^{-\chi_C}$ and avoids the underflow
of a single factor. A sampled bath state (`configuration.s_range`) enters only
through the mean field of the other bath spins.

The two choices of the noise operator can differ by up to 0.04 in $W$ for a ¹³C bath with DFT hyperfine near an NV
(500 G, Hahn echo) and by 1e-3 for a P1 bath. Which is the better approximation is open.

## Model spectra

Without a bath, `spectrum` can be one of two models:

- Ornstein–Uhlenbeck noise, $\langle\delta\nu(t)\delta\nu(0)\rangle = \sigma^2 e^{-|t|/\tau_c}$. $\chi$ is exact for
  any ideal π sequence.
- A power law, $S(f) = A\,(f/1\,\text{kHz})^{-\beta}$ on $[f_\text{min}, f_\text{max}]$, integrated by
  Gauss–Legendre.

Such a run reads no bath, and `order`, `bfield`, `rbath`, `rdip` and `Qubit` may be left out.

## Correlation function and spectrum

`hamiltonian.noise.compute` chooses what is computed from the bath:

- `"coherence"` (default): $\chi$ straight from each cluster's eigenmodes, with no time grid.
- `"spectrum"`: the correlation $C(t) = \operatorname{Re}\operatorname{Tr}[\rho\{B(t), B\}]/2$ on the grid of
  `hamiltonian.noise.fft`, its two-sided spectrum $S(f)$ by FFT, and $\chi$ from $S(f)$. By default the grid ends at the
  run's last time and has 8 points per period of the fastest bath mode kept.
- `"correlation"`: $C(t)$ alone at the run's times. With `s_range` the cluster starts in its sampled pure state, as in
  PyCCE, so this is then a different quantity from the other two.

Without `s_range` the three agree. On the ¹³C bath of `basic/01_first_run` the coherence from `"spectrum"` is within
2.4e-9 of `"coherence"` (`example/tutorial/advanced/09_noise/noiseSpectrum_13C/`).

## Where it runs

- Supported: `cluster.mode` standard, one qubit, ideal π pulses (FID, Hahn, CPMG, Manual).
- With `execution.backend` `cuda`, the sums over the modes run on the GPU and the eigensystems on the CPU. `auto` and the
  model spectra run on the CPU.
- Refused: the outputs that read single clusters (`output.savemode` cluster_coherence and cluster_density_matrix,
  `rdip_ladder`, `output.cluster_contributions`).

## Checks

`test/code_verification/check_noise/` compares the model spectra with the OU closed forms and scipy quadrature, the
three `compute` values with each other, `"correlation"` with PyCCE, and P1 and ¹³C baths with an independent implementation of the mode sum, the cumulant and
the spectrum (FID, Hahn and CPMG; with and without the mean field).

Reference for the correlation function: PyCCE's `CCENoise` (M. Onizhuk and G. Galli, Adv. Theory Simul. 4, 2100254
(2021)).
