# Semiclassical noise (method noise)

`"method": "noise"` treats the bath as a source of classical noise on the qubit's transition frequency and gives the
coherence in the Gaussian approximation, the semi-classical CCE. It computes it from a sum over the bath's modes, from
the second cumulant or from the noise spectrum, and for the correlation function it follows PyCCE's `CCENoise`. The
keys are in [General tags](../user/InputFile_General.md#tag_noise_spectrum), the tutorial is
`example/tutorial/advanced/12_noise/`, and `test/code_verification/check_noise/` compares the method with closed forms,
quadrature, PyCCE and an independent implementation.

1. [The coherence](#noise_coherence): the Gaussian phase and its limits.
2. [The noise from the bath](#noise_bath): the sum over each cluster's modes.
3. [Model spectra](#noise_models): Ornstein–Uhlenbeck and a power law, without a bath.
4. [Correlation function and spectrum](#noise_compute): the three values of `compute`.
5. [Where it runs](#noise_support): what is supported and refused.

<a id="noise_coherence"></a>
## 1. The coherence

Under a sequence of ideal π pulses the qubit picks up the phase $\varphi = \int f(t)\,\delta\omega(t)\,dt$, where
$f(t) = \pm 1$ changes sign at each pulse. For Gaussian noise the coherence is

$$W(T) = e^{-\chi(T)}, \qquad \chi = \langle \varphi^2 \rangle / 2 .$$

This is exact only for Gaussian noise. It leaves out the qubit's back-action on the bath and the non-Gaussian part of
the noise, so it suits many weakly coupled noise sources. Strongly coupled spins need cCCE or gCCE.

<a id="noise_bath"></a>
## 2. The noise from the bath

With `hamiltonian.noise.spectrum` `"bath"` (the default) every cCCE cluster contributes

$$\chi_C = \frac{1}{2D} \sum_{ab} |\langle a|\delta\omega|b\rangle|^2\, |Y(\omega_{ab})|^2, \qquad
Y(\omega) = \int f(t)\, e^{i\omega t}\, dt ,$$

over the eigenstates $|a\rangle$ of $H_\text{bath} + p \sum_i \mathbf{A}_{z,i}\cdot\mathbf{I}_i$ with
$p = (|m_\alpha| + |m_\beta|)/2$, at infinite temperature. $D$ is the cluster's dimension.

- The noise operator is $\delta\omega = (m_\alpha - m_\beta) \sum_i A_{zz,i} I_{z,i}$.
- `pseudo_secular` true adds $A_{zx,i} I_{x,i} + A_{zy,i} I_{y,i}$. The two choices can differ in a bath with nuclei
  close to the qubit.

The $\chi_C$ are summed over the CCE terms. This equals the cCCE product of the $e^{-\chi_C}$ and avoids the underflow
of a single factor. A sampled bath state (`configuration.s_range`) enters only through the mean field of the other bath
spins.

<a id="noise_models"></a>
## 3. Model spectra

Without a bath, `spectrum` can be one of two models:

- Ornstein–Uhlenbeck noise, $\langle\delta\nu(t)\delta\nu(0)\rangle = \sigma^2 e^{-|t|/\tau_c}$. $\chi$ is exact for
  any ideal π sequence.
- A power law, $S(f) = A\,(f/1\,\text{kHz})^{-\beta}$ on $[f_\text{min}, f_\text{max}]$, integrated by
  Gauss–Legendre.

Such a run reads no bath, and `order`, `bfield`, `rbath`, `rdip` and `Qubit` may be left out.

<a id="noise_compute"></a>
## 4. Correlation function and spectrum

`hamiltonian.noise.compute` chooses what is computed from the bath:

- `"coherence"` (default) gives $\chi$ straight from each cluster's eigenmodes, with no time grid.
- `"spectrum"` gives the correlation $C(t) = \operatorname{Re}\operatorname{Tr}[\rho\lbrace B(t), B\rbrace ]/2$ on the grid of
  `hamiltonian.noise.fft`, its two-sided spectrum $S(f)$ by FFT, and $\chi$ from $S(f)$. By default the grid ends at the
  run's last time and has 8 points per period of the fastest bath mode kept.
- `"correlation"` gives $C(t)$ alone at the run's times. With `s_range` the cluster starts in its sampled pure state,
  as in PyCCE, so this is then a different quantity from the other two.

Without `s_range` the three agree. `example/tutorial/advanced/12_noise/noiseSpectrum_13C/` computes the spectrum of
the ¹³C bath of `basic/01_first_run`.

<a id="noise_support"></a>
## 5. Where it runs

- Supported: `cluster.mode` standard, one qubit, ideal π pulses (FID, Hahn, CPMG, Manual).
- With `execution.backend` `cuda`, the sums over the modes run on the GPU and the eigensystems on the CPU. `auto` and the
  model spectra run on the CPU.
- Refused: the outputs that read single clusters (`output.savemode` cluster_coherence and cluster_density_matrix,
  `rdip_ladder`, `output.cluster_contributions`).

<a id="noise_references"></a>
## References

- M. Onizhuk and G. Galli, PyCCE: A Python package for cluster correlation expansion simulations of spin qubit
  dynamics, Adv. Theory Simul. 4, 2100254 (2021), doi:10.1002/adts.202100254 (the correlation function of `CCENoise`).
