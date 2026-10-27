# Continuous bath drive

What the `pulse.bath_drive` block of the input file does: a continuous drive of bath spins during the qubit's own
sequence, in the cCCE (`method` cCCE, `cluster.mode` standard or pcce), the ensemble approach. It is the setting of
García Oliván et al., arXiv 2512.06948: an NV qubit with instantaneous π pulses (Ramsey, Hahn echo, CPMG) and a P1 bath
driven on some of its hyperfine lines. The keys are listed in
[InputFile_Pulse.md](../user/InputFile_Pulse.md#tag_pulse_bath_drive); the kernel is
`src/calculation/simulator_cce_drive.cpp`, a port of pyCCEX `ccex/simulator/drive.py` (feature/spin-squeezing b3ce43e).

## Frames

Bath spin $k$ sits on a line with lab frequency $\omega_0+d_k$ ($\omega_0=-\gamma_k B$ the electron Zeeman frequency,
$d_k$ the detuning of its Defect axis state). A drive component (one tone, one listed line $l$) has the lab carrier
$\omega_0+l+o$ ($o$ the tone's `offset`). Every spin is put in a frame rotating about $z$ at

$$f_k = l+o \ \text{ when a tone drives its line } (\vert d_k-l\vert<\texttt{line\_tolerance}),\qquad f_k=d_k\ \text{ otherwise},$$

and the rotating-wave approximation keeps what is static there:

$$H_k=(d_k-f_k)J^z_k+W(\cos p\,J^x_k+\sin p\,J^y_k),$$
$$H_{kj}=T_{zz}\big[J^z_kJ^z_j-\tfrac12(J^x_kJ^x_j+J^y_kJ^y_j)\big]\ \ (f_k=f_j),\qquad H_{kj}=T_{zz}J^z_kJ^z_j\ \ (f_k\ne f_j),$$
$$H_{qb}=m_s\sum_k A_{zz,k}J^z_k ,$$

$T$ the point-dipole tensor of the pair and $A$ the qubit coupling of the spin. The frame change $e^{if_ktJ^z_k}$ commutes
with $J^z_k$ and is the same in both qubit branches, so with `qubit_transverse` false it drops out of the coherence
$\mathrm{Tr}[U_\alpha U_\beta^\dagger]/d$; with no tone ($W=0$) the result is the lab-frame cCCE less its non-secular terms.

## The qubit's transverse coupling (`qubit_transverse`, default true)

In branch $m_s$ spin $k$ precesses about $\omega_k\hat z+m_sA_k$ ($A_k=(A_{zx},A_{zy},A_{zz})_k$, $\omega_k=\omega_0+d_k$), with
frequency $W_k$ and axis $n_k$. Each branch gets its own axes ($n_k$, $e_{1,k}$ = $\hat x$ off $n_k$, $e_{2,k}=n_k\times e_{1,k}$),
the frame $e^{-i(\omega_0+f_k)t\,n_k\cdot J_k}$ rotates about $n_k$ and the rotating-wave approximation is taken about these
axes: $H_k$ uses $(d_k-f_k+\delta_k)J^n_k$ with $\delta_k=W_k-\omega_k$ and the part of the field off $n_k$, a pair the secular
part of its tensor in the local axes. The two branches' frames differ, so the propagation applies $R_{\rm new}^\dagger R_{\rm old}$
at every π pulse and $R_{\rm end,\beta}^\dagger R_{\rm end,\alpha}$ before the trace. This holds the Van Vleck shift
$m_s^2\vert A_\perp\vert^2/2\omega_k$ to all orders and the first-order pair terms $m_sA_\perp T_\perp/\omega_k$. `false` is
García Oliván Eq. S.17 ($A_{zz}$ only).

## Lab corrections (`shifts`)

Bloch-Siegert $W_c^2/2(\omega_k+\omega_c)$ from every component $c$, AC Stark $W_c^2/2(\omega_k-\omega_c)$ from every other
component outside $3W_c$, the Van Vleck shift above when `qubit_transverse` is false, and the drive-assisted cross-line
flip-flop: tone $k$ tilts an off-line spin $j$ by $W_k/D_{jk}$ ($D_{jk}=d_j-f_k$), which turns the dropped flip-flop into
$-\tfrac12T_{zz}(W_k/D_{jk})(\cos p\,J^1_k+\sin p\,J^2_k)J^n_j$ (pairs with $\vert D_{jk}\vert\ge3W_k$).

## Propagation and output

CPMG-$N$ with intervals $t/2N, t/N, \dots, t/N, t/2N$ ($N=0$: Ramsey), branch $a$ starting in $\alpha$, $b$ in $\beta$, each
π swapping them; $L(t)=\mathrm{Tr}[U_aU_b^\dagger]/d$, conjugated for even $N\ge2$ as the cCCE kernel reports
$2\langle\beta\vert\rho\vert\alpha\rangle$ for every $N\ge1$ (Ramsey $2\langle\alpha\vert\rho\vert\beta\rangle$). Every cluster goes through
this kernel (no CPU batch, no GPU path); the products, the error correction outputs (`_noDiv`) and MPI are the cCCE's.

## Limits

- A carrier within $3W$ of a line it does not target, or two frames within $3W$ (Hartmann-Hahn), break the per-line
  rotating-wave approximation: reported once, not treated.
- When $W$ is close to the NV-P1 coupling a single cluster's coherence can pass through zero (Rabi-synchronous zeros),
  and the CCE-2 division then needs the error correction.
- Not implemented: sampled states with a mean field (`configuration.s_range`), gCCE, a driven qubit, bath spins other
  than electron-like spin-1/2, the qubit's own transverse terms $A^2/(D\pm\omega_0)$.

## Checks

Against pyCCEX `drive.py` (b3ce43e) on a 20 ppm P1 bath at 4000 G (`garciaolivan2025_repro/v2drive/xval_drive.py`):
the four protocols of García Oliván Fig. 3 (none, two resonant tones, hybrid, four tones) × `qubit_transverse` ×
`shifts` with a Hahn echo, and Ramsey, CPMG-2 and CPMG-3, cCCE order 1 and 2 on 54 and 180 spins: the outputs are the same to
the printed digits (pyCCEX gives NaN where a Ramsey coherence underflows in the CCE-2 division; the stable product here
does not).
