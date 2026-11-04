# Dispersive materials and subpixel (CMT-like) interfaces

`stack=<file.json>` describes a planar multilayer stack (normal along y) of dispersive materials. Layer boundaries
may fall anywhere relative to the mesh. `smode=s` (default) averages the materials over each Yee cell, like
Lumerical's conformal mesh (CMT). `smode=c` is the staircase control. Without `stack=` the solver is unchanged
(gate D0-1, gate 0-1a, gate P0-1).

```bash
python3 materials/fit.py                                   # data -> Lorentz fits (materials/models.json)
python3 -c "import disp_tmm as D; D.write_stack('runs/arc.json', 0.6, [('Si3N4', 1.003, 1.128), ('Si', 1.128, 1000)])"
./fdtd3d_oblique bloch=1 Lx=0.1 Lz=0.1 kxn=0.2 nl=40 npml=40 sf=40 tf=120 pol=s src=g fc=1.1 fw=0.2 \
    nfreq=31 fmin=0.75 fmax=1.5 nsteps=8000 energy_every=0 stack=runs/arc.json out=runs/arc
python3 pulse_analysis.py runs/arc                         # R(f) vs the exact lattice reference
python3 tests/disp_gate.py                                 # gate D
```

Units: λ0 = 1 in the solver. `lambda0_um` in the stack file sets the physical scale, so f = λ0/λ. With
λ0 = 0.6 µm, the band 400–800 nm is f = 1.5 … 0.75.

## 1. Materials (`materials/`)

| Material | Data (refractiveindex.info, CC0) | Lorentz poles | RMS rel. ε error, 400–800 nm |
|---|---|---|---|
| Si | M. A. Green, Sol. Energ. Mat. Sol. Cells 92, 1305 (2008) | 3 | 9.3e-3 |
| GaAs | D. E. Aspnes et al., J. Appl. Phys. 60, 754 (1986) | 4 | 8.2e-3 |
| Si3N4 | K. Luke et al., Opt. Lett. 40, 4823 (2015) | 1 | 3.5e-4 |
| SiO2 | I. H. Malitson, JOSA 55, 1205 (1965) | 1 | 2.1e-4 |

The model is ε(ν) = ε∞ + Σ Δε ν0²/(ν0² − ν² − iγν), with time dependence e^{−iωt}. The fit uses variable
projection: ε∞ and Δε are solved by linear least squares, ν0 and γ by Levenberg–Marquardt, with a fixed-seed
multistart. The constraints Δε ≥ 0, γ ≥ 0 and ε∞ ≥ 1 make the model passive and causal, the property Lumerical's
MCM also enforces. The error metric is the one MCM reports: RMS of |ε_model − ε_data|/|ε_data|.

Lorentz poles cannot follow the Si and GaAs absorption edges exactly; the fit error is about 1%. Gate D5 shows the
effect on R by computing the continuous TMM (the Lumerical STACK equivalent) with the model and with the raw data.

## 2. Lorentz ADE

    (1 + a) P^{n+1} = (2 − ω0²Δt²) P^n − (1 − a) P^{n−1} + Δε ω0² Δt² E^n,      a = γΔt/2
    ε∞ (E^{n+1} − E^n) = Δt (∇×H)^{n+1/2} − Σ (P^{n+1} − P^n)

This is the explicit second-order scheme used by Meep. P^{n+1} depends only on E^n, so the update needs no
implicit solve. With D^n = ε∞E^n + ΣP^n, the scheme is exactly the Yee update with the discrete-time permittivity

    ε_FD(ω) = ε∞ + Σ Δε ω0²Δt² / [(1 + a) z − 2 + ω0²Δt² + (1 − a)/z],    z = e^{−iωΔt}

(`disp_tmm.eps_fd`). ε_FD → ε(ω) as O(Δt²). The solver refuses ω0Δt ≥ 2 and ε∞ < 1.

## 3. Subpixel averaging of planar layers

| Component | Node | Averaging interval | Mixing |
|---|---|---|---|
| Ex, Ez (tangential) | primal y_j | dual cell [y_{j−½}, y_{j+½}] | parallel: ε∞,eff = Σ f_m ε∞,m, polarization Σ f_m P_m (all driven by the node E) |
| Ey (normal) | dual y_{j+½} | primal cell [y_j, y_{j+1}] | series: D common, E_m = (D − Σ_k P_mk)/ε∞,m, Ey = Σ f_m E_m |

Tangential E is continuous across the interface, so parallel mixing is exact there. Normal D is continuous, so
series mixing gives the harmonic mean 1/ε = Σ f_m/ε_m at every frequency. This is anisotropic averaging
(Farjadpour et al., Opt. Lett. 31, 2972 (2006); Meep), extended to dispersive materials. Because the series rows
keep one E_m per material, the frequency dependence of each material is handled exactly.

Non-dispersive rows fold into the legacy coefficients: arithmetic ε for tangential rows, harmonic ε for normal rows.
Dispersive rows add `stack_pre_E` and `stack_post_E` around the unchanged `update_E`. Bloch copies keep separate
polarization state. Staircase (`smode=c`) samples the material at each component's own position.

The exact reference `disp_tmm.lattice_rt` solves the same lattice in the frequency domain. It uses the weights from
`stack_used.json`, ε_FD per material, K~t and ω~, with complex ε and a decaying or outgoing wave in the exit medium.
The FDTD must reproduce it up to the CPML reflection and DFT truncation.

## 4. Validation (gate D) and the Lumerical CMT comparison

The test structures cover 400–800 nm (λ0 = 600 nm), with a Bloch cell and kx/k0 = 0.2 (θ = 11.5° at 600 nm):

- **A:** Si3N4 75 nm on Si.
- **B:** GaAs 50 nm on SiO2.
- **C:** (SiO2 103 nm / Si3N4 73.5 nm) × 3 on Si.

| Item | What |
|---|---|
| D0-1 | non-dispersive stack vs legacy medium: same fields |
| D1 | FDTD vs exact lattice, A/B/C, s and p, subpixel and staircase: R, T |
| D2 | convergence vs continuous TMM: subpixel order ≈ 2; staircase reported |
| D3 | thickness sweep at fixed Δ: subpixel follows the continuous curve; staircase steps |
| D4 | how much finer staircase must be to match subpixel (analog of the Lumerical CMT claim) |
| D5 | fit quality and its effect on R (TMM with model vs raw data) |

Lumerical's knowledge-base article "Selecting the best mesh refinement option in the FDTD simulation object" reports,
for its thin-film example, that "CMT at 10 points per wavelength beats staircase at 34" and that staircase would need
"more than 60 points per wavelength" over a 600 nm band. The page could not be retrieved from this environment
(HTTP 403), so its exact structure and numbers are not reproduced here. D3 and D4 are the analogous measurements
for structures A–C. Comparing against that example needs its setup (materials, thicknesses, band, angle, mesh sizes).

### Results (results/disp_gate.json, figures/disp/)

- **D1, FDTD vs the exact lattice (λ0/40 = 15 nm, 300 T0):** max |ΔR| ≤ 6.6e-7 and |ΔT| ≤ 4.2e-7 for A, B and C,
  s and p, subpixel and staircase. This covers the ADE, both mixing rules, Bloch and the pulse together.
- **D2, subpixel convergence to the continuous TMM:** order 2.20 over λ0/160 … λ0/640. That is the upper edge of
  the [1.8, 2.2] window; λ0/80 is still pre-asymptotic, see results/nu_failures.json. Staircase: 1.58 and erratic.
- **D3, thickness sweep (Si3N4 48–102 nm on Si, R at 600 nm):** max error, subpixel vs staircase:
  6.6e-3 vs 0.107 at Δ = 15 nm, and 1.8e-3 vs 0.046 at Δ = 7.5 nm.
- **D4, staircase mesh needed to match subpixel** (RMS over thickness × band):
  - subpixel at 15 nm (12.7 pts/λ in Si3N4 at 400 nm) ≈ staircase at 6.1 nm (31 pts/λ);
  - subpixel at 7.5 nm is not reached by staircase even at λ0/640.

  This is the same order as Lumerical's statement that CMT at 10 points per wavelength beats staircase at 34.
- **D5, effect of the fit:** max |R(model) − R(data)| = 2e-3 (A), 1.1e-2 (B, GaAs), 1.2e-3 (C).
