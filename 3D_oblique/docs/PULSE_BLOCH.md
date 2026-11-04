# Bloch periodic boundaries and broadband pulse

Two additions to `fdtd3d_oblique.c`, both off by default. A run without the new keys follows the legacy code
path unchanged (bit-identical fields; gate 0-1a and gate P0-1).

| Key | Default | Meaning |
|---|---|---|
| `bloch` | 0 | 1: complex Bloch fields; x/z ghost layers carry the phase e^{i kx Lx}, e^{i kz Lz} |
| `kxn`, `kzn` | 2πm/Lx, 2πn/Lz | (bloch=1 only) kx/k0, kz/k0 with k0 = 2π/λ0; any value, the cell need not hold whole periods |
| `src` | c | `c`: legacy CW with ramp; `g`: Gaussian pulse (needs `inc=a`) |
| `fc`, `fw` | 1, 0.15 | pulse centre frequency and spectral standard deviation (units c/λ0) |
| `pt0` | 6/(2π fw) | pulse delay (units T0 = λ0/c) |
| `nfreq`, `fmin`, `fmax` | 0, 0.5, 1.5 | >0: write spectra at `nfreq` equally spaced frequencies in [fmin, fmax] |

Example (half-space ε = 2.25; continuous angle 33.6° at f = 1, 44.4° at 0.8, 26.4° at 1.25):

```bash
./fdtd3d_oblique bloch=1 Lx=0.2 Lz=0.2 kxn=0.37 kzn=0.21 pol=s eps2=2.25 y1=60 tf=120 npml=40 \
    src=g fc=1 fw=0.1 nfreq=51 fmin=0.75 fmax=1.25 nsteps=12000 energy_every=0 out=runs/demo
python3 pulse_analysis.py runs/demo      # R(f), T(f), FDTD vs exact discrete TMM
```

## 1. Bloch periodic boundaries (`bloch=1`)

The y axis (propagation) keeps its CPML at both ends. x and z stay periodic, with a Bloch phase:
F(x + Lx) = e^{i kx Lx} F(x) and F(z + Lz) = e^{i kz Lz} F(z). The complex field F = F_re + i F_im is stored as two real
copies of the legacy state (fields, CPML ψ, TF/SF incident planes). Both copies use the same update kernels and
coefficients. They are coupled only where the ghost layers are filled:

    re[ghost] = cos(kL) re[src] ∓ sin(kL) im[src]
    im[ghost] = ±sin(kL) re[src] + cos(kL) im[src]

The upper sign is for forward ghosts (i = Nx, k = Nz, filled before the H update). The lower sign is for backward
ghosts (i = −1, k = −1, filled before the E update). The TF/SF plane feeds copy 0 with Re[a e^{iφ}] and copy 1 with
Im[a e^{iφ}], where a is the aux-line value.

Consequences:
- **The angle no longer fixes the cell size.** A layered structure needs only a few cells across (the gates use
  4 × 4). A 4 × 4 Bloch cell gives the same spectra as a 40 × 60 real supercell to 2e-15 (P0-3), about 8× faster.
- **Copy 0 is unchanged when the Bloch phase is 1.** If kx Lx/2π and kz Lz/2π are integers (within 1e-12), the
  phase is set to exactly 1. Copy 1 then decouples, and copy 0 equals the legacy run bit for bit (P0-1).
- **Energy and diagnostics.** The conserved energy is the sum over both copies, |F|². The SF maximum and the
  divergence are the worst of the two copies. Snapshots, `dft_<slice>.bin`, `dft_proj.bin` and `fields_final.bin`
  hold copy 0 only; `fields_final_im.bin` holds copy 1.
- **Limits.** Bloch needs uniform x and z (`inc=a`, or `inc=0` for source-free tests). It is refused with `auxref`,
  `mode=e`, `planar`, `init=a`, `inc=m|p|j|n`. The y axis may be nonuniform.

## 2. Gaussian pulse (`src=g`) at fixed (kx, kz)

The aux line is the exact 1D reduction of the main grid for one transverse Floquet mode. It is broadband: the
reduction holds at every ω. In pulse mode the 1D TF/SF source of the aux line is replaced by a soft sheet current at
j_a:

    J(t) = exp(−(t − t0)² / 2τ²) · e^{−i 2π fc t},   τ = 1/(2π fw)

The current points along the tangential part of the s or p unit vector. For a layered medium s and p decouple at
every frequency, so the current excites only the chosen polarization. J is analytic (one-sided spectrum). The
negative-frequency content is ~exp(−(fc/fw)²/2), 1e-22 for the defaults. The main grid is fed through the TF/SF
plane exactly as in CW mode.

**Fixed kx means the angle depends on frequency.** The relation is sin θ(f) = |k_t| c/(2π f) (`pulse_analysis.angle_deg`).
There is a cutoff at f_cut = |k_t|/2π; below it the incident wave is evanescent. To build R(θ) at a fixed angle,
run several kx values and interpolate R(f, kx).

**Cutoff ringing.** At f_cut the 1D system rings with a slowly decaying tail (~t^{−1/2}). If the pulse carries
energy at f_cut, a finite DFT window truncates that tail, and the leakage spreads over the whole band. Measured:
kx/k0 = 0.5, kz/k0 = 1/3, fw = 0.15 gives |ΔR| ≈ 1e-3 near f = 0.7. With fw = 0.08 at the same k, or with
|k_t|/k0 = 0.3, the error drops to ≈ 2e-7. The solver prints a WARNING when |J(f_cut)|/|J(fc)| > 1e-3.
Keep fc − f_cut ≥ 4 fw, and run long enough for the slowest judged frequency to leave the grid (gates: 300 T0).

## 3. Spectra (`nfreq > 0`)

For every step and each y row of the main grid, each component is projected on its own transverse phase:

    a_c(j, t) = s/(Lx Lz) Σ_xz w F e^{−iφ_c}      s = 1 (Bloch complex field) or 2 (real field)

The weights w and the phases φ_c are those of `proj`. The aux line is recorded at the same rows. The spectrum
S(j, ω) = Σ_n a(j, t_n) e^{iωt_n} Δt is accumulated over the whole run, with each component at its own time level.

Output:
- `spec_proj.bin` (main grid) and `spec_aux.bin` (incident field), both complex128 with layout
  `[comp Ex..Hz][freq][j = 0..Ny]`.
- In `meta.json`: `spec`, `pulse`, `bloch`, `bloch_phase_x/z`.

`pulse_analysis.py` computes the discrete-conserved flux through each dual face,
Φ_j = ½ Re[E_z,j H*_x,j+½ − E_x,j H*_z,j+½]. From it:
- R = −Φ(main-grid SF rows) / Φ(aux line), since the SF region holds only the reflected field.
- T = Φ(transmitted rows) / Φ(aux line).

It compares both with `tmm.tmm_discrete`, the exact Yee-lattice R and T at the same ω, K~t, Δt and y grid.
Coordinates come from `grid_used.json`.

## 4. Validation (gate P)

`python3 tests/pulse_gate.py` writes `results/pulse_gate.json` and `figures/pulse/`. The thresholds are in
`tests/thresholds_pulse.json`, separate from the frozen `thresholds_nu.json`, which is not touched.

| Item | What |
|---|---|
| P0-1 | bloch=1 at an integer Bloch phase: copy 0 vs the legacy solver, fields at n = 1…1000 |
| P0-2 | Bloch energy conservation, random init, PEC y walls, non-commensurate k, 20000 steps |
| P0-3 | spectra of a 4 × 4 Bloch cell vs a real 40 × 60 supercell, s and p |
| P1 | half-space, kx/k0 = 0.37, kz/k0 = 0.21, s and p: R, T vs exact discrete; R + T − 1 |
| P2 | normal incidence, s and p |
| P3 | nonuniform-y film F1 (L2_F1_ppw40 y grid, 4 × 4 Bloch cell), s and p |
| P4 | cutoff-leakage control: warning printed; error reported as INFO |
