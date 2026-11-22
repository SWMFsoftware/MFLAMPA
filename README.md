# M-FLAMPA: Multiple Field Line Advection Model for Particle Acceleration

M-FLAMPA is the solar energetic particle (SEP) component (SP) of the Space
Weather Modeling Framework (SWMF), and the particle solver of the SOFIE model
(Solar wind with Field lines and Energetic particles). It solves kinetic
transport equations for energetic protons along a large set of time-dependent
magnetic field lines extracted from the AWSoM(-R) solar corona (SC) and inner
heliosphere (IH) MHD solutions, with the CME-driven shock traced on each line.

The code provides two transport solvers on the same field-line grid, selected
at configuration time:

| Equation | Distribution | Phase space per line | Configuration |
|----------|--------------|----------------------|---------------|
| Parker (diffusive) transport equation | omnidirectional f(s, p, t) | (p³/3, s_L) | `nMu = 1` (default) |
| Focused transport equation (FTE) | gyrotropic, pitch-angle resolved f(s, p, μ, t) | (p³/3, μ, s_L) | `nMu > 1` |

Both run in standalone mode (MHD field-line data read from files) and coupled
mode (field lines advected and fed by SC/IH through the SWMF coupler), and both
have a time-accurate mode and a steady-state (local time stepping) mode.

## Focused transport solver

### Equation solved

For `nMu > 1`, M-FLAMPA advances the gyrotropic distribution f(s_L, p, μ, t)
on each Lagrangian field line (mesh points move with the plasma, so d/dt below
is the Lagrangian derivative), written as a sum of Poisson brackets
(`src/ModAdvance.f90`, `src/ModAdvancePoisson.f90`):

```
df/dt + {f; H_13}_(p³/3, s_L) + {f; H_23}_(μ, s_L) + {f; H_12}_(p³/3, μ)
      = ∂/∂μ (D_μμ ∂f/∂μ)  [ + B ∂/∂s (D_xx/B ∂f/∂s) ]
```

The physical processes carried by each term:

| Process | Term in the code | Control |
|---------|------------------|---------|
| Streaming along the field, μv ∂f/∂s, with relativistic v(p) | `Hamiltonian23_N`, H_23 = (1 − μ²) v / (2B), integrated over the momentum cell | always on |
| Adiabatic focusing (magnetic mirroring) in the diverging field, conserving the magnetic moment | same bracket {f; H_23}_(μ, s_L), through the 1/B dependence | always on |
| Adiabatic cooling and first-order Fermi (shock) acceleration from plasma compression, D ln ρ/Dt, including its μ dependence | time-dependent brackets `dHamiltonian01_FX` {τ, p³/3} and `dHamiltonian02_FY` {τ, μ} with the time-varying control volume Δs/B (time-accurate); `Hamiltonian13_N`, H_13 = −(u/B) p³/3, and the focusing part of `Hamiltonian12_N` (steady state) | always on |
| Betatron acceleration, D ln B/Dt | `Hamiltonian12_N`, ∝ μ(1 − μ²)/2 · p³ · D(1/B)/Dt / (1/B) | `UseBetatron` |
| Inertial force in the non-inertial plasma frame, b · Du/Dt | `Hamiltonian12_N`, ∝ (1 − μ²)/2 · p² γ m_p · D(b·u)/Dt | `UseInertialForce` |
| Pitch-angle scattering, D_μμ = v/λ_μμ · (1 − μ²) \|μ\|^(2/3) (quasi-linear, Kolmogorov turbulence), floored near μ = 0 | `scatter_distribution` in `src/ModDiffusion.f90`, implicit tridiagonal solve in μ for every (p, s_L) | `UseMuScattering` |
| Parallel spatial diffusion (optional alternative closure) | `diffuse_distribution` | `#DIFFUSION` |

The scattering mean free path λ_μμ is tied to the same parallel mean free path
used by the Parker solver (λ_μμ = (14/27) λ_xx, `set_diffusion_coef`), which is
computed from the Alfvén wave turbulence amplitudes (`Wave1_`, `Wave2_`)
delivered by AWSoM along each line, or from the prescribed upstream mean free
path. The Parker and focused solutions therefore use one consistent
turbulence description and can be compared directly on the same event.

### Numerical method

- **Poisson bracket finite volume scheme** (Sokolov et al. 2023, J. Comput.
  Phys., doi:10.1016/j.jcp.2023.111923; `src/ModPoissonBracket.f90`). Each
  term of the FTE is cast as a Poisson bracket {f; H} with a known Hamiltonian
  function, discretized with node-centered Hamiltonians on a (p³/3, μ, s_L)
  control volume. The scheme is conservative (particle number is conserved to
  round-off), second-order accurate in space with TVD limiters, and handles
  multiple brackets plus brackets with respect to time (moving Lagrangian mesh,
  time-dependent cell volume) in one unsplit update (`explicit3`).
- **Shared infrastructure with the Parker solver.** The same module provides
  the single-bracket Parker advance (`advect_via_poisson_parker`) and the
  multi-bracket focused advance (`advect_via_poisson_focused`). Shock
  tracing and steepening (`SP_ModShock`), suprathermal injection at the shock
  and momentum boundary conditions (`SP_ModBc`), sub-cycling so that the shock
  crosses one mesh per sub-step, MPI decomposition over field lines, restart,
  and all output (spectra, energy-channel fluxes, satellite time series) are
  common to both solvers.
- **Grid.** Uniform μ grid on [−1, 1] with `nMu` cells and reflective
  (symmetry) ghost cells at μ = ±1; logarithmic momentum grid with `nMomentum`
  cells; up to `nVertexMax` Lagrangian vertices per line.
- **Time stepping.** Explicit CFL-limited advection over each MHD coupling
  interval, with background ρ and B linearly interpolated between consecutive
  MHD states; implicit scattering/diffusion applied by operator splitting
  after every advection sub-step. Steady-state mode
  (`iterate_poisson_focused`) uses a local time step per (p, μ, s_L) cell.
- **Output.** Distribution plots (`distr*` in `#SAVEPLOT`) retain the μ index
  when `nMu > 1`; fluxes in energy channels are obtained by μ integration, so
  GOES-type and user-defined channel products are available from either solver.

### How to enable

```bash
./Config.pl -g=20000,100,5        # nVertexMax, nMomentum, nMu  (nMu > 1 selects the FTE)
make MFLAMPA                      # standalone; or `make LIB` for the SWMF/SOFIE build
```

```
#ADVECTION
T			UsePoissonBracket

#FOCUSEDTRANSPORT
T			UseBetatron
T			UseInertialForce
T			UseMuScattering
```

The focused solver requires `UsePoissonBracket = T`. Setting `nMu = 1`
recovers the Parker solver with no other change to PARAM.in.

## Parker transport solver

### Equation solved

For `nMu = 1`, M-FLAMPA advances the omnidirectional (pitch-angle averaged)
distribution f(s_L, p, t), which applies when pitch-angle scattering keeps
the particles nearly isotropic in the plasma frame. Each field line is a thin
magnetic flux tube with cross section ∝ 1/B, and its mesh points move with the
plasma, so D/Dt below is the Lagrangian derivative:

```
Df/Dt + (1/3) (D ln ρ/Dt) ∂f/∂ln p = B ∂/∂s (D_xx/B ∂f/∂s)
```

along each Lagrangian line, with either the legacy upwind scheme in ln p
(`src/ModAdvanceAdvection.f90`) or the conservative Poisson bracket scheme
(`src/ModAdvancePoisson.f90`, `#ADVECTION`).

The physical processes carried by each term:

| Process | Term | Control |
|---------|------|---------|
| Advection with the plasma along the line | Df/Dt: the mesh moves with the plasma | always on |
| Adiabatic cooling where the plasma expands, D ln ρ/Dt < 0 | (1/3) (D ln ρ/Dt) ∂f/∂ln p | always on |
| First-order Fermi (diffusive shock) acceleration where the plasma is compressed, D ln ρ/Dt > 0, mainly across the CME-driven shock; particles gain momentum each time diffusion carries them back across the shock, giving f ∝ p^(−3r/(r−1)) for compression ratio r | the same term, together with parallel diffusion | always on |
| Parallel spatial diffusion, D_xx = v λ_xx / 3, along the diverging flux tube | B ∂/∂s (D_xx/B ∂f/∂s), `diffuse_distribution` in `src/ModDiffusion.f90` | `#DIFFUSION` |

The parallel mean free path λ_xx is computed from the Alfvén wave turbulence
amplitudes (`Wave1_`, `Wave2_`) delivered by AWSoM along each line, or from
the prescribed upstream mean free path. Particles are injected from the
suprathermal tail at the shock (`SP_ModBc`). Optional perpendicular diffusion
across lines (`src/ModDiffusionPerp.f90`) and self-generated turbulence
(`src/ModTurbulence.f90`) are available.

In the steady-state mode (`#TIMEACCURATE F`) the MHD background is fixed and
the plasma velocity u along the line enters explicitly. Advection and the
adiabatic momentum change then form one Poisson bracket on the (p³/3, s_L)
plane, with the reduced Hamiltonian H̃ = −(u/B) p³/3:

```
∂f/∂t + B {f; H̃}_(p³/3, s_L) = B ∂/∂s (D_xx/B ∂f/∂s)
```

since B {f; H̃} = u ∂f/∂s − (1/3) (∇·u) p ∂f/∂p and ∇·u = B ∂(u/B)/∂s in the
flux tube. The control volume is Δ(p³/3) Δs/B. The difference δH̃ of H̃
across a cell face gives the flux through it: along s_L it is advection with
the flow, along p³/3 it is the adiabatic momentum change, upward where the flow
is compressed and downward where it expands. The steady state is reached by
iterating in pseudo-time.

### Numerical methods

- **Time accurate** (`#TIMEACCURATE T`, the default): explicit, CFL-limited
  (`#CFL` ≤ 1) advection in momentum with either the legacy upwind scheme in
  ln p (`src/ModAdvanceAdvection.f90`, `#ADVECTION F`) or the conservative
  Poisson bracket scheme (`src/ModAdvancePoisson.f90`, `#ADVECTION T`);
  parallel diffusion is applied by operator splitting after every advection
  sub-step.
- **Steady state, explicit** (`#TIMEACCURATE F`): `iterate_poisson_parker`
  always uses the Poisson bracket scheme, whatever `#ADVECTION` says. It takes
  a local pseudo-time step in each (p, s_L) cell at `Cfl` ≤ 1 and then applies
  the diffusion by operator splitting.
- **Steady state, implicit** (`#IMPLICIT T`): one backward-Euler step that
  couples advection, the adiabatic momentum change and parallel diffusion,
  with local pseudo-time steps far above the explicit limit (below).

### Implicit steady-state solver

Each iteration is one backward Euler step on the (p³/3, s_L) plane of each
line, written here for the increment δf = f^(n+1) − f^n:

```
(1 + Σ_k b_ck) δf_c − Σ_k b_ck δf_n(k) = τ_c R_c(f^n),    b_ck = τ_c (|δH̃_k⁻| + D_k) ≥ 0
```

- R_c(f^n) is the steady-state residual of cell c at the old state: upwind
  fluxes |δH̃_k⁻| (f_n(k) − f_c) through the inflow faces (δH̃_k < 0),
  diffusion fluxes D_k (f_n(k) − f_c) through the faces along s_L, and the
  second-order correction with the Koren limiter.
- τ_c = Δt / V_c, with one pseudo-time step for the line,
  Δt = Cfl · min_c (V_c / Σ_k |δH̃_k⁻|), the cell volume V_c = Δ(p³/3) Δs/B
  and `Cfl` from `#CFL`.
- The matrix is solved by Gauss–Seidel iterations over the momentum rows,
  each row a tridiagonal (Thomas) solve along s_L, with at most 10 iterations
  and a relative tolerance of 10⁻⁵ (`implicit2` of `ModPoissonBracket`).
- δf = 0 exactly when R(f^n) = 0, so the converged solution is the steady
  state of the second-order scheme for any `Cfl`; `Cfl` sets only how fast it
  is reached.

The kernel is `implicit2` in `src/ModPoissonBracket.f90`, the driver
`iterate_implicit` in `src/ModAdvancePoisson.f90`, and the diffusion face
conductances `get_diffusion_face` in `src/ModDiffusion.f90`.

### How to enable the implicit solver

```bash
./Config.pl -g=20000,100,1        # nVertexMax, nMomentum, nMu  (nMu = 1: Parker transport)
make MFLAMPA
```

```
#TIMEACCURATE
F			IsTimeAccurate

#IMPLICIT
T			UseImplicit

#CFL
9800.0			Cfl
```

`#IMPLICIT` (default `F`) acts only in the steady-state Parker solver;
time-accurate runs and the focused transport solver (`nMu > 1`) are
unchanged. `Cfl` is the Courant number of the local
pseudo-time steps and may far exceed 1 here; with the default 0.9 the
iterations are valid but converge slowly, so a large value such as 9800 should
be set. With `#VERBOSE` above 1, every iteration prints the number of sweeps,
the Gauss–Seidel error and the range of Δt_c.

### Steady-state GCR runs

For the steady-state galactic cosmic ray (GCR) spectrum along the lines, no
particles are injected: the spectrum enters through the outer end of each line
and is modulated inward by advection, adiabatic cooling and diffusion.

```
#TIMEACCURATE
F			IsTimeAccurate

#IMPLICIT
T			UseImplicit

#CFL
9800.0			Cfl

#MOMENTUMBC
grad			TypeMomentumMinBc
grad			TypeMomentumMaxBc

#LOWERENDBC
float			TypeLowerEndBc

#UPPERENDBC
lism			TypeUpperEndBc
usoskin			TypeLisBc
T			UseModulationPot
534.7			ModulationPot

#STOP
100			nIterMax
-1.0			TimeMax
```

- `#MOMENTUMBC grad grad`: at both ends of the momentum grid the ghost cell
  continues the local power law, f₀ = f₁²/f₂ and f_(nP+1) = f_nP²/f_(nP−1).
  `float` instead sets a zero gradient, f₀ = f₁; `inject` (the default at the
  lower end) adds the suprathermal seed population of SEP runs.
- `#UPPERENDBC lism`: the outer end of each line holds the local interstellar
  spectrum (`usoskin`, Usoskin et al. 2005, or `corti`), modulated by the
  potential `ModulationPot` in MV when `UseModulationPot` is T, e.g. the Oulu
  neutron monitor value for the date of the run.
- `#LOWERENDBC float`: the inner end of each line copies the spectrum of the
  first vertex; `escape` sets the background spectrum there.
- `#STOP` `nIterMax` is the number of implicit iterations; with `#VERBOSE`
  above 1 the Gauss–Seidel error of every iteration shows the convergence.

The nightly steady-state decks (`Param/PARAM.in.test.steady*`) use the same
`#UPPERENDBC` with the explicit iteration and the default momentum boundary
conditions.

## Build, run, test

M-FLAMPA builds inside an installed SWMF tree (`SWMF/SP/MFLAMPA`).

```bash
./Config.pl -s                                   # show grid configuration
make MFLAMPA                                     # bin/MFLAMPA.exe (standalone)
make rundir RUNDIR=run_test STANDALONE=YES SPDIR=`pwd`
make test                                        # nightly test: upwind, Poisson bracket, steady-state/restart
make test_poisson_bracket                        # unit and analytical tests of the Poisson bracket schemes
make test_poisson_bracket_exe                    # build test_poisson.exe only
./test_poisson.exe                               # nightly Poisson bracket tests (also explicitorig)
./test_poisson.exe Comprho4Shock                 # steady shock acceleration test, explicit and implicit
```

Parameters are documented in `PARAM.XML`; examples are in `Param/`; the user
manual is `Doc/Tex/USERMANUAL.tex`.

## References

- Liu, W., Sokolov, I. V., Zhao, L., et al. (2025), Physics-based simulation
  of the 2013 April 11 solar energetic particle event, Astrophys. J., 985, 82,
  doi:10.3847/1538-4357/adc4e3.
- Liu, W., Liu, X., Lario, D., et al. (2026), Counterintuitive magnetic
  connectivity and energetic particle flux differences among nearby
  spacecraft during the 2023 February 24 solar energetic particle event,
  Astrophys. J., 1008, 6, doi:10.3847/1538-4357/ae9097.
- Sokolov, I. V., et al. (2023), High resolution finite volume method for
  kinetic equations with Poisson brackets, J. Comput. Phys., 476, 111923,
  doi:10.1016/j.jcp.2023.111923.
- Borovikov, D., Sokolov, I. V., Roussev, I. I., Taktakishvili, A., and
  Gombosi, T. I. (2018), Toward a quantitative model for simulation and
  forecast of solar energetic particle production during gradual events. I.
  Magnetohydrodynamic background coupled to the SEP model, Astrophys. J., 864, 88.
- Sokolov, I. V., Roussev, I. I., Gombosi, T. I., et al. (2004), A new field
  line advection model for solar particle acceleration, Astrophys. J. Lett.,
  616, L171.
- Zhao, L., Sokolov, I., Gombosi, T., et al. (2024), Solar Wind with Field
  Lines and Energetic Particles (SOFIE) model: Application to historical solar
  energetic particle events, Space Weather, 22, e2023SW003729.

## License

Apache License 2.0, see `LICENSE.txt`. Copyright Regents of the University of
Michigan.

