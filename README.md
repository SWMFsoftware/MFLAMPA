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

For `nMu = 1` the code solves

```
∂f/∂t + (1/3) (D ln ρ/Dt) ∂f/∂ ln p = B ∂/∂s (D_xx/B ∂f/∂s)
```

along each Lagrangian line, with either the legacy upwind scheme in ln p
(`src/ModAdvanceAdvection.f90`) or the conservative Poisson bracket scheme
(`#ADVECTION`). Optional perpendicular diffusion across lines
(`src/ModDiffusionPerp.f90`) and self-generated turbulence
(`src/ModTurbulence.f90`) are available. 

## Build, run, test

M-FLAMPA builds inside an installed SWMF tree (`SWMF/SP/MFLAMPA`).

```bash
./Config.pl -s                                   # show grid configuration
make MFLAMPA                                     # bin/MFLAMPA.exe (standalone)
make rundir RUNDIR=run_test STANDALONE=YES SPDIR=`pwd`
make test                                        # upwind, Poisson bracket, steady-state/restart
make test_poisson_bracket                        # unit tests of ModPoissonBracket
```

Parameters are documented in `PARAM.XML`; examples are in `Param/`; the user
manual is `Doc/Tex/USERMANUAL.tex`.

## References

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

