# ATHAD — Atmosphere of the Earth in the Hadean Eon

An atmospheric general-circulation model of the Earth as it was in the Hadean
(~4.4 Ga): a ~250 bar, water-vapour-dominated atmosphere over a molten or quenching
surface with no known topography.

**Forked from `ATOM_Precipitation` @ `1e3f319` (2026-07-28), atmosphere half only.**
Started 2026-08-11. There is no hydrosphere, no paleogeography, and no time-slice
series — ATHAD is one epoch.

## Build, run, test

```bash
make had                                        # -> cli/had
make test                                       # IAPWS self-test, run this first
cd python && OMP_NUM_THREADS=8 ../cli/had config_athad.xml
```

Output lands in `python/output_Hadean/` — the ATOM line's convention
(`output_<name>/`, underscore, relative to the run directory; the giant-planet
siblings hyphenate instead).

`make` regenerates the parameter bindings from `param.py` first. The generated files
(`atmosphere/*.inc`, `python/atmosphere_pxd.pxi`, `python/pyathad.pyx`,
`cli/config_athad.xml`, `python/config_athad.xml`) are **tracked on purpose**, so a
`param.py` change shows its full effect in the diff. Regenerate and commit them together.

Removing a parameter from `param.py` deletes a C++ member, so it must be removed
together with its uses or the build breaks.

## Atmospheric composition

Mole fractions are the input; the model works in mass fractions. The residual 7 % is
split evenly across the five trace gases.

| Species | Mole frac. xᵢ | Mᵢ [g/mol] | Mass frac. qᵢ | Rᵢ [J/(kg·K)] |
|---|---|---|---|---|
| H₂O | 0.800 | 18.015 | 0.6724 | 461.5 |
| CO₂ | 0.100 | 44.010 | 0.2053 | 188.9 |
| N₂  | 0.030 | 28.014 | 0.0392 | 296.8 |
| CH₄ | 0.014 | 16.043 | 0.0105 | 518.3 |
| NH₃ | 0.014 | 17.031 | 0.0111 | 488.2 |
| H₂  | 0.014 |  2.016 | 0.0013 | 4124.2 |
| CO  | 0.014 | 28.010 | 0.0183 | 296.8 |
| SO₂ | 0.014 | 64.066 | 0.0418 | 129.8 |

- **M_mean = 21.434 g/mol**, **R_mix = 387.9 J/(kg·K)** (dry air is 286.9)
- **Background** (everything except H₂O and CO₂): M_bg = 26.207 g/mol,
  **R_bg = 317.3 J/(kg·K)**. This is what the `R_Air` parameter now means — it is not air.
- **p_surf = 250 bar**, **T_surf = 1500 K** (prescribed), **ρ_surf = 42.97 kg/m³**
- cp ≈ 2040 J/(kg·K), strongly T-dependent across 300–1500 K — roughly 2× Earth's.

Only **H₂O (`c`) and CO₂ (`co2`) are prognostic**, both as **mass fractions** (not ppm —
at 20 % by mass ppm is meaningless). The other six are a fixed well-mixed background
entering R_mix, cp_mix and the opacity.

CO₂ is prognostic *in fact* only since README item 12: `co2Atmosphere()` used to re-impose
a uniform field inside the time loop and discard the transported one. It is now the initial
condition only, and `ThermoAtm::co2Column()` monitors the global mass-weighted mean, which
is conserved (no CO₂ source or sink exists). The field still comes out uniform — with no
gradients a passive tracer has nothing to transport — but that is now computed rather than
asserted.

## Where the physics lives

| File | What it owns |
|---|---|
| `MixtureAtm.h` | Composition → mass fractions, `R_of`, `cp_of` (Shomate), `M_of`, `M_nonwater`, water's critical point |
| `SaturationH2O.h` | IAPWS saturation + sublimation curves, Watson `latentHeat(T)`, exact `saturationMassFraction`, `dewPoint` (bisection), `dqSatdT` |
| `MultiLayerRadiation.h` | Grey optical depth from column mass with pressure broadening; surface energy balance |
| `ThermoAtm.h` | Densities and the hydrostatic column, anchored to `p_0`; `printColumnProfile` / `printLevelSummary` diagnostics; the `ATM_PROGNOSTIC_T` knob |
| `ConvectiveAdjustment.h` | Dry Manabe–Strickler, generalised to the stretched grid and to a local `cp_of()` — the family's shared version assumes neither |
| `test/saturation_selftest.cpp` | IAPWS reference-point checks — `make test` |

## Four invariants — do not silently break these

1. **There is no topography, and the planet is hemispherically symmetric.** `h ≡ 0`,
   `i_topography ≡ 0` everywhere, so `AtomUtils::is_land()` is false at every point.
   Nothing in the model can sustain a north-south asymmetry either: the surface
   temperature is a symmetric parabola, the insolation is explicitly mirrored, and there
   is no obliquity and no seasonal cycle. So any asymmetry in the output is a defect —
   see README item 13, where the Hadley cells differed by 32 % because of one initial
   velocity coefficient inherited from Earth's land-sea-driven ITCZ offset. The Hadean surface is unknown; a
   featureless global surface is the deliberate choice, not a missing data file. Do not
   reintroduce a bathymetry read, and do not "fix" the dead land branches.
   `LandOceanFraction()` throws if a land point ever appears.

2. **Water is supercritical below ~177 km.** Critical point 647.096 K / 220.64 bar.
   Every condensation path must be a genuine no-op there, not a clamp. Use
   `SaturationH2O.h`; **never reintroduce Magnus** (calibrated to ~320 K, returns 1.2e7 hPa
   at 1500 K, which flips the sign of any `p − E` denominator) and never the dilute
   `q_sat = ep·E/(p−E)` (water is the bulk gas, so there is no small parameter).

3. **Radiation runs in mode 2** (direct σT⁴). Modes 0/1/3/4/5 all lean on the Scotese
   snapshot or the 280 ppm CO₂ reference; neither exists at 4.4 Ga. Radiation must *set*
   the profile, not nudge it toward a prescribed one — **and by default it still does not**:
   `ThermoAtm::densities()` re-imposes the adiabat on `t` every iteration, overwriting
   what the dynamics and the radiation computed. The OLR is a real integral over that
   prescribed profile. `ATM_PROGNOSTIC_T=1` switches the overwrite off (README item 19) and
   `ConvectiveAdjustment.h` supplies the mixing the adiabat was silently providing (item 20),
   so the flip is now a measurement rather than a rewrite. **It has not been made, and as of
   item 23 it makes the imbalance worse**: −226 W/m² prognostic against −66 prescribed, all
   of it in the top five levels, where the prescription pins an isothermal `t_skin` = 263 K
   that is itself an assumption. Fixing the prescription is the open task, not a licence to
   restore a prescribed target — but note that the better-looking of the two numbers is the
   more assumed one.

4. **The column is on its own adiabat, integrated not fitted.** `dT/dz = −g/cp` with local
   cp, hydrostatic on the layer-mean T, isothermal above `t_skin`. Do **not** restore the
   COSMO `T = T₀√(1−coeff·h)` form: it is a sqrt in height, so matching its near-surface
   slope to the adiabat does not make it an adiabat — it reaches zero at 156 km, inside
   the domain.

## Assumptions vs. results — read this before quoting any number

The model reproduces its design targets exactly, and run long enough its energy balance
closes to -1.26 W/m2 — **which is the least trustworthy number in this file, not the most**,
because closing is what the `t_skin` fixed point guarantees (item 25). None of this makes
the outputs predictions. These are inputs, in rough order of how much they move the answer:

| Parameter | Value | Status |
|---|---|---|
| `kappa_H2O` / `kappa_CO2` / `kappa_bg` | 0.01 / 0.001 / 1e-6 m²/kg | Factor-of-2 uncertain. **Whether they are a lever on the *converged* OLR at all is untested** — README item 25 |
| `geothermal_flux` | 150 W/m² | Open. The ≥195 W/m² argument is retracted, and item 25 makes it worse: it enters the `t_skin` fixed point, so it helps set the very flux it was being compared against |
| `t_surf_equator` / `t_surf_pole` | 1500 / 1450 K | **Prescribed, not solved** |
| `t_skin` | 254.0 K start, relaxes to 262.96 | **Now the prime suspect** (item 25): it is a fixed point of σT⁴ = absorbed, the prescribed profile's top is isothermal at it, and the converged OLR falls onto it |
| insolation | 0.71 S₀ | Faint young Sun at 4.4 Ga |
| `omega` | 3.17e-4 (5.5 h day) | Estimates range 4–6 h |
| `cosmo_lapse_fraction` | 1.0 (dry adiabat) | Justified: nothing condenses in the deep column |

A grey scheme also cannot represent the window regions that set the real runaway limit.

## What the model currently says

**The model now computes an outgoing longwave flux.** Everything below is measured; see
README items 9-11.

- `MultiLayerRadiation` is two-stream flux sweeps, not the inherited tridiagonal solve:
  `up[i] = up[i-1](1-eps_i) + eps_i sigma T_i^4` upward, the mirror downward, and radiative
  equilibrium `sigma T_i^4 = (up[i-1] + dn[i+1])/2` in which **eps cancels**. Nothing
  divides by eps, so an optically thin top is exact rather than fatal. This removed the
  ceiling on the shell.
- **Shell 300 km**, 61 levels; top 3.8e-4 bar with lid eps = 0.0000, isothermal skin
  resolved from 256 km. **OLR decoupled from sigma*T_lid^4 — a real column integral.** At
  the old 230 km the two were equal and the OLR was an input.
- **THE OPACITY CLAIM IS WITHDRAWN, and what replaces it is a question about the boundary.**
  Every "the atmosphere radiates away more than it takes in" figure this file has carried —
  581 W/m2 through items 11-17, 679.8 and -409 at 400 iterations in item 18, -66.1 at 20
  iterations after item 22 — was measured before the transient had decayed. Run to 200
  iterations the imbalance falls monotonically to **-1.26 W/m2** (README item 25). **But it
  closes by the OLR dropping onto a sigma*T_lid^4 that never moves**: t_skin goes 262.91 ->
  262.96 K across the whole run, and t_skin is itself the fixed point of sigma*T_skin^4 =
  absorbed. Once the effective radiating level migrates into the isothermal skin, OLR =
  absorbed is arithmetic, not a result. **Do not claim anything about `kappa_H2O` = 0.01
  m2/kg, in either direction, until a kappa scan has been run to 200 iterations** and shown
  whether the converged OLR moves with it. If it does not, item 10's "the OLR is an input"
  survived item 11's rewrite and merely hid until iteration 20 had passed.
- **Not grid-converged**: 519 W/m2 at 260 km against 581 at 300 km with `im` fixed at 61
  (both at 100 iterations, both pre-item-22 — the check has to be redone, not just extended).
- Mean planetary albedo 0.4981 ≈ `albedo_cloud`; the reflectivity saturates the moment any
  condensate exists, so the model reports that parameter.
- Insolation is now TOA (mean 241.56 W/m2); `t_skin` is a fixed point against the model's
  own albedo, converging to 262.85 K.
- `radiation.x` is the **upward long-wave flux**, not sigma*T^4; `bcRadius` no longer pins
  the radiation lid.
- The **geothermal >= 195 W/m2** claim of Phase 7 is retracted: it rested on an OLR that
  was `sigma*t_skin^4` plus a stale lid pin (README item 10).

**Reproducibility, stated precisely (README item 18).** Run-to-run **bit-identical at a
fixed thread count**; **thread-count dependent at ~1e-8**. The old claim of "bit-identical at
1, 4 and 8 threads" had silently stopped being true and is not fully restored:

- **Fixed — the data races.** The Poisson loop wrote `p_dyn` in place under
  `collapse(2) schedule(dynamic,4)` over the two indices its stencil reads across, and
  `UtilsAtm::findResiduumAtm` wrote the shared `m.residuum_old` from inside every thread.
  The same binary at the same thread count gave different answers run to run; it no longer
  does. Red-black colouring and a single post-reduction write.
- **Not fixed — floating-point reduction order.** OpenMP combines partial sums in a
  thread-count-dependent order and `+` is not associative, so 1, 4 and 8 threads still
  differ in the last digit (`residuum_atm` 0.74743479 / 0.74743479 / 0.74743481) while the
  printed wind extrema agree exactly. The feedback path is a global mean that re-enters the
  physics; `t_skin` is the prime suspect, being a reduction the whole column is then rebuilt
  from. Curing it needs ordered reductions.

The pressure-solver race is the third time the family has found that defect and the first
time here, despite it being listed below under *traps already solved elsewhere*. **A
cross-reference is not a check.**

Text diagnostics print every 10 iterations for
short runs (`nm ≤ 100`), every 100 for longer ones; `diagnostic_stride` overrides.

**`atom_log.txt` is reserved by the model — never redirect a run into it.**
`lib/Utils.cpp:45` opens `atom_log.txt` in the run directory with `std::ofstream::out`,
which truncates it. A shell redirect to the same name loses the whole console output when
the model's handle closes; a 100-iteration run's printouts were destroyed this way. Use any
other filename.

## Relationship to the family

Siblings live beside this directory: `ATOM_Precipitation` (modern Earth), `ATJUP`,
`ATSAT`, `ATURAN`, `ATNEPT` (giants), `ASTIM` (impacts).

C++ class, file and function names are kept **identical to `ATOM_Precipitation`** so fixes
cherry-pick in both directions; only the outer shell is renamed (`libathad.a`, `cli/had`,
`config_athad.xml`, `pyathad`). Preserve that.

**Nineteen defects found in the inherited code so far, all latent on Earth and live here.**
The pattern is consistent and worth expecting: *Earth's numbers as bare literals inside
physics kernels, each with a comment justifying it by Earth's conditions.* Examples —
`dr = 0.025` silently tied to `im = 41`; a 333.15 K cap written back into the prognostic
temperature; `287.0` J/(kg·K) as the density gas constant; convective triggers as absolute
hPa; `p_stat` cubically extrapolated at the lid. When something behaves oddly, look for a
constant that was true at 1 bar and 288 K.

The three most recent are worth stating because they show the pattern's worst form — an
Earth-only regime written as a *fallback branch*, so it never runs at home and is never
tested: `SaturationAdjustment::clampAndFade` returned `q_sat = ep*1e-5` for superheated
vapour when the correct answer is 1, and so condensed the entire water column in the one
place where nothing can condense; the same file's Newton loop kept the dilute form the
entry point had already been fixed away from; and `AtmMixture::M_nonwater` took only the
CO2 fraction, so the renormalisation "to exclude H2O" its comment promised was
arithmetically a no-op. A fourth of the same shape: `init_tropopause_layers` converted a
height to a level index as `round(h / L_atm)`, which is only an index on a uniform grid —
this one is exponentially stretched, so the pole's convective top was placed at level 12
(13 km) instead of 52 (196 km), and `VelocityInitializer` built the entire initial wind
structure inside the bottom 4 % of the atmosphere.

**The two newest are about order, not about literals, and they are the ones this file
under-rated.** `initTemperatureData` reads `c` and `co2` at every level, but
`initWaterWapour` and `co2Atmosphere` ran *after* it, so the whole initial column was
integrated with R_bg = 317.3 instead of R_mix = 387.9 and a background-only cp — and
`initCloudIce` and the first `SaturationAdjustment` then laid the cloud deck on that column,
50 km too low, in a band where this atmosphere cannot condense. Halving the OLR when
corrected (README item 22). And `p_stat.x[0]` was re-anchored every iteration to
`r_air·R_mix·T_surf`, which holds the surface *density* fixed and lets the mass of a "250 bar
atmosphere" follow the surface temperature (item 19). Both were dismissed in writing here as
harmless because `densities()` rebuilds the column afterwards. **The lesson generalises: an
initialisation defect is not excused by a later overwrite until you have listed everything
that runs in between.**

**Fixes worth porting back upstream** (not yet applied to ATOM_Precipitation as of
2026-08-11): the `t.x[-1]` out-of-bounds in `MoistConvection::findCloudBaseLFS`; the
`m_node_weights` OpenMP race in `GetMean_2D/3D`; the UB in `get_temperatures_from_curve`;
and `-MMD -MP` header dependencies in the Makefile.

Traps already solved elsewhere in the family — check before re-deriving:
Coriolis/centrifugal signs (ATURAN `8b284cb`, `4201957`; ATNEPT `024c37f`, `e412b1b` —
ATHAD's dynamics already agree, its *diagnostics* did not); mass- not mole-weighted mixture
properties (ATNEPT `c116d71`); in-place Gauss–Seidel as a threading defect (ATURAN
`ffd0e0e`); report failures and limits in the README (ATURAN `74b4ded`, ATNEPT `34286b8`).

## Open risks

- **The profile is still prescribed by default.** `densities()` overwrites `t` with the
  adiabat every iteration, so the OLR is a real integral over a profile the radiation did not
  choose. `ThermoAtm::printPlanetaryBalance` prints the lid temperature and emissivity next
  to the OLR, and flags the case where the two coincide. `ATM_PROGNOSTIC_T=1` plus
  `ConvectiveAdjustment` is the intended replacement and is not yet good enough to be the
  default — see README items 19-20 and 22.
- **The surface temperature is prescribed, not solved.** Every result is conditional on it.
- **Boussinesq.** The solver rests on the Boussinesq buoyancy approximation, but density
  varies by ~2 orders of magnitude across the column. This may force an anelastic or
  compressible formulation. The family's partial answer is the ATJUP hydrostatic split
  (ported in ATURAN `302a51e`) — and it did not cure the giants' problem. **Now testable
  here**: an anelastic projection (`∇·(ρ̄u) = 0`, base state, matching Poisson stencil,
  zero mass flux at the walls) is implemented behind `ATM_ANELASTIC` and is **on by
  default** since README item 18. It cuts the anelastic residual 22 % and halves the
  spurious radial wind in the initial projection; it does *not* change the tracer mass
  budget, because that was never a transport error (item 17); and over 400 iterations it
  tracks the Boussinesq baseline to 0.03 K in mean T, 0.002 % in KE and 0.6 % in Ψ_max,
  with the `p_dyn_cap` clamp binding in 0 of 3 791 399 cells. `ATM_ANELASTIC=0` forces the
  Boussinesq path for A/B.
- **The column air mass is conserved now** (README item 19). `p_stat.x[0] = p_0`, a constant,
  because the surface pressure of an atmosphere is the weight of the air above it and the
  surface density is what follows. Drift −0.2128 % → −0.0011 % over 20 iterations, and the
  water "creation" `waterBudget()` was reporting went with it. Do not re-anchor the column to
  `r_air·R_mix·T_surf`.
- **Deep convection is inactive.** Its trigger thresholds (1000/970/900/800 hPa) are
  absolute Earth surface pressures and never fire at 250 bar. They need to become
  fractions of surface pressure.
- `time_start/end/step` remain because the time-slice loop is still structural, though only
  one slice ever runs.
