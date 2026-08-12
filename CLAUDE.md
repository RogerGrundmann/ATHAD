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

## Where the physics lives

| File | What it owns |
|---|---|
| `MixtureAtm.h` | Composition → mass fractions, `R_of`, `cp_of` (Shomate), `M_of`, `M_nonwater`, water's critical point |
| `SaturationH2O.h` | IAPWS saturation + sublimation curves, Watson `latentHeat(T)`, exact `saturationMassFraction`, `dewPoint` (bisection), `dqSatdT` |
| `MultiLayerRadiation.h` | Grey optical depth from column mass with pressure broadening; surface energy balance |
| `ThermoAtm.h` | Densities and the hydrostatic column; `printColumnProfile` / `printLevelSummary` diagnostics |
| `test/saturation_selftest.cpp` | IAPWS reference-point checks — `make test` |

## Four invariants — do not silently break these

1. **There is no topography.** `h ≡ 0`, `i_topography ≡ 0` everywhere, so
   `AtomUtils::is_land()` is false at every point. The Hadean surface is unknown; a
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
   the profile, not nudge it toward a prescribed one — **and it currently does not**:
   `ThermoAtm::densities()` re-imposes the adiabat on `t` every iteration, overwriting
   whatever the dynamics and the radiation computed. That is why the OLR is an input.
   Fixing it is the open task, not a licence to restore a prescribed target.

4. **The column is on its own adiabat, integrated not fitted.** `dT/dz = −g/cp` with local
   cp, hydrostatic on the layer-mean T, isothermal above `t_skin`. Do **not** restore the
   COSMO `T = T₀√(1−coeff·h)` form: it is a sqrt in height, so matching its near-surface
   slope to the adiabat does not make it an adiabat — it reaches zero at 156 km, inside
   the domain.

## Assumptions vs. results — read this before quoting any number

The model reproduces its design targets exactly and its energy balance closes. That does
**not** make its outputs predictions. These are inputs, in rough order of how much they
move the answer:

| Parameter | Value | Status |
|---|---|---|
| `kappa_H2O` / `kappa_CO2` / `kappa_bg` | 0.01 / 0.001 / 1e-6 m²/kg | **Biggest lever on OLR**, factor-of-2 uncertain |
| `geothermal_flux` | 150 W/m² | See below — the model now argues against this value |
| `t_surf_equator` / `t_surf_pole` | 1500 / 1450 K | **Prescribed, not solved** |
| `t_skin` | 254.0 K | From energy balance, but clear-sky albedo — not a fixed point |
| insolation | 0.71 S₀ | Faint young Sun at 4.4 Ga |
| `omega` | 3.17e-4 (5.5 h day) | Estimates range 4–6 h |
| `cosmo_lapse_fraction` | 1.0 (dry adiabat) | Justified: nothing condenses in the deep column |

A grey scheme also cannot represent the window regions that set the real runaway limit.

## What the model currently says

**The model does not currently produce an outgoing longwave flux.** That is the headline,
and it retracts the previous one. Everything below is measured; see README items 9 and 10.

- **OLR ≡ σ·T_lid⁴.** The topmost layer is optically thick (τ ≈ 6), so the model emits at
  its lid temperature — which is prescribed, not solved. Verified by experiment: setting
  `t_skin` to 254 K gives a "measured" OLR of 236.01 W/m², setting it to 240 K gives
  188.13, and σT⁴ is 236.01 and 188.13. The Phase 7 result "OLR = 236.0 W/m², energy
  balance closes" was that identity plus a **stale lid pin**: `bcRadius` held the lid at a
  snapshot taken from `initTemperatureData`, which builds its adiabat before the water and
  CO₂ fields exist, so `densities()` could rebuild the column every iteration and never
  move the lid. `densities()` now refreshes the snapshot, and the lid reads 344 K.
- Current state: shell 230 km, top **0.29 bar — three times above the 0.1 bar radiating
  level**; lid 344 K with ε = 1.0; OLR 802 W/m² against 271 W/m² absorbed + geothermal, an
  imbalance of **−531 W/m²**. This is a truer description than the closed budget it
  replaces.
- **The shell cannot be deepened yet.** 260 km and 300 km both build a sane initial state
  and then NaN across the whole field in the first `MultiLayerRadiation` call. The Thomas
  assembly degenerates as ε → 0, which is the condition at the top of any domain that
  reaches the radiating level. **The radiation scheme is the next task.**
- **Mean planetary albedo 0.4999 = `albedo_cloud`.** The reflectivity saturates the moment
  any condensate exists, so the model reports that parameter rather than computing an
  albedo.
- Insolation fixed: the parameters were Earth's *surface*-absorbed fluxes used as TOA
  insolation, lighting the planet at 45 % of 0.71 S₀. Now a zero-obliquity annual-mean fit
  constrained to S/4; measured mean **241.56 W/m²**.
- `t_skin` is now a fixed point against the model's own albedo, converging 254 → 262.85 K.
- The **geothermal ≥ ~195 W/m²** claim was derived from the OLR identity above and does not
  survive it. The model currently makes no claim of its own.

Bit-identical at 1, 4 and 8 OpenMP threads. Text diagnostics print every 10 iterations for
short runs (`nm ≤ 100`), every 100 for longer ones; `diagnostic_stride` overrides.

## Relationship to the family

Siblings live beside this directory: `ATOM_Precipitation` (modern Earth), `ATJUP`,
`ATSAT`, `ATURAN`, `ATNEPT` (giants), `ASTIM` (impacts).

C++ class, file and function names are kept **identical to `ATOM_Precipitation`** so fixes
cherry-pick in both directions; only the outer shell is renamed (`libathad.a`, `cli/had`,
`config_athad.xml`, `pyathad`). Preserve that.

**Sixteen defects found in the inherited code so far, all latent on Earth and live here.**
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
arithmetically a no-op.

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

- **The radiation scheme is the blocking task.** Its tridiagonal solve degenerates for
  optically thin layers, so the domain cannot reach the radiating level; and the profile
  it is supposed to determine is overwritten by `densities()` every iteration. Until both
  are fixed the OLR is an input. `ThermoAtm::printPlanetaryBalance` prints the lid
  temperature and emissivity next to the OLR so the identity stays visible.
- **The surface temperature is prescribed, not solved.** Every result is conditional on it.
- **Boussinesq.** The solver rests on the Boussinesq buoyancy approximation, but density
  varies by ~2 orders of magnitude across the column. This may force an anelastic or
  compressible formulation. The family's partial answer is the ATJUP hydrostatic split
  (ported in ATURAN `302a51e`) — and it did not cure the giants' problem. **Untested here.**
- **Deep convection is inactive.** Its trigger thresholds (1000/970/900/800 hPa) are
  absolute Earth surface pressures and never fire at 250 bar. They need to become
  fractions of surface pressure.
- `time_start/end/step` remain because the time-slice loop is still structural, though only
  one slice ever runs.
