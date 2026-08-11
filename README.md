# ATHAD

**Atmosphere of the Earth in the Hadean Eon.**

An atmospheric general-circulation model of the Earth at ~4.4 Ga: a finite-difference
Navier–Stokes solver on a spherical shell with RK4 time integration and vertical
coordinate stretching, applied to a ~250 bar, water-vapour-dominated atmosphere over a
molten or quenching surface with no topography.

ATHAD is the atmosphere half of [ATOM_Precipitation](https://github.com/RogerGrundmann/ATOM_Precipitation)
(forked at `1e3f319`) re-based onto Hadean conditions. It keeps that model's turbulence
closures (k-ε, k-ω, k-ω SST), saturation adjustment, Zero/One/Two/Three-Category ice
schemes, moist convection and multi-layer radiation, and replaces everything that was
calibrated to a 1 bar, 288 K, N₂/O₂ Earth.

## Repository layout

```
atmosphere/   C++ source — the atmosphere model
cli/          Command-line interface source and config file
lib/          Shared utilities (arrays, geometry, …)
python/       Python bindings (Cython)
test/         Self-tests
tinyxml2/     Embedded TinyXML-2 XML parser
param.py      Parameter definitions (generates .inc, .pxi and .xml at compile time)
Makefile      Top-level build file
CLAUDE.md     Project constants, invariants and open risks
```

## Prerequisites

| Dependency | Notes |
|---|---|
| g++ ≥ 7 or clang++ ≥ 6 | C++17 required |
| OpenMP | Typically ships with the compiler |
| Cython, Python ≥ 3.6 | Only for the Python interface |

## Build

```bash
make had          # the CLI executable, cli/had
make all          # CLI + Python extension
```

## Run

```bash
cd cli && OMP_NUM_THREADS=1 ./had config_athad.xml
```

Output lands in the directory named by `<output_path>` (default `output-Hadean/`)
as ParaView `.vtk` slices and `.vts` panoramas.

## Configuration

`cli/config_athad.xml` is generated from `param.py` and carries every parameter with
its default and a description. You need only include the entries you want to change;
everything else falls back to the compiled default.

## Hadean conditions

Surface pressure 250 bar, surface temperature ≈ 1500 K, mean molar mass 21.43 g/mol,
R_mix = 387.9 J/(kg·K), surface density 43 kg/m³. Composition by mole fraction:
80 % H₂O, 10 % CO₂, 3 % N₂, and 1.4 % each of CH₄, NH₃, H₂, CO and SO₂. H₂O and CO₂
are prognostic; the rest are a fixed well-mixed background.

Full constants, the three model invariants, and the assumptions that are still open
are in [CLAUDE.md](CLAUDE.md).

## Status and limitations

Under construction. This section records what has actually been measured — including
the measurements that did not work out — rather than what is intended.

1. **Bootstrap (done).** The tree builds and links as `libathad.a` + `cli/had`.
   Hydrosphere, paleogeography and the NASA/Scotese data pipeline are removed.

2. **Flat surface (done).** `init_topography()` prescribes `h ≡ 0`,
   `i_topography ≡ 0`; the run reports 65341 surface points, all water. The surface
   temperature is prescribed from `t_surf_equator` / `t_surf_pole` (1500 K / 1450 K,
   area mean 1483 K) instead of read from a paleo curve. A 2-iteration run completes
   and writes 20 ParaView files, **bit-identically at 1, 4 and 8 OpenMP threads**.

   Three defects had to be fixed to get there, all inherited and all latent on Earth:

   - `MoistConvection::findCloudBaseLFS` indexed `t.x[-1][j][k]`. Its scan for
     `p_stat <= 1000 hPa` leaves the index at its `-1` sentinel when no level
     qualifies; the two other consumers of that index guard the sentinel, this one
     did not. Latent on Earth (p_surf ≈ 1013 hPa always fires by level 1), live on a
     column that sits entirely above 1000 hPa.
   - `AtomUtils::GetMean_2D/3D` built the global `m_node_weights` lazily with
     `clear()` + `push_back()`. `printDataAtm()` calls them from ~9 concurrent OpenMP
     sections, so the first call raced several threads through the same reallocation
     and corrupted the heap. It was masked in the parent only because an earlier
     single-threaded `GetMean_2D(temperature_NASA)` happened to run first.
   - `get_temperatures_from_curve()` dereferenced `begin()` and decremented `end()`
     *before* its own `size() < 2` guard — undefined behaviour on an empty map, which
     is what every call became once the Scotese curves were gone.

   The Makefile also gained `-MMD -MP` header dependencies; without them, edits to
   the headers where nearly all the physics lives did not trigger a rebuild.

   **What is not yet right.** The physics is still Earth physics on a 16 km shell.
   The column mean settles near 332 K and the prescribed 1500 K survives only at the
   domain top, because radiation still relaxes toward an Earth-like target, the
   saturation formula is still Magnus (capped ~101 °C), and the shell is ~20× too
   shallow for a 250 bar atmosphere. Deep convection is inactive: its trigger
   thresholds (1000/970/900/800 hPa) are absolute Earth surface pressures and never
   fire here. Phases 2–5 address these.

3. **Composition and mixture thermodynamics (done).** `MixtureAtm.h` derives the mass
   fractions, mean molar mass and gas constants from the configured mole fractions and
   supplies the *local* mixture properties `R_of(c, co2)`, `cp_of(c, co2, T)` (Shomate
   fits, 298–2000 K) and `M_of(c, co2)`. The constant `R_Air` / `cp_l` were replaced at
   the sites that set the pressure, density and Poisson exponent. CO₂ became a well-mixed
   mass fraction (ppm is meaningless at 20 % by mass) and water is initialised well-mixed
   rather than as a fraction of a saturation value that does not exist.

   Measured at the equator, initial state — every target hit exactly:

   | quantity | target | measured |
   |---|---|---|
   | surface pressure | 250 bar | 249.947 bar |
   | surface temperature | 1500 K | 1499.5 K |
   | surface density | 42.97 kg/m³ | 42.9695 kg/m³ |
   | mixture gas constant | 387.9 J/(kg·K) | 387.915 |
   | mean molar mass | 21.434 g/mol | 21.4337 |
   | H₂O mass fraction | 0.6724 | 0.6724 |

4. **Domain and non-dimensionalisation (done).** The shell went from 16 km to 300 km
   (`L_atm` is the stretch *amplitude*, so the shell is `(exp(zeta)−1)·L_atm`), `zeta`
   3.715 → 3.0 and `im` 41 → 61. The top cell now spans 1.65 local scale heights instead
   of 2.92. The column reaches **0.0237 bar** at the equator and 0.0051 bar at the pole,
   both below the 0.1 bar radiating level. `dt_visc` was scaled by the `dr²` CFL ratio,
   `abl_height` by the scale-height ratio, and the COSMO lapse parameter `beta` is now
   *derived* (`cosmo_lapse_fraction · R_mix · T_surf / cp`) rather than Earth's 42 K —
   which would have given a 0.71 K/km, essentially isothermal, 300 km column.

   Diagnostics added: a per-column profile and a per-level global summary
   (`ThermoAtm::printColumnProfile` / `printLevelSummary`), printed at init and every
   checkpoint. The level summary is what located two of the defects below.

   Five more inherited defects surfaced, all latent on Earth:

   - **`dr` was hard-coded to 0.025**, silently tied to `im = 41` (0.025 × 40 = 1). With
     61 levels the radial span became 1.5, one rad.z unit was read as 932 km instead of
     300, and the domain top landed at 1399 km. Now derived as `1/(im−1)`.
   - **`SaturationAdjustment` capped the prognostic temperature at 333.15 K and wrote the
     cap back** — "well above any physical surface temperature" is an Earth statement, and
     it collapsed 250 bar to 30 bar in one iteration. Replaced by a physical bound plus
     the correct statement: above 647.096 K water is supercritical and there is nothing to
     condense, so the adjustment is a genuine no-op there.
   - **The COSMO profile drives T → 0 at finite height** (305 km at the 1500 K equator,
     285 km at the 1450 K pole, since the coefficient goes as 1/T²). Temperature was
     floored at `t_00` but pressure kept falling at a rate its own temperature no longer
     justified. Now continues isothermally above the floor.
   - **`initCloudIce` manufactured cloud from a negative `q_sat`.** Magnus `E_sat` at
     1500 K (~1.2 × 10⁷ hPa) exceeds the 250 bar column, so `ep·E_sat/(p−E_sat)` goes
     negative and `c − H_crit·q_sat` *adds* water. It produced a condensate mass fraction
     of 0.47, inflating density by 1.87×. Same root cause drove the surface evaporation
     scheme to write a water mass fraction of 21.
   - **`bcRadius` extrapolated `p_stat`, `r_humid` and `r_dry` cubically at the lid.** The
     file already documents this stencil overshooting through zero for velocity and
     amplifying concavity for turbulence — but the hydrostatic quantities were left on it.
     On a 300 km shell it drove `p_stat` to −36 hPa across ~11 500 top-level cells. They
     now use log-linear extrapolation, which is exact for an isothermal layer and positive
     by construction.

   Also fixed: `co2Atmosphere()` ran *after* `densities()`, so the density was built with
   `R_of(c, 0) = 414.2` instead of 387.9 — a 7 % error through the whole column. Harmless
   on Earth, where CO₂ was ppm and never entered the density.

   **Stability.** A 20-iteration run at 8 threads completes cleanly: no NaN, no
   non-positive pressure at any of the 5 checkpoints, 34 ParaView files written. The
   equatorial surface settles rather than drifting away —

   | iter | T [K] | p [bar] | ρ [kg/m³] |
   |---|---|---|---|
   | 0  | 1499.5 | 249.947 | 42.9695 |
   | 5  | 1497.1 | 249.615 | 42.9222 |
   | 10 | 1496.0 | 249.360 | 42.8784 |
   | 15 | 1495.6 | 249.307 | 42.8692 |
   | 20 | 1495.5 | 249.280 | 42.8646 |

   with the per-checkpoint change decaying 2.4 → 1.1 → 0.4 → 0.1 K. The domain top rises
   from 0.024 to 0.049 bar over the run and stays below the radiating level.

   **What is not yet right.** Radiation is still the inherited Earth scheme (Bignami /
   Atwater–Ball emissivity, `radiation_mode = 5` relaxing toward a target that no longer
   exists), so the small residual drift above has no physical meaning — it is an Earth
   parameterisation being evaluated far outside its calibration, not a Hadean climate
   settling. Deep convection remains inactive: its trigger thresholds are absolute Earth
   surface pressures. Between 373 K and 647 K the saturation curve is still Magnus, far
   outside its validity. These are Phases 4 and 5.

5. **Saturation physics (done).** `SaturationH2O.h` replaces the Magnus formula
   everywhere — all ~25 call sites across 8 files — with the IAPWS correlations:

   - Wagner & Pruss (IAPWS-95) saturation line, 273.16 K to the critical point;
   - IAPWS (2011) ice-Ih sublimation curve below the triple point;
   - Watson latent heat, correctly vanishing at the critical point;
   - the **exact** mass-fraction conversion `q = x·M_H₂O / (x·M_H₂O + (1−x)·M_other)`,
     replacing the dilute `ep·E/(p−E)` that assumed water is a trace;
   - dew point by bisection (IAPWS has no closed-form inverse) and dq_sat/dT from
     Clausius–Clapeyron rather than the Magnus fit's own derivative.

   `make test` runs `test/saturation_selftest.cpp` against published reference points.
   All pass: the liquid curve to ~1e-5 relative (6.11657 hPa at 273.16 K, 1013.25 hPa at
   373.124 K, 220640 hPa at 647.095 K), sublimation to ~6e-3, plus monotonicity,
   [0,1]-boundedness, and agreement with the dilute form in the dilute limit.

   Writing the test first paid for itself: it caught a sign error in the sublimation
   formula (the summand is `a_i·θ^b_i`, not `a_i·(1−θ^b_i)` — the wrong form gives
   49 hPa instead of 0.76 hPa at 250 K) and two wrong expectations of my own.

   **Result: there is no condensation anywhere in the model.** The equatorial column runs
   1499.5 K at the surface to 676.7 K at 243 km and never crosses 647.096 K, so every
   level is supercritical — no cloud, no rain, no latent heat. The column printout now
   reports the saturation state per level and says so explicitly.

   That is physically coherent for a runaway greenhouse, but it is *not yet a result*: it
   follows from the prescribed lapse rate, and `cosmo_lapse_fraction = 0.5108` is still
   Earth calibration. At the dry adiabat (4.81 K/km) the column would cross 647 K near
   177 km and condensation would begin there. Which is right is for the radiation to
   decide, not for a prescribed lapse — see Phase 5.

6. **Radiation (done).** The Bignami/Atwater–Ball emissivity is replaced by a grey
   optical depth built from the absorber mass each layer actually holds, with pressure
   broadening:

   ```
   tau_i = SUM_s kappa_s * q_s * (dp_i/g) * (p_i/p_ref)
   ```

   The pressure-broadening factor is what no Earth-calibrated emissivity fit contains and
   what matters most here — at 250 bar it is a factor of 250 over the 1 bar reference.
   The old fit could not work in principle: `eps = 0.684 + 0.0056*e_surf` saturates to
   0.999 once `e_surf` exceeds ~56 hPa, and here it is ~2×10⁵ hPa, so every layer became
   a blackbody and the scheme carried no information about composition at all.

   Also changed: `radiation_mode` 5 → **2** (direct σT⁴ heating; modes 0/1/3/4/5 all lean
   on the Scotese snapshot or the 280 ppm CO₂ reference, neither of which exists at
   4.4 Ga); insolation rescaled to the faint young Sun (0.71 S₀); a **geothermal flux**
   added to the surface energy balance, since a molten surface supplies heat from below;
   albedo made a single cloud-deck value, Earth's ice/ocean contrast having no subject
   here; and the hard-coded `287.0` J/(kg·K) in the cloud-density calculation replaced by
   the local mixture gas constant.

   **The headline measurement.** Outgoing longwave flux at the top of the atmosphere:

   | | value |
   |---|---|
   | OLR | **324.5 W/m²** |
   | σ·T_surf⁴ | 286 629 W/m² |
   | suppression | **×883** |
   | transmitted fraction | 0.000 |

   and it stays at 324.5 W/m² at *every* checkpoint while σ·T_surf⁴ drifts from 286 629 to
   286 462. That near-total independence from surface temperature is the defining
   signature of a runaway greenhouse — the emission level sits inside the optically thick
   water column, not at the ground. The value is ~5 % above the canonical Nakajima /
   Komabayashi–Ingersoll limit of 280–310 W/m², which is close for a grey scheme whose
   three κ coefficients were chosen a priori rather than fitted.

   **What this does and does not establish.** It establishes that the opacity is the right
   order of magnitude and that the scheme reproduces the correct qualitative behaviour. It
   does not establish that 324.5 W/m² is the right number: `kappa_H2O`, `kappa_CO2` and
   `kappa_bg` carry roughly a factor-of-two uncertainty, and they are the biggest single
   lever on the answer. A grey scheme also cannot represent the window regions that set
   the real limit. Treat the value as a consistency check that passed, not as a
   prediction.

   Still open: the column remains supercritical throughout (Phase 4), so there is no
   cloud deck to justify the 0.4 albedo that is being assumed — the two are inconsistent,
   and resolving it needs the upper atmosphere to cool enough to condense.

*(Further entries are added as each phase is measured.)*
