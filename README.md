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
cd python && OMP_NUM_THREADS=8 ../cli/had config_athad.xml
```

Output lands in `python/output_Hadean/` as ParaView `.vtk` slices and `.vts`
panoramas — the ATOM line's convention (`output_<name>/`, relative to the run
directory). Run from `cli/` instead and it writes `cli/output_Hadean/`; the
`<output_path>` entry in the config controls it either way.

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

7. **Lapse rate and albedo made self-consistent (done).** This closed the inconsistency
   Phase 5 left standing — a 0.4 cloud albedo asserted over a column that condensed
   nothing.

   The root cause was `cosmo_lapse_fraction = 0.5108`, carried over from Earth. That
   number is not a tuning constant: Earth's ~5 K/km sits half way to its 9.8 K/km dry
   adiabat *because of latent heat release*. ATHAD's deep column condenses nothing, so
   there is nothing to flatten the lapse, and the consistent value is the dry adiabat.
   Carrying Earth's value across was self-fulfilling — the flattened lapse was justified
   by condensation that the flattening then prevented.

   Setting the fraction to 1.0 alone broke the model (surface NaN, 30 bar instead of 250),
   which exposed a deeper problem: the inherited COSMO profile is `T = T₀√(1−coeff·h)`, a
   **sqrt** in height. Matching its near-surface slope to the adiabat does not make it an
   adiabat — it plunges to zero at 156 km, inside the domain. So the fitted profile was
   replaced by the physics it was standing in for, integrated layer by layer:

   ```
   dry adiabat     dT/dz = -g/cp          (cp local: follows composition and T)
   hydrostatic     dp/dz = -p*g/(R*T)     (R local, on the layer-mean T)
   isothermal top  T = t_skin             where the adiabat falls below it
   ```

   No tuned constant, exact for a constant-cp adiabat, and it cannot produce the zero
   temperature the sqrt form did.

   Consequences measured:

   - The atmosphere is far more compressed on the true adiabat, so the domain came back
     from 300 km to **230 km**; 300 km put the top 90 km into near-vacuum (6e-8 bar).
   - **A cloud deck now forms**, from ~207 km (0.09 bar) upward, where p_H₂O finally
     exceeds p_sat. The top two levels report `CONDENSING`. The albedo is no longer
     asserted: `MultiLayerRadiation` uses the dark molten-surface value (0.08) and the
     existing condensate-driven cloud bump raises it where the model actually makes cloud.
     (The `albedo_pole`/`albedo_equator` parameters turned out to be **inert** — the
     radiation module builds its own albedo and never read them, so the earlier 0.4 was
     never in effect.)
   - `t_skin` is now derived from the **energy budget**, σT⁴ = (1−α)·SW + geothermal =
     236.0 W/m² → 254.0 K, replacing a value taken from a previously measured OLR, which
     was circular.

   **Energy balance closes:** OLR = 236.0 W/m², exactly the absorbed SW plus geothermal,
   with the surface suppressed by ×1214.

   **And that is the interesting result.** 236 W/m² is *below* the 280–310 W/m² runaway
   limit. At 0.71 S₀ with the assumed 150 W/m² geothermal flux, the planet does not absorb
   enough to sustain a runaway greenhouse — so the 1500 K surface is being held by fiat
   (it is prescribed), not by the budget. For a 1500 K surface to be self-consistently in
   runaway, the geothermal flux would have to be **≥ ~195 W/m²** (absorbed ≥ 280). That is
   not implausible for a genuine magma ocean, but it is a prediction the model now makes
   rather than an assumption it was given, and `geothermal_flux` should be revisited
   against magma-ocean cooling estimates.

   Still not closed: `t_skin` is a one-shot estimate using the clear-sky albedo, while the
   cloud deck raises it (albedo 0.4 would give 245.5 K). Making it a true fixed point means
   iterating `t_skin` against the model's own albedo. `initComposition()` prints both values
   and warns when they diverge.

8. **Rotation, paleo-time removal, Python bindings (done).**

   **Rotation.** ω = 3.17e-4 rad/s — a 5.5 h Hadean day, 4.35× modern, reported at
   startup. Checked against the family's four Coriolis/centrifugal sign fixes rather than
   re-derived. Two findings:

   - The **dynamical** Coriolis signs in `RHS_Atm_Turb.cpp` already agree with ATURAN
     `8b284cb` / ATNEPT `024c37f` once the opposite storage convention is accounted for
     (they store −a because their RHS subtracts; this file stores +a because its RHS adds).
     ATOM_Precipitation had fixed this independently. No change needed.
   - The **diagnostic** in `ThermoAtm::forces()` did not match. It carried a
     +2Ω·sinθ·u term the momentum equations drop under the traditional approximation, so
     the ParaView "Coriolis force" field showed a force that was never applied. It now
     follows the same `ATOM_CORIOLIS_NONTRAD` switch as the dynamics.
   - The centrifugal diagnostic used `(1 + |sinθ|)`, which is **maximal at the pole** where
     the true value is zero. Replaced by the correct decomposition about the rotation axis,
     `a_r = ω²r sin²θ`, `a_θ = ω²r sinθcosθ` (ATURAN `4201957` / ATJUP `8649675`).
     Diagnostic only — ATHAD's RHS carries no centrifugal term, the force being curl-free
     and absorbed by the pressure projection, which is what ATURAN found.

   **Paleo-time removal.** Nineteen parameters deleted (topography grids, NASA surface
   fields, Scotese curves, the pygplates reconstruction script, `Ma_switch`, `Ma_max`,
   `t_paleo_max`, `co2_paleo`, …), plus `read_Hydrosphere_SST()` and its two parameters —
   the atmosphere half of an atm↔ocean Picard loop with no ocean to couple to.
   `time_start/end/step` remain: the time-slice loop is still structural, and a single
   slice is all that runs.

   **Python bindings.** `pyatom` → `pyathad`, hydrosphere halves stripped. `model.py` was
   driving *both* spheres and — at the bottom — actually running the ocean rather than the
   atmosphere; it now drives the one model that exists.

   **Diagnostic cadence** (separate from `checkpoint`, which writes ParaView files):
   every **10** iterations for a short run (`nm ≤ 100`), every **100** for a longer one,
   overridable with `diagnostic_stride`. The chosen value is reported at startup.
   Verified: nm=20 prints at iters 10 and 20; nm=400 selects the 100 stride.

   A 20-iteration run completes clean — no NaN, no non-positive pressure, OLR steady at
   236.0 W/m².

9. **The saturation conversion finished, and 60 km of manufactured cloud removed (done).**

   Phase 5 replaced Magnus with IAPWS everywhere, which fixed the saturation *pressure*.
   It did not finish the *conversion* from that pressure to a mass fraction. About twenty
   sites still carried the dilute form `q_sat = ep·E/(p − E)` that invariant 2 in
   CLAUDE.md forbids, and the equatorial column printout was showing the consequence
   without anyone reading it that way: water vapour oscillating between 0 and 0.88 above
   140 km, and a density that *rose* with height at 142 km.

   Three distinct defects, all latent on Earth:

   - **`SaturationAdjustment::clampAndFade` had the superheated branch inverted.** It read
     `q_sat = (p > E_sat) ? ep·E_sat/(p − E_sat) : ep·1e-5`. When p_sat(T) exceeds the
     local pressure the vapour is superheated and **cannot condense at all**, so the
     saturation limit is 1 — all of the water stays vapour. The fallback asserted the
     opposite, a limit of ~7e-6, and the block below it therefore dumped the **entire
     0.67 water mass fraction into cloud** and released its latent heat, in exactly the
     layers where nothing can condense. On this column that is every level between ~373 K
     and the critical point: a **60 km slab of manufactured cloud from ~140 to ~200 km**,
     which set the planetary albedo and, smeared downward by `damp_wiggles`, inflated the
     density in the supercritical layers beneath it through the `(1 − cloud − ice)`
     loading term. Latent on Earth because no terrestrial cell is ever above 373 K, so the
     branch never ran — the same shape of defect as the 333.15 K temperature cap.
   - **The fix from Phase 5 was applied to the entry `q_sat` but not to the copy inside
     the same function's Newton loop**, which went on pulling q_v toward `ep·1e-5`.
   - **`AtmMixture::M_nonwater` could not do what its name and comment claim.** It took
     `(co2, M_background)` and set `q_b = 1 − q_c`, so its renormalisation "to exclude
     H₂O" was identically 1 and the water's share of the mass was silently handed to the
     background gas. It returned 28.58 g/mol where the reference composition gives 35.11 —
     the carrier 19 % too light, hence a q_sat some 23 % too large in every exact
     conversion in the model. Invisible on Earth, where water is ~1 % of the mass. It now
     takes the local water mass fraction as well.

   All remaining dilute sites went with them: `initCloudIce` (guarded above the critical
   point but not in the subcritical superheated band, where its q_sat goes negative and
   the `c − H_crit·q_sat` subtraction manufactures cloud), all four ice schemes,
   `RHS_Atm_Turb`'s `q_Rain`/`q_Ice` thresholds — where a negative threshold switches the
   latent-heat term on unconditionally — and the dead `init_vapour_cloud`, which is one
   uncomment away from being live.

   Diagnostics added, because the defect was visible for weeks in output nobody could
   read: the column profile now prints `q_cld` and `q_ice` and flags a level `CLOUD?!`
   when it carries condensate that cannot exist there, and a new
   `ThermoAtm::printPlanetaryBalance` prints the cos(latitude)-weighted mean albedo,
   absorbed shortwave, OLR and the implied skin temperature every diagnostic step.

   **Measured, 20 iterations, 8 threads, clean (no NaN, pressure positive everywhere):**

   | | before | after |
   |---|---|---|
   | condensate at 142–200 km | up to 0.42 | ≤ 0.05, and falling |
   | q_H₂O at 186 km | 0.0000 | 0.9177 |
   | T at 218 km | 254.0 K | 427.3 K |
   | domain top | 0.028 bar | 0.295 bar |
   | mean planetary albedo | — (not measured) | 0.4999 |

   The upper column is **much warmer**, and that is the point: with the bogus condensation
   gone it stays steam, keeps steam's high cp, and cools along a shallower adiabat. The
   234 km shell no longer reaches the radiating level.

   **Two things this breaks, both now open:**

   - **The shell is too shallow again.** The top is at 0.295 bar, above the 0.1 bar
     radiating level. The 300 → 230 km reduction in Phase 7 was made against the profile
     the manufactured latent heat produced, so it has to be re-derived.
   - **The planetary albedo is pinned at 0.4999**, which is the hard-coded thick-cloud
     asymptote `alpha_cloud = 0.50` in `MultiLayerRadiation`. Once any condensate is
     present the reflectivity `tau/(tau+2)` saturates, so the model is not *computing* an
     albedo, it is *reporting a constant* — and that constant is a bare literal inside a
     physics kernel, not a parameter. It has to be exposed and justified before the OLR or
     the runaway claim that depends on it means anything.

   With that albedo the budget no longer closes: absorbed SW + geothermal = 203.8 W/m²
   against an OLR of 236.0, an imbalance of **−32.3 W/m²**, because `t_skin` is still the
   configured 254 K while the energy balance now implies 244.8 K. That is the fixed point
   below, and it is no longer optional.

10. **The insolation, the albedo parameters, the t_skin fixed point — and what they
    exposed: the OLR is not an output (done, and it is bad news).**

    Three fixes, and then the thing they uncovered.

    **The insolation was Earth's surface flux used as top-of-atmosphere insolation.**
    `rad_equator_short` / `rad_pole_short` were 116 / 71 W/m², which is
    ATOM_Precipitation's 163.3 / 100.0 scaled by 0.71 — and those are Earth's
    absorbed-at-the-*surface* shortwave, already reduced by Earth's albedo and its
    atmosphere's absorption. `MultiLayerRadiation` uses them as the flux incident at the
    top and then applies ATHAD's own albedo, counting Earth's albedo twice. The
    cos(latitude)-weighted mean was **107.5 W/m² where 0.71 S₀ delivers 241.6**: the planet
    was lit at 45 % of its own insolation. Replaced by a fit of the model's parabola to the
    annual-mean insolation of a zero-obliquity planet, `Q(φ) = S/π·cos φ`, constrained so
    the weighted mean is exactly S/4 — 298.0 at the equator, 0.0 at the pole, RMS 6 W/m².
    Measured mean after the change: **241.56 W/m²**. (Hadean obliquity is unknown. Zero is
    the assumption and also the only shape the parabola can follow; at Earth's 23.44° the
    true curve flattens toward the pole in a way it cannot.)

    **The albedo constants became parameters.** `alpha_cloud = 0.50` and the molten-surface
    0.08 were bare literals in `MultiLayerRadiation` while the config carried an
    `albedo_pole`/`albedo_equator` pair that nothing read. They are now `albedo_cloud` and
    `albedo_surface`, and the inert pair is gone. This matters more than tidiness:
    the cloud deck's condensate path is ~10⁵ g/m² against a `cwp_tau` of 100, so
    `refl = τ/(τ+2)` saturates and **every cloudy column returns `albedo_cloud` to four
    decimals**. The measured mean planetary albedo is 0.4999. The model is not computing an
    albedo; it is reporting that constant.

    **`t_skin` is now an iterated fixed point** against the model's own albedo,
    `σ·t_skin⁴ = (1−ᾱ)·S̄ + F_geo`, relaxed by `t_skin_relax` (0.25) after each radiation
    call. It converges from the configured 254 K to **262.85 K** against a target of
    262.88 K within 20 iterations. Setting `t_skin_relax = 0` restores the old fixed value.

    **What all of that exposed.** The model's outgoing longwave flux is not a computed
    quantity. It is σT⁴ of the topmost layer, which is prescribed.

    | run | lid T | lid ε | reported OLR | σ·T_lid⁴ |
    |---|---|---|---|---|
    | `t_skin` = 254 K | 254.0 K | 1.0000 | 236.01 W/m² | 236.01 |
    | `t_skin` = 240 K | 240.0 K | 1.0000 | 188.13 W/m² | 188.13 |
    | after these fixes | 343.98 K | 1.0000 | 802.32 W/m² | 793.84 |

    The first two rows are a deliberate experiment: changing one configured number moved
    the "measured" OLR to σ times that number to five significant figures. The lid is
    optically thick (τ ≈ 6 at 0.29 bar), so the emission level sits on the domain boundary
    and the model reports the boundary condition back.

    **The 236.0 W/m² of Phase 7, and the "energy balance closes" that went with it, was
    this.** It was worse than a coincidence: `bcRadius` pins the lid temperature to a
    snapshot `t_top_init` taken once from the initial condition, and that snapshot came
    from `initTemperatureData`, which builds its adiabat *before* the water and CO₂ fields
    exist and therefore with a background-only cp and a much steeper lapse. `densities()`
    rebuilt the whole column on the correct mixture cp every iteration — and could never
    move the lid, because bcRadius pinned it straight back to the stale value. The lid
    stayed at 254.0 K, the OLR stayed at σ·254⁴, and the agreement with the energy budget
    was the budget being read back out of the number it had been used to set. `densities()`
    now refreshes the snapshot, which is what makes the 343.98 K above visible.

    **So the honest current state is a −531 W/m² imbalance and no meaningful OLR**, and
    that is a better description of the model than the closed budget it replaces.

    **The shell cannot currently be deepened to fix it.** 230 km tops out at 0.29 bar,
    three times above the 0.1 bar radiating level. Going deeper — 260 km (top 0.017 bar)
    and 300 km (top 3.2e-4 bar) — produces a finite, sane initial state and then **NaN
    across the entire field in the first `MultiLayerRadiation` call**: the field is finite
    in the diagnostic immediately before it and NaN in the one immediately after, in both
    runs. The scheme's tridiagonal (Thomas) assembly is built from products of the layer
    emissivity and its rows degenerate as ε → 0 — which is exactly the condition at the top
    of any domain that reaches the radiating level. The radiation has to be reformulated
    before the shell can grow, so 230 km stays.

    Also fixed in passing: `python/PythonStream.cpp` still included `pyatom.h`, so the
    default `make` target had not built the Python bindings since the fork. It does now.

## Remaining work

**The next task is the radiation scheme**, and everything else waits on it. Two coupled
defects:

- **The tridiagonal solve degenerates for optically thin layers** (item 10), which is why
  the domain cannot reach the radiating level. A formulation that stays well conditioned
  as ε → 0 — a direct up/down flux integration rather than the Thomas inversion — would
  remove the ceiling.
- **The temperature profile is prescribed, not solved.** `ThermoAtm::densities()`
  re-imposes the adiabat + isothermal top on `t` every iteration, so whatever the dynamics
  and the radiation compute is overwritten before it can matter. Invariant 3 in CLAUDE.md
  says radiation must *set* the profile rather than nudge it toward a prescribed one; in
  the current code the prescription wins. Until it does not, the OLR is an input however
  deep the shell is, because the model will always emit σT⁴ at whichever prescribed level
  first becomes opaque.

After that:


- **`albedo_cloud = 0.50` IS the model's planetary albedo** and is an assumption. It is
  now the second-biggest lever after the opacities. A deep, cold, slowly sedimenting
  Hadean deck could plausibly be brighter.
- **`geothermal_flux` is the open number.** The ≥ ~195 W/m² figure was derived with the
  clear-sky albedo and the too-low insolation, so it has to be redone once the radiation
  can produce an OLR worth comparing to. Check against magma-ocean cooling estimates.
- **`initTemperatureData` builds its column before the composition exists** — no water, no
  CO₂, so a background-only cp and too steep a lapse. `densities()` overwrites it, so the
  only thing that ever escaped was the lid snapshot (item 10), which is now refreshed. It
  should still be built on the real mixture.
- **`initCloudIce`'s H_crit parabola is keyed to absolute pressure** (`p_crit = 1000` hPa,
  `p_mid = 550`), an Earth surface pressure. It should be a fraction of the local surface
  pressure, like the deep-convection triggers.
- **The three κ opacities** carry factor-of-two uncertainty and are the biggest lever on
  the OLR. A grey scheme cannot represent the window regions that set the real limit.
- **The surface temperature is prescribed, not solved.** Everything above is conditional
  on that.
- **Boussinesq** remains untested against a column whose density spans two orders of
  magnitude (see CLAUDE.md).
