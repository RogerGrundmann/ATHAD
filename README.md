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

11. **The radiation scheme rewritten as two-stream flux sweeps; the shell deepened to
    300 km; the OLR is a computed quantity for the first time (done).**

    The inherited scheme assembled a tridiagonal system whose every entry was a product of
    the layer emissivity — sub-diagonal `ε_{i-1}σT_{i-1}⁴`, diagonal `−2ε_iσT_i⁴`,
    super-diagonal `ε_{i+1}σT_{i+1}⁴` — with a right-hand side built from differences of
    cumulative transmitted sums, solved by Thomas for a correction to `σT⁴`. Optical depth
    goes as p² through the pressure broadening, so near the top adjacent rows differ by
    orders of magnitude while the right-hand side is a difference of two nearly equal
    sums. It did not survive a domain that reaches the radiating level (item 10).

    Replaced by the plane-parallel grey transfer it was standing in for, integrated
    directly:

    ```
    up[i] = up[i-1]·(1 − ε_i) + ε_i·σT_i⁴        surface → top,  up[0] = σT_surf⁴
    dn[i] = dn[i+1]·(1 − ε_i) + ε_i·σT_i⁴        top → surface,  nothing enters from space
    ```

    Radiative equilibrium of a layer is then the statement that it absorbs what it emits,
    `ε_i(up[i-1] + dn[i+1]) = 2ε_iσT_i⁴`, **and the emissivity cancels**:

    ```
    σT_i⁴ = (up[i-1] + dn[i+1]) / 2
    ```

    Nothing divides by ε. The ε → 0 limit is exactly right rather than merely survivable —
    a transparent layer passes both streams and takes the mean of what goes by — and at
    the top, where `dn → 0`, it reduces to `σT⁴ = up/2`, the classical skin temperature
    that this model has until now been *prescribing* as `t_skin`. The temperature and the
    fluxes are iterated against each other (Lambda iteration, 4 passes); convergence is
    slow in optically thick layers, the known weakness of the method, but there the update
    is nearly a no-op anyway, and the fluxes are exact for the current profile at every
    pass.

    **Measured, 20 iterations, 8 threads, no NaN, pressure positive everywhere:**

    | shell | top p | lid ε | OLR | σ·T_lid⁴ | |
    |---|---|---|---|---|---|
    | 230 km | 0.29 bar | 1.0000 | 795 W/m² | 787 W/m² | OLR *is* the lid — an input |
    | 260 km | 0.017 bar | 0.0610 | 519 W/m² | 271 W/m² | decoupled |
    | 300 km | 3.8e-4 bar | 0.0000 | 581 W/m² | 271 W/m² | **a real column integral** |

    The 230 km row is the regression check: where the old scheme worked, the new one
    reproduces it (795 against the old 802 W/m²). The other two rows are runs that
    previously NaN'd in the first radiation call.

    **`L_atm` is now 300 km.** The column tops out at 3.8e-4 bar with a transparent lid,
    the isothermal skin is resolved from 256 km up, condensation begins at 242.8 km, and
    the outgoing flux is no longer the boundary temperature read back out.

    **And the first thing the model says with it is that its own opacity is too low.**
    OLR = 581 W/m² against 271 W/m² absorbed plus geothermal: the atmosphere radiates away
    more than twice what it takes in, so it cannot hold the prescribed 1500 K surface. That
    is now a statement about `kappa_H2O` = 0.01 m²/kg rather than an artefact of the
    boundary. Note also that the OLR is **not yet grid-converged** — 519 W/m² at 260 km
    against 581 at 300 km, because `im` is fixed at 61 and a deeper shell is a coarser one.

    Two consequential changes came with it:

    - **`radiation.x` is now the upward long-wave flux at the top of each layer**, not
      `σT⁴` of that layer. It is a diagnostic field (nothing feeds it back into the
      dynamics), it is continuous across the surface by construction — so the 1-2-1
      smoothing pass that used to hide the surface kink is gone — and `radiation.x[im-1]`
      is the OLR, which is what the mode-5 cloud diagnostic already assumed it was.
    - **`bcRadius` no longer pins the radiation lid** to `σ·(t·t_0)⁴`. That pin was
      justified by "radiation.x = σ·(t·t_0)⁴", which has stopped being true; keeping it
      would have thrown away the one number the column integration exists to produce.

    **A separate defect, found by asking whether the vertical stretch reached everything:**
    `init_tropopause_layers` converted a height to a level index with `round(h / L_atm)`.
    That is a level index only for uniformly spaced layers, and the grid is exponentially
    stretched — `height(i) = (exp(zeta·i/(im−1)) − 1)·L_atm`, so `L_atm` is the *amplitude*
    of the stretch, not a spacing. The exact inverse is
    `i = (im−1)·ln(1 + h/L_atm)/zeta`. At 300 km the pole's 195 km convective top sits at
    level **52**; the old formula returned **12**, which is 13 km. `VelocityInitializer`
    builds its entire jet profile between the surface and that level and applies a linear
    taper to zero from it up to the domain top, so the initial wind structure was
    compressed into the bottom 4 % of the atmosphere and the remaining 96 % got the taper.
    Wrong on Earth too (28 of 41 levels = 4.9 km, not the intended 11 km), but wrong there
    in a way that still landed inside the troposphere, so it never showed.

12. **CO₂ made genuinely prognostic, and a 100-iteration run (done).**

    The CO₂ in the ParaView output was *exactly* constant — `min co2 = max co2 = 0.205300`
    in every cell at every checkpoint. `ThermoAtm::co2Atmosphere()` fills the whole field
    with the uniform `co2_0·co2_scale`, and it was being called **inside the time loop**,
    immediately before `densities()`. Meanwhile the model does carry a full CO₂ transport
    equation — `RHS_Atm_Turb` builds `rhs_co2`, `RungeKutta_Atm_Turb` carries it through
    all four stages — and every iteration the result was discarded. CLAUDE.md's "H₂O and
    CO₂ are prognostic" was false for CO₂.

    `co2Atmosphere()` is now the initial condition only. In its place in the loop,
    `ThermoAtm::co2Column()`:

    - fills `co2_total`, which was declared ("areas of higher co2 concentration") and
      **never written**, so its min/max and the `co2_average` derived from it both read
      0.000. It now carries the column CO₂ mass path in kg/m², the CO₂ analogue of
      precipitable water: 519 108 kg/m² global mean, 505 205 (pole) to 522 343 (equator);
    - reports the drift of the global mass-weighted mean q_CO₂. There is no CO₂ source or
      sink anywhere in the model, so that mean is conserved and any drift is transport
      error. Measured over 20 iterations: **−0.0000 %**.

    Unit labels corrected with it: the 3-D field is kg/kg not ppm, the column is kg/m², and
    the startup banner's `co2_0=0.205 ppm` is kg/kg.

    **The field still comes out uniform — and that is now the answer rather than the
    assumption.** A passive tracer with no sources and no gradients has `∇q = 0`, so
    advection and diffusion both vanish and uniform is the exact solution. The difference
    is that the model now computes it, monitors it, and will transport any structure that
    does appear instead of erasing it every step.

    **100-iteration run**, 300 km shell, 24 threads, no NaN, ParaView every 10 iterations.
    The model relaxes monotonically toward radiative balance:

    | iter | 10 | 20 | 40 | 60 | 80 | 100 |
    |---|---|---|---|---|---|---|
    | OLR [W/m²] | 617.8 | 581.0 | 522.0 | 477.0 | 441.7 | 413.4 |
    | imbalance [W/m²] | −347.0 | −310.2 | −251.2 | −206.2 | −170.9 | −142.6 |

    The imbalance decays by about 9 % per 10 iterations and the rate is slowing, so
    reaching balance (OLR → 271 W/m²) needs of order 400 iterations. Note the OLR passes
    *through* the 280–310 W/m² Nakajima / Komabayashi–Ingersoll band on the way down —
    which is the first time this model has approached that limit from a computed flux
    rather than a prescribed one. `t_skin` has converged to 262.88 K and the lid emissivity
    is 0.0000, so the outgoing flux is a column integral throughout.

13. **The Hadley cells were not symmetric, and nothing here can make them so (done).**

    Spotted in the iteration-0 ParaView output. `VelocityInitializer::compute()` set the
    meridional wind at 15°N to a surface coefficient of **4.0** and at 15°S to **3.0**:

    ```cpp
    init_v_or_w(m.v,  75, -3.0,  4.0);   // lat:  15N
    init_v_or_w(m.v, 105, -3.0,  3.0);   // lat:  15S
    ```

    Every other mirror pair in the whole u/v/w initialisation is identical — 0/180, 15/165,
    30/150, 45/135, 60/120 — and the `form_diagonals` spans are mirrored too, so this was
    the only asymmetry in the velocity initial condition. The commented-out lines that sat
    directly above each showed the values had once been the other way round (3.0 north,
    4.0 south), so it had never been symmetric; someone had swapped which hemisphere won.

    Measured in `meridional_streamfunction_10.csv`, at the level of maximum |Ψ|:

    | lat | Ψ(+lat) | Ψ(−lat) | sum (0 if symmetric) |
    |---|---|---|---|
    | 30° | 74 187 | −72 049 | 2 138 |
    | 15° | −278 138 | 366 213 | **88 075** |
    | 10° | −198 131 | 260 698 | 62 567 |

    (10⁹ kg/s.) The southern cell is stronger than the northern by 366/278 = **1.32**
    against the coefficients' 4.0/3.0 = 1.33, and the global antisymmetry error is 11.4 %.
    At 30°, away from the injected asymmetry, the mismatch is 2.9 %. It does not wash out:
    still 8.3 % at iteration 100.

    On Earth this asymmetry is physical — the ITCZ sits north of the equator because of the
    land–sea distribution. **ATHAD cannot have it.** There is no land, no topography, the
    prescribed surface temperature is a symmetric parabola, the insolation profile is
    explicitly mirrored (`short_wave_radiation[j] = short_wave_radiation[j_max-j]`), and
    there is no obliquity and no seasonal cycle. Nothing in the model can sustain a
    hemispheric asymmetry, so all of it was inherited from these two numbers.

    Both are now **3.5**, their mean, which removes the asymmetry and leaves the total
    initial Hadley mass flux unchanged. Measured after the fix, same diagnostic:

    | lat | Ψ(+lat) | Ψ(−lat) | sum |
    |---|---|---|---|
    | 30° | 73 094 | −73 138 | −44 |
    | 15° | −322 169 | 322 181 | 12 |
    | 10° | −229 399 | 229 418 | 19 |

    **Global antisymmetry error 11.38 % → 0.0152 %**, the residue being floating-point and
    OpenMP reduction noise rather than structure.

14. **No condensed phase where none can exist (done) — and what enforcing it exposed.**

    The cloud deck in the ParaView output sat at 166–230 km, and every level of it was
    already flagged `CLOUD?!` by the model's own column diagnostic: condensate in cells
    whose `q_sat` is 1.0000, which is what "no saturation limit" means. The genuine
    condensation level is 242.8 km.

    Cause: **none of the four ice schemes, nor MoistConvection, contained a single
    critical-point check** — `grep T_CRIT` across all five files returned nothing.
    `TwoCatIceScheme` ran its full warm and cold microphysics from the ground up, and the
    precipitation maxima landed at i=50 (175.8 km) and i=54 (218.2 km), inside the band.
    Condensate formed legitimately at 243 km, sedimented into air at 350–700 K where it
    should flash to vapour instantly, and the scheme kept it and shuffled it between rain,
    snow, cloud and ice. Latent on Earth, where no cell is near 647 K or above its own
    saturation pressure, so the question never arises.

    `IceSchemeCommon::canCondense()` now states the two conditions — supercritical
    (T ≥ 647.096 K, where liquid and vapour are one phase) and superheated (p_sat(T) > p,
    where the vapour cannot reach saturation whatever its abundance) — and
    `evaporateWhereImpossible()` enforces them: condensate goes back to the vapour with its
    latent heat (no latent heat above the critical point, where there was never a separate
    phase), and the sources and precipitation fluxes are cleared. It is called from all
    four ice schemes and from `SaturationAdjustment::clampAndFade`, whose "supercritical
    cells carry no condensate and need no fade" `continue` was exactly the wrong response —
    they carry none only if something removes it, and nothing did.

    **Result: zero `CLOUD?!` levels anywhere in the run**, from every level between 166 and
    230 km before.

    **What it exposed.** With the condensate returned to the vapour, `q_H2O` in the
    sub-cloud band went to **1.27** — a mass fraction above one. It turned out the field had
    been unphysical all along, just less visibly: before the fix it sat at 0.9971 against
    `co2` = 0.2053, a composition summing to 1.20. Nothing bounded it. The three mass
    fractions must sum to one and the background is carried as the remainder `1 − c − co2`,
    so `c > 1 − co2` is a **negative background mass** — and `AtmMixture::split()`
    renormalises defensively, so the only symptom was a gas constant pinned at 415.1 instead
    of moving with the composition. `UtilsAtm::valueLimitationAtm` carries such a bound at
    the wrong ceiling of 1, and is commented out at both of its call sites.

    The ceiling `c ≤ 1 − co2` is now enforced every iteration, outside the moist-physics
    block (which runs only on moist iterations, while the RK4 transport of `c` runs on all
    of them). The column is physical again: `c` = 0.7947 through the band, R = 405.6.

    **But the ceiling bites in 282 302 cells, and that is an open problem, not a fix.**
    Water is pumped downward out of the single condensing level by sedimentation and
    evaporates into the superheated band with no return path, so `c` there grows until the
    clamp stops it — which deletes water. The count is printed at every diagnostic step.
    A run in which it stays large is not to be trusted, and this one does.

15. **The water is not lost, it is created — and the moist physics was not running (done:
    diagnosis and instrumentation).**

    ⚠️ **The headline of this item is wrong, and item 17 corrects it.** The water was not
    being created: the mass-weighted mean divides by an air mass that `densities()` rebuilds
    every iteration, and it was the air mass that moved. Measured against frozen weights the
    water field drifts +0.0011 % where this item reports +0.17 %. The instrumentation and
    the two experiments below stand; the attribution to the transport does not. Read this
    item for how the measurement was built and item 17 for what it turned out to measure.

    Item 14 left the `c ≤ 1 − co2` ceiling deleting water in 282 302 cells per iteration,
    with the working hypothesis that sedimentation pumps water out of the condensing level
    and nothing returns it. That hypothesis is wrong.

    `ThermoAtm::waterBudget()` now reports the global mass-weighted mean of
    vapour + cloud + ice + graupel, against its initial value. ATHAD has no water source
    and no water sink — the surface is supercritical, so `waterVapourEvaporation()` returns
    at once — so that mean is conserved exactly, and any drift is scheme error. Nothing had
    ever measured it. Measured over 20 iterations:

    ```
    iter  0:  q = 0.672751   drift +0.0000 %   deleted by the ceiling 0.000000
    iter 10:  q = 0.673903   drift +0.1712 %   deleted by the ceiling 0.007173
    iter 20:  q = 0.673893   drift +0.1698 %   deleted by the ceiling 0.014293
    ```

    Water is **created** at ~0.0007 kg/kg per iteration — about 0.11 % of the total per
    iteration — and the ceiling removes almost exactly as much as appears. The clamp is not
    starving the band; it is bailing out a leak.

    Two experiments locate it:

    - **`CategoryIceScheme = -1`** (no ice scheme at all): identical, +0.1718 % against
      +0.1712 %. Not the microphysics.
    - **`moist_phys_start_iter = 300`.** The moist physics — SaturationAdjustment, the ice
      schemes, MoistConvection — **does not run at all until iteration 300.** Every
      20-iteration diagnostic in this README was therefore a *dry* run, and the water is
      created with no moist physics executing. What remains is the RK4 tracer transport.

    **And the CO₂ conservation of item 12 does not contradict this — it is worthless as a
    test.** CO₂ drifts +0.0000 % because it is *uniform*, and a uniform tracer is conserved
    by any consistent advection scheme, however non-conservative, since `u·∇q = 0`. Water
    carries structure (put there by the one-off moist physics in the initialisation), and
    only a field with structure can expose the error. The earlier claim that the CO₂ result
    validated the transport was wrong.

    **Where the drift is.** Per-level attribution in `waterBudget()` puts it at levels
    42–47, **113 to 149 km** — deep in the supercritical column, far below both the cloud
    deck at 243 km and the band the `c` ceiling clamps at 185–230 km. It is static between
    iterations 10 and 20 (+2.83e-3 against +2.82e-3), so the redistribution happens early
    and then settles into a steady state in which creation balances the ceiling's deletion.
    Water moved *downward* out of the low-density band into air one to two orders of
    magnitude denser, and the mass-weighted total rose accordingly.

    **The cause is that the flow does not satisfy the continuity the tracer equation
    assumes.** `RHS_Atm_Turb` integrates `∂q/∂t = −u·∇q + …`, which conserves `∫ρq dV` only
    when `∇·(ρu) = 0`. The pressure projection enforces `∇·u = 0`, and ATHAD's density spans
    five orders of magnitude across the column.

    **The flux-form correction proposed here first is the wrong fix, and was not applied.**
    Expanding `∂(ρq)/∂t + ∇·(ρuq) = 0` gives `q[∂ρ/∂t + ∇·(ρu)] + ρ[∂q/∂t + u·∇q] = 0`, and
    the first bracket *is* continuity — so for a **mass fraction** the advective form is
    already exactly right, conditional on continuity. Adding `−(q/ρ)·∇·(ρu)` against a fixed
    `ρ` and a velocity field with `∇·(ρu) ≠ 0` would restore the global integral by making a
    *uniform* tracer develop structure: CO₂ would stop being well mixed, which is precisely
    the failure the uniformity test in the scope below exists to catch. It would trade a
    measured 0.17 % global error for an unphysical local one.

    So the defect is upstream, in the flow, and the fix belongs in `PressureSolverAtm`:
    anelastic continuity, `∇·(ρu) = 0`. This is the **Boussinesq risk in CLAUDE.md arriving
    in concrete form**. It changes the dynamical core and every result in this README, which
    is why it is scoped below rather than applied.

    Two further consequences worth knowing: the cloud deck discussed in items 9 and 14 is
    an *initialisation-time* state that the first 300 iterations only advect, and the
    ParaView output of a 400-iteration run has active moist physics only over its last 100
    iterations.

16. **The cloud came back: an undamped Newton step, and precipitation mass returned
    (done). And step 1 of the anelastic scope, measured and falsified.**

    The output showed water vapour everywhere and **no cloud, ice or rain anywhere in the
    domain** — while the column at 243 km sat ninefold supersaturated, `q_H2O` = 0.72
    against `q_sat` = 0.076, flagged `CONDENSING`. Two separate causes.

    **The precipitation mass was being destroyed.** `evaporateWhereImpossible()` (item 14)
    zeroed `P_rain` and `P_snow` without returning their mass to the vapour. The deck
    condensed, the ice scheme autoconverted it to rain, the rain fell one level into the
    superheated band, and it vanished. A downward flux `P` falling at speed `v` corresponds
    to a mass fraction `P/(v·ρ)` in the air it passes through; that is now converted before
    the flux is cleared, using the fall speeds the schemes' own residence times are built
    from (1.6 m/s rain, 0.96 m/s snow).

    **The first Newton step of the saturation adjustment was undamped**, and that is what
    actually suppressed the cloud. `adjustSaturation` initialised `q_v_hyp = q_sat`, a jump
    straight to the saturation value, and applied its `omega = 1/(1+Gain)` damping only
    from the second pass. On Earth this is harmless — condensing the ~0.01 kg/kg a
    terrestrial parcel holds releases about 12 K. Here the undamped step condenses
    **0.65 kg/kg in one go and releases 800 K** of latent heat. `T` is then clamped at the
    critical point, where `p_sat` = 220 bar against a local 0.085 bar, so `q_sat` becomes 1,
    the target inverts, and the next pass evaporates everything back. The iteration
    flip-flops between fully condensed and fully evaporated and finishes at zero.

    This one is not an inherited constant, it is an inherited **assumption**: that latent
    heating is a perturbation. At 67 % water by mass, condensation is a bulk phase change of
    the atmosphere. Starting the iteration from `q_v_b` makes the first pass compute a
    properly damped target. The equilibrium it should find is modest — condensing ~0.04
    kg/kg warms the parcel ~47 K, after which the vapour is superheated and nothing more
    can condense.

    **Measured after the fix** (20 iterations, moist physics on from iteration 1):

    ```
    max cloud water = 15.687 g/kg  @ 73°S, 242 774 m
    max cloud ice   = 14.400 g/kg  @ 68°S, 256 027 m
    ```

    A cloud deck at the condensation level, in both hemispheres, at the high latitudes
    where the column is coldest — and none anywhere it cannot exist.

    **Step 1 of the anelastic scope, measured.** `PressureSolverAtm` now reports the
    divergence of the *actual* velocity after the projection — nothing had ever measured
    it, so "the projection enforces ∇·u = 0" was an assumption about the code rather than
    an observation of it. A/B on the existing `ATM_POISSON_METRIC_FIX` knob:

    | | `div(u)` rms | water drift, 20 iters | ceiling deletions |
    |---|---|---|---|
    | metric fix off | 2.625e-02 | +0.1698 % | 0.014293 |
    | metric fix on | **2.153e-02** | +0.1698 % | 0.014293 |

    The consistent metric reduces the residual divergence by 18 % — real, and worth
    keeping — but the water drift is **bit-identical**. So the metric inconsistency is not
    the cause of the mass error, and step 1 of the scope is falsified as an explanation
    while remaining valid as a repair. Note also the absolute number: an rms `∇·u` of
    2.6e-02 is not a small residual. The projection is leaving a great deal of divergence
    behind, which is consistent with the anelastic diagnosis and makes steps 2–7 the
    remaining candidate.

17. **The anelastic projection, built and measured — and the water was never being created
    (done).**

    Steps 2–7 of the scope below are implemented behind `ATM_ANELASTIC` (default 0,
    bit-identical when unset): a one-dimensional base state `ρ̄(z)` rebuilt by
    `ThermoAtm::densities()` alongside the profile it averages; the divergence source
    `D = ∇·u* + u*_r·dlnρ̄/dr`; the same `dlnρ̄/dr` as a first-derivative term in the Poisson
    stencil, so source and operator stay adjoint — the lesson of step 1; and `ρ̄u_r = 0` at
    the surface and the lid in place of the `c43/c13` extrapolation, which permitted a
    through-wall mass flux. Step 7 needed no change: `t_ref_level[i]` is already the same
    horizontal mean, of the same prescribed profile, that `ρ̄` is.

    **It works, and it does not fix the water.** A/B at 20 iterations, metric fix on in
    both, `ATM_ANELASTIC` the only difference:

    | | off | on |
    |---|---|---|
    | `∇·u` rms | 2.153e-02 | 2.670e-02 |
    | `∇·(ρ̄u)/ρ̄` rms | 3.359e-02 | **2.607e-02** |
    | max radial wind after the initial projection | 0.1937 m/s | **0.0972 m/s** |
    | max meridional wind, iter 20 | 3.562967 m/s | 3.652312 m/s |
    | water drift, 20 iters | +0.0630 % | **+0.0630 %** |
    | drift by level, top 6 | i=44 +2.82e-03, i=45 +2.80e-03, … | **identical to 3 s.f.** |
    | CO₂ drift; min vs max | 0.0000 %; equal | 0.0000 %; equal |

    The anelastic residual falls 22 %, the Boussinesq one rises — it is no longer the
    enforced quantity — and half the spurious radial wind the old projection left in the
    initial field is gone. **Step 6 turns out to be unnecessary**: the `p_dyn_cap` source
    clamp, which the scope suspected of clipping the projection before it could act, binds
    in **0 of 3 791 399 fluid cells**. That is now printed every solve rather than assumed.

    And the water drift does not move — the third repair in a row to leave it identical.
    A tracer error indifferent to a velocity field that changed by a factor of two in `u`
    is not an advection error, which rules out the flow exactly as the algebra of item 15
    ruled out the tracer equation.

    **The drift is in the diagnostic's denominator.** `waterBudget()` reports water mass
    over air mass, and *both* come from `p_stat`, which `densities()` re-integrates
    hydrostatically every iteration on the local `R` and `cp` — which follow the
    composition and the surface temperature. Nothing separated the two. Measured against
    weights frozen at the reference time, with the column air mass those weights carry
    reported beside it:

    ```
    iter 10:  q_mean drift +0.0621 %   q against FROZEN weights +0.0011 %   column air mass -0.1258 %
    iter 20:  q_mean drift +0.0630 %   q against FROZEN weights +0.0025 %   column air mass -0.2128 %
    ```

    The water field moved by **+0.0025 % over 20 iterations**, twenty-five times less than
    the number this README has been quoting. What moved is the air: the column is losing
    about 0.01 % of its mass per iteration, steadily and without sign of stopping, so
    water-over-air rises. The per-level attribution at 113–149 km is the same artefact —
    those are the levels where `dp` changed most, not where water arrived.

    So **item 15's headline is wrong and is corrected here**: water was not being created
    at 0.11 % per iteration. The transport error is ~0.0001 % per iteration, and what the
    budget was measuring is the hydrostatic column being re-weighed. Two consequences:

    - The **column air mass is not conserved**, because nothing in this model makes it a
      conserved quantity. `p_stat.x[0]` is re-anchored every iteration to
      `r_air·R_mix·T_surf`, so the mass of a "250 bar atmosphere" follows the surface
      temperature — which the model evolves (1496.6 K against the prescribed 1500) — and
      the whole column is then re-integrated on the new anchor. It is a prescription
      defect, not a transport defect, and it is the open one.
    - The `c ≤ 1 − co2` ceiling is still deleting water (0.000164 kg/kg over 20 iterations),
      and that deletion is real. It is now the larger of the two mass errors.

    **A fourth repair, correct and not the cause.** `ATM_TRACER_DIFF_FLUX` (default 0) adds
    the missing term of the diffusive flux: for a mass fraction the conservative form is
    `∂(ρq)/∂t = ∇·(ρK∇q)`, i.e. `∂q/∂t = K∇²q + K∇lnρ̄·∇q`, and the second term was absent —
    the tracer diffusion was conserving `∫q dV` rather than `∫ρq dV`, in the same way the
    advective form does but *without* continuity to repair it. Unlike the advective
    flux-form correction of item 15 this one is proportional to `∇q`, so a uniform tracer
    stays uniform and CO₂ stays well mixed. Worth having and kept; measured effect at 20
    iterations is the fifth decimal of `residuum_atm` and nothing else. Note also that
    `diff_prec_re_inv` is **not** multiplied by `diffusion_ramp`, so moisture diffusion runs
    at full strength through a spin-up in which heat and momentum diffusion are ramped from
    zero. That asymmetry is undocumented and probably unintended.

    **`ATM_ANELASTIC` is left off by default.** It is the right continuity for this column
    and it measurably improves the projection, but 20 iterations is not evidence of
    stability, and the scope's own warning — every result in this README moves with it —
    stands. Flip it after a 400-iteration run, the way the metric fix was flipped.

## Scope: the anelastic continuity fix

**Status: steps 1–7 are implemented and measured (items 16 and 17). The premise below —
that the tracer mass error comes from the velocity field — is false; the error was in the
budget's denominator, not in the flow. The scope is kept because the anelastic projection
is right on its own terms and the reasoning is what the measurements were made against.**

The tracer mass error of item 15 comes from a velocity field that does not satisfy
continuity for the prescribed density. `PressureSolverAtm` projects onto `∇·u = 0`; a
column whose density spans five orders of magnitude needs `∇·(ρ̄u) = 0`. This is the scope
of that change, written before touching it because it is the dynamical core and every
result in this README moves with it.

**What the solver does now.** `PressureSolverAtm::run()` builds a provisional velocity
(`aux_u/v/w`), takes its divergence with single-power metric factors, solves a Poisson
equation for `p_dyn` with a variable-coefficient Laplacian, and corrects the velocity by
the pressure gradient. The file documents its own defect: the Laplacian coefficients use
`exp_2_rm`, `inv_rm2`, `inv_rm2sinthe2` while the divergence source and the gradient
correction are single-power, so `div` and `grad` are **not discretely adjoint** and the
projection does not exactly remove the divergence it measured. Three stabilisers
(`p_dyn_cap = 2.0`, the `p_dyn_ceiling`, the topography Dirichlet pins) were calibrated
against that inconsistency.

**The change, in order of dependency.**

1. **Resolve the metric inconsistency first.** Until `div` and `grad` are adjoint, no
   projection — Boussinesq or anelastic — removes the divergence. This is the file's own
   open TODO and it is independently testable: measure `∇·u` before and after the
   projection and require it to drop. *Do this step alone first and re-measure the water
   drift; it may account for much of it.*
2. **Introduce a base-state density `ρ̄(z)`** — one-dimensional, time-independent, the
   horizontal average of the prescribed profile. It must be a function of height only:
   a fully three-dimensional ρ makes the elliptic operator time-varying and costs
   solvability. Rebuild it when `densities()` rebuilds the profile.
3. **Anelastic divergence source.** `D = (1/ρ̄)∇·(ρ̄u*) = ∇·u* + u*_r · dln ρ̄/dr`. Because
   ρ̄ depends on r only, this is one extra term on the radial component.
4. **Anelastic Poisson operator.** `∇·(ρ̄∇φ) = ρ̄∇²φ + (dρ̄/dr)(∂φ/∂r)`, so the existing
   7-point stencil gains a first-derivative term in `r` proportional to `dln ρ̄/dr`. No new
   solver is needed — the change is to the stencil coefficients.
5. **Boundary conditions.** Impose `ρ̄u_r = 0` at the surface and the lid. `aux_u` is
   currently `c43/c13`-extrapolated at both ends, which does not impose zero normal mass
   flux.
6. **Re-derive the stabilisers.** `p_dyn_cap` and the ceiling were tuned against the
   inconsistent operator and will otherwise clip the anelastic projection before it acts.
7. **Check the buoyancy reference.** `BuoyancyForce = coeff_buoy·(t − t_ref_level[i])` is
   already base-state-referenced; confirm the reference is the same ρ̄.

**Verification — the point of the diagnostics already in place.**

| test | required | result (item 17, `ATM_ANELASTIC=1`, 20 iters) |
|---|---|---|
| `waterBudget()` drift, 20 iters | ≈ 0 | +0.0630 %, **unchanged** — and the wrong test: see item 17 |
| water drift against frozen weights | ≈ 0 | +0.0025 %, the number that was wanted all along |
| water-vapour ceiling hits | ≈ 0 | still deleting 0.000164 kg/kg / 20 iters |
| CO₂ max − min (uniformity) | still 0 | 0 exactly |
| `∇·(ρ̄u)/ρ̄` after projection | measured, and small | 3.359e-02 → **2.607e-02** |
| `∇·u` after projection | measured | 2.153e-02 → 2.670e-02 (no longer the enforced one) |
| source clamped at `p_dyn_cap` | not clipping the projection | **0 of 3 791 399 cells** |
| Ψ_max over 400 iters | bounded | not yet run at 400 |
| bit-identical at 1/4/8 threads | still passes | not re-checked |

The CO₂ uniformity test is the one that catches an over-correction: any scheme that makes a
well-mixed tracer develop structure is wrong, whatever it does for the mass budget. It is
what ruled out the advective flux-form correction, and it passes here because the anelastic
change is in the flow rather than in the tracer equations.

**Effort and risk.** The code is modest — a few dozen lines across `PressureSolverAtm` and
the divergence source. The verification is the work, and every number in this README
changes. Step 1 is separable, independently valuable, and should be measured before
steps 2–7 are started.

That estimate held: about forty lines, and the verification took four 20-iteration runs.
What it did not anticipate is that all of it would be measured against a diagnostic whose
denominator was moving. Three repairs — step 1, steps 2–7, and the diffusive flux term —
each left the water drift identical, and *that* is what finally identified the diagnostic
rather than any one of them. A repair that changes nothing measurable is evidence about
the measurement.

18. **The 400-iteration run, the anelastic default, and a threading defect the family had
    already solved twice (done).**

    Item 17 left `ATM_ANELASTIC` off with one condition: flip it after a 400-iteration run,
    the way the metric fix was flipped. Both runs were made at `29ca2f9` with identical
    configuration, 24 threads, the only difference the knob.

    | at iteration 400 | Boussinesq | anelastic |
    |---|---|---|
    | mean T | 800.631 K | 800.630 K |
    | mean KE | 34.1103 m²/s² | 34.1147 m²/s² |
    | KE drift per window | 4.358 % | 4.360 % |
    | Ψ_max | 1 251 136 | 1 259 044 (1e9 kg/s) |
    | OLR | 679.78 W/m² | 678.79 W/m² |
    | imbalance | −408.96 W/m² | −407.97 W/m² |
    | water drift, frozen weights | +0.0502 % | +0.0502 % |
    | column air mass drift | −0.3421 % | −0.3421 % |
    | deleted by the `c` ceiling | 0.002364 kg/kg | 0.002364 kg/kg |
    | `p_dyn_cap` clamp | — | 0 of 3 791 399 cells |

    **No instability, no clipping, nothing to argue with.** The two formulations track each
    other to 0.03 K in mean temperature, 0.002 % in kinetic energy and 0.6 % in Ψ_max over
    400 iterations, and the three water numbers are identical to four digits — which is
    item 17's conclusion restated from a different direction, since a budget that does not
    move when the velocity field's continuity constraint changes is not a transport budget.
    `ATM_ANELASTIC` is now **on by default**; the env var still forces it off for A/B.

    **What the run also showed is that nothing is converged, and why.** Ψ_max grows
    linearly in both formulations — 658k, 851k, 1053k, 1251k at iterations 100/200/300/400
    — with no sign of turning over, and after iteration 200 its maximum sits at the surface.
    The `[vbudget]` series says why in one line:

    ```
    iter  20   vbar=-1.133   pgf=-0.00000   cor=+0.00615
    iter 200   vbar=+0.012   pgf=-0.00005   cor=+0.00622
    iter 400   vbar=+1.259   pgf=-0.00011   cor=+0.00618
    ```

    The meridional wind is in **free acceleration**: `dv/dt` equals the Coriolis term to
    three figures, straight through zero, for 400 iterations. Advection, diffusion and the
    Rayleigh surface drag are each below 1e-5 — the drag is ~300× too weak to matter, and
    the only force that can balance a Coriolis torque is the pressure gradient, which is at
    **1.8 % of it** and growing linearly. Extrapolating that rate puts geostrophic
    adjustment of order 10⁴ iterations away, not 400. So this README's "of order 400
    iterations are needed" was wrong by a factor of ~25, and 400 iterations is a stability
    check, not a convergence check. **This is common to both formulations and therefore
    says nothing about the anelastic default** — but it is the largest open question about
    the dynamics.

    **And the determinism check that was never re-run had stopped passing.** At 20
    iterations:

    ```
    ATM_ANELASTIC=1, 1 thread    residuum_atm = 0.75321620
    ATM_ANELASTIC=1, 4 threads                = 0.75328193
    ATM_ANELASTIC=0, 1 thread                 = 0.75451284
    ATM_ANELASTIC=0, 4 threads, run A         = 0.75455454
    ATM_ANELASTIC=0, 4 threads, run B         = 0.75458447
    ```

    Both paths diverge, at the same place — the first `project_initial_velocity`, before any
    anelastic code runs — so the anelastic work is exonerated. And the last two lines are
    the sharp ones: **the same binary at the same thread count gives different answers run
    to run.**

    The site is this file's own Poisson loop:

    ```
    p_dyn[i][j][k] = (p_dyn[i±1][j][k]*num1 + p_dyn[i][j±1][k]*num2
                    + p_dyn[i][j][k±1]*num3 + ... - div_src) * inv_denom
    ```

    written **in place** under `#pragma omp parallel for collapse(2) schedule(dynamic, 4)`
    — over the very two indices the stencil reads across. Cell (i,j,k) was read by the
    thread owning (i+1,j) or (i,j+1) while its owner was writing it, and `schedule(dynamic)`
    meant which thread got which chunk varied with timing.

    **This defect has now been found three times in this family.** ATOM's shared
    `PressureSolver.h` records fixing it by red-black colouring; ATURAN `ffd0e0e` found it
    again in its own solver and cured it by serialising; and CLAUDE.md has carried
    "in-place Gauss–Seidel as a threading defect (ATURAN `ffd0e0e`)" in its list of *traps
    already solved elsewhere — check before re-deriving* the whole time. It was listed and
    not checked. A cross-reference is not a check.

    The fix here is **red-black colouring**, not ATURAN's serialisation: ATURAN could
    serialise because its `computePressure` is 0.003 s of a 5.3 s step, while here
    `project_initial_velocity` alone is 26.7 s at one thread. Each solve is two passes over
    a checkerboard of `(i+j+k)`; every cell of one colour has all six stencil neighbours in
    the other, so within a pass nothing is read while it is written. Red-black rather than
    Jacobi because `k` runs serially inside a thread, so the loop was reaching for
    lexicographic Gauss–Seidel — correct in serial, broken only by the (i,j) parallelism —
    and Jacobi would have cost the convergence rate. The colour is selected with a
    `continue` rather than by striding `k`, because the `k` loop carries a sliding window
    over the land mask that assumes consecutive `k`.

    **Every number measured before this commit moves in its last digits**, including the
    400-iteration table above: red-black is a different sweep order from lexicographic
    Gauss–Seidel, so it reaches the same solution by a different path. Unlike ATURAN's
    serialisation, this fix does *not* reproduce the old single-thread answer, and the
    honest statement is that the comparison table is a comparison of two runs that were each
    individually irreproducible.

    **A second race was hiding under the first**, and only became visible once the pressure
    solver stopped drowning it out. With red-black in place `residuum_atm` came back
    identical at 1, 4 and 8 threads — but 50 log lines still differed, all of them
    downstream of one statement in `UtilsAtm::findResiduumAtm`:

    ```cpp
    if(res > local_max.val) {
        local_max = {res, i, j, k};
        m.residuum_old = res;        // shared model member, written from every thread
    }
    ```

    Each thread wrote its own running maximum into a shared member with no synchronisation,
    and whichever finished last survived — 0.75315196 at one thread against 0.28779950 at
    four. It is read three times in the printout, including the test that decides whether
    the log says *"absolute error declining"* or *"absolute error is too high"*, so **a raced
    value drove the line a human reads to judge convergence** — and it was never the previous
    iteration's residuum that the name promises and that test needs. Now captured once at
    entry and written once at exit.

    The same block also chose its reported error *location* by thread arrival order: a plain
    `>` inside `omp critical` lets the first thread win an equal maximum, giving lon 18 at
    one thread against lon 20 at four, same latitude, same residual. Ties are not an accident
    here — invariant 1 makes the model hemispherically symmetric, so **equal maxima are the
    expected case**. Tie-broken on the smallest `(i,j,k)`.

    **What is fixed, and what is not.** Measured at 20 iterations after both repairs:

    | comparison | before | after |
    |---|---|---|
    | 4 threads, run A vs run B | differ (0.75455454 / 0.75458447) | **identical** |
    | 1 vs 4 threads | differ, 334 lines | differ, 54 lines |
    | 1 vs 8 threads | differ | differ |

    The **races are gone** — that is what the run-to-run comparison proves, and it was the
    serious defect. What remains is ~1 ulp and a different animal: `residuum_atm` comes out
    0.74743479 / 0.74743479 / 0.74743481 at 1/4/8 threads while the printed wind extrema are
    identical between 1 and 4. That is floating-point summation order, not a race — OpenMP
    combines partial sums in a thread-count-dependent order and `+` is not associative — and
    it re-enters the physics through a global mean, `t_skin` being the prime suspect since
    `densities()` rebuilds the whole column from it. **So the accurate claim is: run-to-run
    bit-identical at a fixed thread count, thread-count dependent at ~1e-8.** Ordered
    reductions are a separate item.

    **And the elliptic solve was under-converged, but that is not why nothing balances.**
    `run()` did exactly one Gauss–Seidel sweep per physics iteration, which moves information
    one cell — ATURAN's shared solver has carried the note *"one sweep per call is not an
    elliptic solve"* for its whole history, and ATHAD never had the knob. `ATM_PRESS_SWEEPS`
    adds it, default 1 and bit-identical. The scan, at iteration 20:

    | | 1 sweep | 10 | 50 |
    |---|---|---|---|
    | `pgf` | −0.00000 | −0.00003 | −0.00006 |
    | `cor` | 0.00615 | 0.00614 | 0.00614 |
    | pgf/cor | ~0 % | 0.5 % | **1.0 %** |
    | `vbar` | −1.13245 | −0.98020 | −0.83438 |
    | `∇·(ρ̄u)/ρ̄` rms | 2.653e-02 | 2.335e-02 | 2.319e-02 |
    | wall clock | 9:07 | 6:39 | 11:32 |

    More sweeps do build more pressure-gradient force, and the meridional wind accelerates
    measurably less — so the solve **was** under-converged. But it saturates: the residual
    improves 12 % from 1 to 10 sweeps and 1 % more from 10 to 50, while the initial
    projection goes 18 s → 73 s → 357 s, and at 50 sweeps `pgf` is **still 1 % of the
    Coriolis term**. A projection's pressure removes divergence instantaneously; the balanced
    field is a slow mode the flow must build over many steps. Fifty sweeps buys a factor of
    about two in how fast, not the factor of fifty that would close the gap. **The ~10⁴
    iteration estimate stands, and the solver is ruled out as its cause** — which is the
    result worth having, since it was the obvious suspect.

19. **The column anchored to a mass, and the prescribed profile given a switch (done, and
    instrumented).**

    `p_stat.x[0]` was re-anchored every iteration to `r_air·R_mix·T_surf`. That holds the
    surface *density* fixed and lets the surface *pressure* follow the surface temperature,
    so every time the prognostic `T_surf` drifted off its prescribed 1500 K the whole
    250 bar column was re-integrated from a different anchor and the atmosphere gained or
    lost mass. The surface pressure of an atmosphere is the weight of the air above it: if
    no mass enters or leaves it is a constant, and the surface density is what follows.
    `p_prev = p_0`.

    | column air mass drift | old anchor | new anchor |
    |---|---|---|
    | iteration 10 | −0.1258 % | −0.0011 % |
    | iteration 20 | −0.2128 % | −0.0011 % |

    The old one grows every 10 iterations; the new one takes one rounding-level step and
    holds. Water against frozen weights improves with it, from a monotone +0.0025 % to an
    oscillation about zero (−0.0005, +0.0008 %) — which is what item 17 predicted, since
    the water field was never the problem, the denominator was.

    The conservation references also move to **after** the first `densities()` call. Both
    are mass-weighted by `p_stat`, so capturing them before it recorded a state the model
    never integrates from. That was harmless while the anchor moved with the surface
    temperature — the offset was buried in the drift it caused — and became visible the
    moment the anchor stopped moving, as a fixed +0.62 % reported as drift.

    **`ATM_PROGNOSTIC_T` (default 0, bit-identical) stops `densities()` overwriting `t`
    with its adiabat**, so the temperature the dynamics and the radiation computed survives
    the iteration. At 20 iterations it did not blow up, no NaN, run completes — and the
    energy imbalance fell from −433 to −19 W/m². **That −19 W/m² is retracted in item 20**;
    it was the removal of a constraint reporting itself, not the physics. What was worth
    keeping is the profile it produced: the lower column cooling at 5.97 K/km against a
    4.19 K/km dry adiabat — **super-adiabatic, convectively unstable** — and a grid-scale
    sawtooth in the top 50 km, 366.8 / 333.0 / 252.7 / 263.6 / 318.2 / 231.0 K on adjacent
    levels.

    Both follow from one fact: **this adiabat is the model's only convective adjustment.**
    ATHAD had no `ConvectiveAdjustment` at all, and `MoistConvection`'s triggers are
    absolute Earth pressures that never fire at 250 bar, so switching the prescription off
    does not give a radiative–convective profile — it gives a radiative one, which in an
    optically thick atmosphere is super-adiabatic by construction. Hence item 20 before the
    flip.

20. **A convective adjustment, ported and generalised — and yesterday's −19 W/m² retracted
    (done).**

    `ConvectiveAdjustment.h`, dry Manabe–Strickler, from the family's shared
    `planet/ConvectiveAdjustment.h` (ATURAN/ATJUP/ATSAT). **Not copied verbatim, and the
    two reasons are the usual ones:**

    1. The shared file assumes a **uniform grid** — `dz_m = L_atm*1e3/(im-1)` is a layer
       thickness only when the layers are equal. ATHAD's radial coordinate is exponentially
       stretched, `dz` running 0.81 km at the surface to ~15 km at the lid, a factor of 19,
       so one critical drop would be 19× too strict at the bottom. Same defect class as
       `init_tropopause_layers`' `round(h / L_atm)`. The critical drop is now per layer and
       the segment adiabat a cumulative sum rather than the linear `dT_ad·(q−a)`.
    2. It assumes a **constant `cp_mix`**. ATHAD's cp varies twofold across 300–1500 K and
       follows the composition, so cp comes from `AtmMixture::cp_of()` locally — which is
       what the shared file's own comment warns to do.

    Enthalpy drift 1.6e-15 (prescribed) and 2.3e-14 (prognostic): the generalisation
    conserves what the algorithm promises. With `ATM_PROGNOSTIC_T=1` the sawtooth of item 19
    is gone — 609.1 / 548.9 / 484.3 / 417.9 / 347.2 / 272.2 K, monotone, and by construction
    nowhere super-adiabatic. The adjustment does real work there: max ΔT 291.6 K in a single
    column, 3.9 M layers mixed, against 3.18 K with the profile prescribed, where
    `densities()` overwrites `t` on the next line anyway. It is wired in immediately before
    `densities()`, so with the prescription on it costs one pass and reports the residual
    instability the dynamics generate per iteration — which is itself a measurement of how
    much work the prescription has been doing.

    **And it retracts item 19's headline.** The −19 W/m² came from the sawtooth: a
    convectively unstable column radiates from the wrong levels, and the apparent balance
    was the instability flattering the OLR. With the column forced back onto its adiabat:

    | at 20 iterations | OLR | imbalance | σT_lid⁴ |
    |---|---|---|---|
    | prescribed | 704.24 W/m² | −433.42 W/m² | — |
    | prognostic + adjustment | 634.46 W/m² | −363.65 W/m² | 169.29 W/m² |

    A 16 % improvement, not a factor of twenty. The imbalance is real, item 11's claim that
    the opacity is too low to hold the surface stands, and **nothing in the κ question is
    repaired by making the profile prognostic**. The lesson is the one this README keeps
    relearning: a number that improves dramatically when a constraint is removed is usually
    reporting the removal, not the physics. (Both rows predate item 22 and are superseded in
    magnitude; the comparison between them is not.)

21. **The water ceiling localised, and the mechanism its own comment asserted refuted
    (measured — and the repair deliberately not made).**

    The `c ≤ 1 − co2` ceiling deletes water: 0.000165 kg/kg per 20 iterations, 0.002364 per
    400. The plan was to make it conservative. Measuring first was right, and stopped it.

    The comment above the ceiling asserted a mechanism — water is pumped down out of the one
    condensing level by sedimentation and evaporates into the superheated band below with no
    return path, so `c` there grows until something stops it. A new per-level report tests
    it. At 20 iterations with `moist_phys_start_iter = 300`, so `SaturationAdjustment`, the
    ice schemes and **all sedimentation are off and never run**, the ceiling still bites
    16 606 cells:

    ```
    by level:  i=51 (185 km)  5 054      i=52 (195 km)  11 552
    ```

    There is no sedimentation running, so sedimentation is not the cause. And the count is
    **identical at iterations 10 and 20**, with the same per-level split, so it is not an
    accumulation either. Both halves of the asserted mechanism fail.

    What it looks like instead: a standing excess left at i = 51–52 by the one saturation
    adjustment that *does* run, at initialisation, re-clipped to the same value every
    iteration by something that either nudges `c` back above the ceiling or nudges the
    ceiling below `c` — `1 − co2` moves when the transported CO₂ does, and CO₂ has been
    prognostic since item 12. Which of the two has not been established.

    **So the ceiling is left deleting.** A conservative ceiling that redistributes instead
    of deleting would convert a bounded 0.35 %-per-400-iteration loss into an unbounded
    pile-up somewhere else, and until the mechanism is known the deletion is the only thing
    holding that band at a physically possible composition. The diagnostic is the deliverable
    here; the repair is not, and saying so is the point.

22. **The composition laid down before the column that reads it — and the OLR halves
    (done).**

    `initTemperatureData` reads `c` and `co2` at every level — `R_of(q_v, q_c, R_bg)` for the
    hydrostatic step, `cp_of(q_v, q_c, T, M_bg)` for the adiabat — but `initWaterWapour` and
    `co2Atmosphere` ran **after** it, so both fields were still zero. With `q_v = q_c = 0`
    the gas constant collapses to the background, R_bg = 317.3 instead of R_mix = 387.9, an
    18 % short scale height, and cp is the background-only value so the adiabat is too steep.
    **The whole initial column was built for an atmosphere this model does not have.** Both
    fills are unconditional uniform assignments with no dependence on `t` or `p_stat`, so
    they simply move ahead of it; `initCloudIce` stays behind, since it genuinely needs the
    profile.

    This README said the defect was survivable because `densities()` rebuilds the column
    straight afterwards, and that the only escapee was the lid snapshot — refreshed in
    item 10 rather than repaired at source. **That accounting was incomplete.**
    `initCloudIce` *and* the initial `SaturationAdjustment` both run between
    `initTemperatureData` and `densities()`, so they laid the initial cloud and ice fields on
    the wrong column. At 20 iterations:

    | | before | after |
    |---|---|---|
    | max cloud water | 196 km | 243 km |
    | max cloud ice | 186 km, 36.6 g/kg | 243 km, 29.2 g/kg |
    | planetary albedo | 0.4998 | 0.4981 |
    | OLR | 704.24 W/m² | **337.39 W/m²** |
    | imbalance | −433.42 W/m² | **−66.13 W/m²** |

    The cloud deck moves up about 50 km and the OLR halves. The consistency check that says
    which state is right is the model's own profile report: *"condensation possible from
    i = 55 (230.2 km) upward"*, and the old deck sat at 186–196 km — **in a band where this
    atmosphere cannot condense at all.** The new one sits at 242.8 km, inside it. The deck
    was in an impossible place because it was built on a column that did not exist.

    **This also means the −409 to −433 W/m² imbalance quoted since item 11 was carrying a
    second artefact** on top of the one item 7 looked for and did not find. It does *not*
    retract item 11's conclusion — 337 W/m² against 271 absorbed is still an atmosphere
    radiating away more than it takes in — but the margin is a third of what was claimed, and
    every κ argument resting on the old figure needs redoing against this one. **Every
    radiation number in item 18's 400-iteration table predates this fix and is superseded**;
    the anelastic-versus-Boussinesq comparison in it is not, since both columns carried the
    same artefact.

23. **The prognostic column re-measured on the corrected baseline — it is now the worse of
    the two — and item 21's open question answered by item 22 (done).**

    Two 20-iteration runs at this commit, identical configuration (`nm = 20`, 8 threads),
    the only difference `ATM_PROGNOSTIC_T`. The prescribed column reproduces item 22's
    numbers to the last digit, so the pair is comparable:

    | at iteration 20 | prescribed (default) | prognostic + adjustment |
    |---|---|---|
    | OLR | 337.39 W/m² | 497.61 W/m² |
    | imbalance | −66.13 W/m² | −226.35 W/m² |
    | planetary albedo | 0.4981 | 0.4981 |
    | lid T / σT_lid⁴ | 262.95 K / 271.10 W/m² | 259.58 K / 257.47 W/m² |
    | deleted by the `c` ceiling | 0.000000 kg/kg, 0 cells | 0.000000 kg/kg, 0 cells |
    | column air mass drift | −0.0000 % | −0.0024 % |

    **The sign of item 20's comparison has reversed.** Before the composition-ordering fix,
    making the profile prognostic *improved* the imbalance by 16 % (704 → 634 W/m²) and that
    was the argument for flipping the default. With the cloud deck in its right place the
    prescribed column falls to 337 W/m² and the prognostic one only to 498, so **the
    prescription is now better by a factor of 3.4** — and the whole difference is in the top
    five levels:

    | height | prescribed | prognostic |
    |---|---|---|
    | 242.8 km | 263.0 K | 309.0 K |
    | 256.0 km | 263.0 K | 352.3 K |
    | 270.0 km | 263.0 K | 337.5 K |
    | 284.6 km | 263.0 K | 246.0 K |
    | 300.0 km | 263.0 K | 257.1 K |

    Everything above ~243 km is pinned at the isothermal skin `t_skin` = 263.0 K when the
    profile is prescribed. **That skin was holding the OLR down**, and it is an assumption
    rather than a result: `t_skin` is a fixed point against the model's own albedo (item 10),
    not a temperature the radiation computed at those levels. So the −66 W/m² of item 22 is
    the better-looking number and the more constrained one; the prognostic −226 W/m² is what
    this radiation scheme actually produces when nothing pins its top.

    **`ConvectiveAdjustment` is right not to touch that top.** 243 → 256 km is a temperature
    *increase*: an inversion is convectively stable, and a dry adjustment leaves it alone by
    construction. It is not idle elsewhere — 100 % of columns, 6.47 M layers, max ΔT 127.7 K,
    every iteration, and **steady from iteration 2 to 20**, which says the radiative–dynamical
    tendency regenerates the instability as fast as the adjustment removes it. What the
    prognostic top needs is not more mixing, and item 20's structural improvement (the
    sawtooth is gone; the lower column is monotone) survives this result intact.

    **The water ceiling now never fires, and item 21's question is answered — by item 22, and
    by none of the three mechanisms that had been proposed for it.** Zero cells and
    0.000000 kg/kg deleted in both runs, against 16 606 cells and 0.000165 kg/kg per 20
    iterations before the ordering fix. The tell was in the refuted report all along: the
    clipping sat at **i = 51 (185 km) and i = 52 (195 km)**, and the misplaced cloud deck sat
    at **186–196 km**. It was the same object. Neither the sedimentation mechanism the
    comment asserted, nor the `c`-rises/`1−co2`-falls pair item 21 proposed in its place, was
    the cause; the deck had been built on a column the model does not have, in a band where
    this atmosphere cannot condense, and the ceiling was reporting that. **Two mechanisms
    asserted, both wrong, and the measurement that settled it was a location.** The ceiling
    stays as a diagnostic — the untested case is a run past `moist_phys_start_iter = 300`,
    where sedimentation finally runs — and it is not to be made conservative while it has
    nothing to delete.

    **One thing this pair raises that is not yet explained**: in the prescribed run the
    equatorial column radiates **271.1 W/m²** while the global mean is 337.4, so the hottest
    surface on the planet sits under the *least* emitting column. With a surface suppression
    of ×1049 and a transmitted fraction of 0.000 the OLR is set entirely by the cloud and
    skin structure aloft, not by the surface beneath it — which is expected in an optically
    thick atmosphere, but the latitudinal sign of it has not been checked and belongs in the
    κ scan.

24. **41 levels against 61: the conclusion survives, the radiation detail does not (done).**

    The question was whether the shell could go back to `ATOM_Precipitation`'s `im = 41`.
    **It is worth separating two things that look like one.** The 300 km shell is not a
    resolution choice: the vertical extent of an atmosphere is `H = R·T/g`, which is 59.3 km
    here against Earth's 8.4 — a factor of 7.1, of which the composition supplies 1.35 (M is
    21.43 g/mol against 28.96, so R is *larger*) and the surface temperature supplies 5.2.
    At 15.9 km, where ATOM's lid sits, this column is still at **189.8 bar**, with 76 % of
    the atmosphere's mass above the boundary. Nor can the lapse rate compensate: bringing
    1500 K down to a ~263 K radiating temperature within 16 km needs 77 K/km, hence
    `cp = g/Γ = 127 J/(kg·K)`, against this mixture's 2040 and pure xenon's 158. Earth's own
    6.5 K/km applied to a 1500 K surface still needs 190 km. The measured 4.2 K/km reaches
    `t_skin` at ~295 km, which is where the 300 km shell came from.

    So the shell stays and `im` alone moves. Two 20-iteration runs, same configuration, same
    `dt_visc = 4e-5`, machine idle for each:

    | at iteration 20 | `im = 41` | `im = 61` |
    |---|---|---|
    | wall clock | **2:01.2** | 2:54.6 |
    | RK4 per call | 3.79 s | 5.75 s |
    | peak RSS | 2.58 GB | 3.82 GB |
    | surface / top layer | 1.22 / 22.8 km | 0.81 / 15.4 km |
    | OLR | 349.73 W/m² | 337.39 W/m² |
    | imbalance | −78.67 W/m² | −66.13 W/m² |
    | **lid emissivity** | **0.0068** | **0.0000** |
    | σT_lid⁴ / lid T | 270.90 W/m² / 262.91 K | 271.10 W/m² / 262.95 K |
    | top of domain | 0.00032 bar | 0.00028 bar |
    | planetary albedo | 0.4988 | 0.4981 |
    | Ψ_max | 492 263 @ 44.9 km | 543 769 @ 48.0 km |
    | condensation begins | i = 37, 236.4 km | i = 55, 230.2 km |
    | equatorial column OLR | **403.3 W/m²** | 271.1 W/m² |
    | `c` ceiling | 0 cells | 0 cells |
    | column air mass drift | −0.0037 % | −0.0000 % |
    | `residuum_atm` | 0.64277089 | 0.74743459 |

    **The cost saving is real and unremarkable: 1.44× in wall clock, 0.68× in memory**, both
    close to the 61/41 cell ratio. No superlinear gain, and none expected.

    **The physical conclusion is robust to the vertical resolution.** The OLR moves 3.7 % and
    the imbalance 19 %, in the same direction, and "this atmosphere radiates away more than
    it absorbs" holds at both resolutions. 12 W/m² is far inside the factor-of-two κ
    uncertainty, so nothing in item 22's argument depends on the grid.

    **Two things degrade, and one of them is the thing item 11 bought.** The lid emissivity
    goes 0.0000 → 0.0068: with 22.8 km top layers the top cell is ~2.2 local scale heights
    (H = 10.4 km at 263 K) and is no longer transparent. It has not cost the decoupling — OLR
    349.7 against σT_lid⁴ 270.9 is still a genuine column integral, not the boundary
    temperature read back out — but that is the margin a deeper shell or a coarser grid eats
    first, and it should be checked, not assumed, in any configuration that changes either.

    The larger effect is on **where the cloud deck lands**. Condensation begins at 236.4 km
    against 230.2, and the equatorial column's OLR comes out **403.3 W/m² against 271.1** —
    49 % in a single column, on exactly the sensitivity item 22 exposed. Ψ_max differs 9.5 %
    with its maximum one layer lower. **A grid that locates the deck to ±11 km is not a grid
    to quote an OLR from**, even though it gives the same answer about the sign of the
    imbalance.

    This is also the **second grid-convergence data point** and the first clean one. The
    519-versus-581 W/m² pair still quoted below is at two different shell depths with `im`
    fixed, so it confounds depth with resolution, and it predates item 22 besides. This pair
    holds the shell fixed and moves only `im`.

    Incidentally, the `im = 61` run above was made fresh and reproduced item 22's
    337.39 / −66.13 **exactly**, which is item 18's "run-to-run bit-identical at a fixed
    thread count" holding on a different day and a different working tree.

    **The CFL headroom is real, and it is the actual prize.** `dt_visc` was cut from the
    inherited `1e-4` to `4e-5` when `im` went 41 → 61, because the explicit diffusion limit
    goes as the *physical* surface spacing squared and that spacing went 1.22 → 0.81 km, a
    factor of 2.25. Going back to 41 levels should hand it back. Two 200-iteration runs at
    `im = 41`, identical but for the step, **both completed, exit 0, no NaN, and both cleared
    iteration 174** — the iteration at which `5e-4` is documented to have blown up the
    near-surface cells:

    | | `dt = 1e-4` | `dt = 4e-5` |
    |---|---|---|
    | iteration 20 | OLR 333.09, imb −62.03 | OLR 349.73, imb −78.67 |
    | iteration 100 | OLR 279.18, imb −8.07 | OLR 294.35, imb −23.25 |
    | iteration 200 | OLR 272.37, imb −1.26 | OLR 280.75, imb −9.64 |
    | Ψ_max at 200 | 1 307 027 @ 43° | 710 833 @ 44° |

    Those columns are **not** at the same physical time — the point of a bigger step is that
    they cannot be. Compared where the physical time *does* match, `dt=1e-4` at iteration 80
    against `dt=4e-5` at iteration 200, both at 8.0e-3 non-dimensional:

    | at t = 8.0e-3 | `dt = 1e-4`, iter 80 | `dt = 4e-5`, iter 200 |
    |---|---|---|
    | OLR | 283.55 W/m² | 280.75 W/m² |
    | imbalance | −12.45 W/m² | −9.64 W/m² |

    **1.0 % in the OLR for 2.5× fewer iterations.** With `im = 41`'s 1.44× per iteration on
    top, that is ~3.6× less wall clock to reach a given physical state, which is the number
    that matters against item 18's ~10⁴-iteration estimate for geostrophic adjustment.

    Two limits on that result. It validates `1e-4` **at `im = 41` only** — the CFL argument
    says the same step is ~2.2× over the limit at 61 levels, so it must not be carried back to
    the fine grid. And 200 iterations with `moist_phys_start_iter = 300` is a dry test: the
    stiff microphysics has never seen this step.

25. **The energy imbalance is a spin-up transient, and the converged OLR may be the skin
    temperature read back out (measured; the discriminating test is named, not yet run).**

    The 200-iteration runs of item 24 were meant to test a time step. What they showed is
    that **the imbalance this README has been quoting since item 11 decays monotonically
    toward zero**:

    ```
    dt = 4e-5:  -78.67  -53.66  -38.83  -29.48  -23.25  -18.88  -15.66  -13.19  -11.23  -9.64
    dt = 1e-4:  -62.03  -33.96  -20.01  -12.45   -8.07   -5.39   -3.68   -2.55   -1.78  -1.26
                (iterations 20, 40, 60, 80, 100, 120, 140, 160, 180, 200)
    ```

    No flattening, no sign change, no bounce — a clean decay in both, at a rate that scales
    with the step. **Every imbalance figure in items 11 through 24 was measured at iteration
    20 and is therefore a snapshot of a decaying transient**, including item 22's −66.13 and
    item 23's comparison of the prescribed against the prognostic profile. The claim built on
    them — *this atmosphere radiates away more than it takes in, so `kappa_H2O` is too low* —
    does not survive at 200 iterations, where the imbalance is −1.26 W/m².

    **But the way it closes is exactly the failure mode item 10 identified.** Over the same
    200 iterations `t_skin` moves 262.91 → 262.96 K and σT_lid⁴ sits at 271.1 W/m²,
    essentially constant; what moves is the OLR, 349.7 → 280.8, *downward onto them*. The
    isothermal top of the prescribed profile is at `t_skin`, and `t_skin` is a fixed point
    solved from σT_skin⁴ = absorbed SW + geothermal. So once the transient cloud and
    temperature structure decays and the effective radiating level migrates up into that
    isothermal skin, **OLR = absorbed is guaranteed by construction, not discovered**. A
    converged column *must* satisfy it, which is why the coincidence is not by itself proof
    of anything — and why it cannot be used as evidence that the opacities are right either.

    **The discriminating test is a κ scan run to 200 iterations, not to 20.** If the converged
    OLR lands on σT_skin⁴ for every `kappa_H2O`, the outgoing flux is still an input on long
    runs and the grey scheme is only reporting its boundary condition — item 10's defect
    surviving the item 11 rewrite, in a form that only shows up after the transient dies. If
    instead the converged OLR moves with κ, the column integral is doing real work and the
    balance is a result. **Nothing about the opacity should be claimed in either direction
    until that run exists.** The same test settles item 24's other loose end, since it is the
    converged OLR, not the 20-iteration one, that has to be grid-converged.

    Two smaller things the long runs also settle. The `c ≤ 1 − co2` ceiling deleted
    **0.000000 kg/kg over 200 iterations**, so item 23's result is not an artefact of stopping
    at 20; and the column air mass drift is −0.0036 % over 200 against −0.0011 % over 20, so
    item 19's anchor holds over ten times the length. Both still with the moist physics off.

26. **The Hadley cell was not decaying — the diagnostic was losing it, and the initial state
    was never balanced (done).**

    The question was why the cell structure visible at iteration 0 fades over the next
    hundred. Three things turned out to be involved, and only the third is the cause.

    **The streamfunction was a volume flux wearing a kg/s label.**
    `write_meridional_streamfunction` carried `const double rho = r_air` — one constant, the
    *surface* density, applied from the ground to the lid. On Earth's 16 km shell ρ spans
    about 6× and the error is a stretch factor; here it spans four orders of magnitude, so a
    cell at the lid was weighted ~40 000× too heavily and the thick upper layers — 23 km each
    at `im = 41` — dominated an integral that belongs to the bottom few kilometres. ρ now
    comes from `r_humid` and is averaged zonally *together* with v, keeping the ⟨ρ'v'⟩
    correlation a warm rising branch carries. Both fields from the same 20-iteration run:

    | | Ψ_max | latitude | height |
    |---|---|---|---|
    | old (ρ const) | 492 341 | +45° | 44.9 km |
    | new (ρ(z) inside) | 155 762 | **+15°** | 0 km |

    **It does not merely rescale the plot.** The peak moves from mid-latitudes into the
    tropics, and at 15° the corrected Ψ *reverses sign* at ~60 km — −1 849 at 66 km, +7 823 at
    55 km, +155 762 at the surface — while the old curve is single-signed at every height with
    its largest magnitude at ~95 km. A sign reversal with height at fixed latitude is a closed
    overturning cell. Share of |dΨ| below 55 km: **old 30.6 %, new 79.1 %.**

    At iteration 100 of a moist run the same comparison gives old 879 470 at +44° against new
    134 484 at +15°, and the vertical profiles say what the picture is:

    ```
    lat 15N:  0km:134484  36km:28904  55km:5810  79km:-8261  113km:-10632   <- reverses
    lat 30N:  0km:107785  36km:77154  55km:61287 79km:42319  113km:22401    <- never
    lat 45N:  0km:127060  36km:99210  55km:80874 79km:57298  113km:31015    <- never
    lat 70N:  0km:47754   36km:25743  55km:18280 79km:11062  113km:5044     <- never
    ```

    **Exactly one closed cell exists**, tropical, below 60 km. Everything poleward of ~22° is
    single-signed from surface to lid: mass moving poleward through the entire depth with no
    return branch, which is not a cell at all. The 45° surface maximum is the peak of that
    drift, not a Ferrel cell — and it is what the old diagnostic was reporting as the Hadley
    circulation. **Every Ψ_max in items 18, 24 and 25 is a volume-flux number**: the growth
    those series show is real and so is the extremum migrating to the boundary, but the
    magnitudes and the latitude and height they report are not.

    **The radial Shapiro filter is acquitted, and this is item 18's lesson in reverse.**
    ATOM_Precipitation (`1e59daa`) and ASTIM (`cf43bfd`) both measured the radial Shapiro pass
    on u,v,w as a dominant momentum sink, and this file's inherited comment says it *"erodes
    the Hadley/Ferrel cells' two-branch v structure"* and carries *"~85 % of Mdot_net at
    30–60°"* — ATURAN `7288cda` measured the same 85 % from the temperature side. A/B here
    with `ATM_RADIAL_SHAPIRO_STRENGTH=0`:

    | iter | Ψ_max, filter on | filter off |
    |---|---|---|
    | 2 | 157 220 | 157 263 |
    | 14 | 154 463 | 154 708 |

    **0.16 %, and the mid-latitude cell dies at the same iteration either way.** The
    annotation is true where it was written — Earth, sharp jets on a 16 km shell — and false
    on 41 levels spread over 300 km, where the vertical structure is smooth enough that a
    4th-order filter in `i` barely touches it. Item 18 recorded the family naming a defect and
    nobody checking; this is the family naming a suspect and the check acquitting it. **The
    conclusion is grid-dependent and does not travel; the method does.**

    **What kills the cell is that nothing was ever balanced.** At 45° the surface
    streamfunction crosses zero at iteration ~5 and then grows linearly at ~1 335 per
    iteration — constant dv/dt — and the budget says why in one line: `cor = 0.01532`,
    `pgf = 0.00000`, advection, diffusion and drag all `0.00000`. `VelocityInitializer`
    imposes an analytic profile onto a field with no pressure structure to support it, so the
    Coriolis torque on the imposed wind is unopposed from iteration 0 and the two-branch
    structure is buried within five iterations. Not damped — **overwhelmed**.

    **`initBalancedState`, ported from ASTIM `cf43bfd` and re-derived against this model's own
    equation.** Setting u = v = 0 and ∂/∂φ = 0 and requiring `rhs_v = 0`:

    ```
    dp_dyn/dthe = w^2*cotanthe + 2*force_nd*costhe*w*metricRadius(rm)
    ```

    It writes `p_dyn` because **`p_stat` appears nowhere in the momentum equations** (only in
    the Held-Suarez σ) — the entire meridional pressure-gradient force in this model is
    `p_dyn`, and a 50 K equator-to-pole contrast over a 250 bar column exerts none of it
    directly. It runs *after* `project_initial_velocity`, which zeroes `p_dyn`; `run()` then
    relaxes `p_dyn` in place rather than replacing it, and the balanced field is measured to
    persist exactly — `pgf` identical at iterations 1, 2 and 3. This is also why item 18 found
    `pgf` stuck at 1 % of Coriolis after 400 iterations: a balanced pressure is nearly
    divergence-free, so it sits in the null space of what the projection solves each step.

    **The radius is `metricRadius(rm)`, not `rad.z[i]`, and getting it wrong cost a factor of
    20.** `inv_rm` carries the *planetary* radius (~21 in shell units) while `rad.z[i]` ~ 1,
    and the Coriolis term carries no radius at all, so for `−dpdthe*inv_rm` to cancel it,
    `dp/dθ` must carry `metricRadius`. `RungeKutta_Atm_Turb.cpp` warns about precisely this
    two lines from where the metric is built — *"a solver that disagrees with the RHS about
    the metric is the defect class this repo keeps finding"* — and the first version of this
    port reproduced it:

    | | pgf / cor | dv_dyn at iter 1 |
    |---|---|---|
    | unbalanced | 0.000 | 0.01532 |
    | balanced, `rad.z[i]` | 0.051 | 0.01532 |
    | balanced, `metricRadius(rm)` | **0.997** | **−0.00009** |

    **And it collides with an Earth-tuned clamp.** The balance peaks at **455** while
    `p_dyn_ceiling` is 10 during spin-up and 3 after iteration 300 — a backstop sized against
    *"the accumulated steep-orography value (~7.7)"* in a model with no orography. **94.7 % of
    (i,j) columns exceed it**, so `ATM_BALANCED_INIT=1` alone yields a 95 %-clipped field, and
    unlike a uniform scaling, clipping distorts the *shape*. `ATM_P_DYN_CEILING` overrides it;
    unset is bit-identical. `initBalancedState` now reports the fraction the clamp will
    truncate, because a balance the model then clips is not a balance.

    **Measured, 20 iterations, no NaN:**

    | iter | Ψ_max | return branch @45° | Ψ₀ @45° |
    |---|---|---|---|
    | 2 | 157 443 | 28 319 | −4 858 |
    | 8 | 157 373 | 28 265 | −4 833 |
    | 20 | 157 403 | **28 195** | −4 813 |

    The return branch is still there at iteration 20, having vanished by iteration 6 in every
    run without this. Ψ₀ at 45° is frozen instead of crossing zero and accelerating. And Ψ_max
    dips 0.04 % to iteration 8 and then **recovers** — it oscillates about the prescribed value
    instead of decaying, which is how ASTIM described its cured jet. Both knobs default off
    and the build is bit-identical: the balance-off run reproduces the pre-port run in every
    traced column.

    **`im = 41` becomes the default** on item 24's measurement. The shell depth is not a
    resolution choice — the vertical extent is `H = R·T/g`, 59.3 km here against Earth's 8.4 —
    but the level count is, and 41 buys 1.44× wall clock and 0.68× memory for 3.7 % on the OLR
    with the sign of the imbalance unchanged. What it costs is at the top, so an absolute OLR
    still wants re-measuring at 61.

27. **The balanced initial state, measured over 200 moist iterations and made the default —
    and the radiation does not notice (done).**

    Item 26 cured the spin-down over 20 dry iterations. That is not enough to move a default,
    so the A/B was repeated at κ_H2O = 0.010 with the moist physics running from iteration 0,
    `im = 41`, `dt_visc = 1e-4`, 200 iterations, identical but for the balance. **It survives
    the microphysics**: no NaN, 34 m 50 s, despite `p_dyn` running at ~455 against an
    inherited ceiling of 10.

    | | unbalanced (control) | balanced |
    |---|---|---|
    | Ψ_max at iter 20 | 153 301 @ 15° | 157 584 @ 15° |
    | Ψ_max at iter 100 | 134 483 @ 15° | 157 054 @ 15° |
    | Ψ_max at iter 200 | **255 670 @ 45°** | 157 555 @ 15° |
    | return branch @45°, iter 200 | **none** | 30 770 |
    | Ψ₀ @45°, iter 200 | **+255 670** | −504 |

    **The control loses the circulation and grows an artefact in its place.** Its tropical cell
    decays 27 %, and by iteration 120 the global Ψ_max has moved to 45° — to the one-signed
    drift, which then doubles to become twice the strength of the cell it replaced. The
    balanced run instead **dips 0.34 % by iteration 100 and returns to its starting value by
    200**: an oscillation about the prescribed state, not a decay, which is how ASTIM described
    its cured jet. The mid-latitude return branch, which never exists in the control, gains
    9.4 % over the run.

    **And the radiation is completely indifferent to all of it:**

    | at iteration | OLR unbalanced | OLR balanced | imbalance unbal. | bal. |
    |---|---|---|---|---|
    | 20 | 329.23 | 329.22 | −58.16 | −58.15 |
    | 100 | 283.97 | 284.27 | −12.89 | −13.20 |
    | 200 | 273.34 | 273.43 | −2.26 | −2.36 |

    σT_lid⁴ = 271.31 in both. **Two circulations differing by 500× in their column-integrated
    drift and 27 % in cell strength give the same OLR to 0.03 %.** This is item 25's suspicion
    confirmed from a direction it did not use: the κ scan showed the outgoing flux does not
    move with the opacity, and this shows it does not move with the dynamics either. It
    relaxes onto `t_skin` regardless of what the atmosphere underneath is doing.

    **What this does not fix.** The *tropical* return branch erodes in both runs — −11 969 →
    −8 689 balanced against −11 787 → −9 670 unbalanced, so marginally **worse** with the
    balance. The 15° counter-flow aloft is weakening for some other reason, and that is now
    the open question about the cell rather than the drift that was masking it.

    **Defaults flipped, on the strength of the above.** `ATM_BALANCED_INIT` is **on** (=0
    restores the old path), and `p_dyn_ceiling` goes from the inherited phase-dependent
    10.0 / 3.0 to a flat **2000**. Both inherited values are keyed to terrain this model does
    not have — 10.0 against "the accumulated steep-orography value (~7.7)", 3.0 against a
    Tian Shan/Pamir pegging past moist onset — and at 455 the balance was 45× over the first
    and 94.7 % of columns were clipped. **The phase step is removed deliberately**: at
    iteration 301 it would have dropped to 3.0 and clipped the balance away mid-run.

    Two caveats, stated because flipping a stabiliser deserves them. **2000 is untested beyond
    200 iterations** — the Earth values existed because unbounded `p_dyn` blew this solver up,
    that failure was orography-driven and should not exist over a featureless surface, but
    "should not" is not a measurement. And nothing is bit-identical any more; every run from
    here differs from every run before it.

    A smaller thing caught in the act: the ceiling briefly lived in two places, and the
    balance diagnostic went on warning about a clip that no longer happened. It is now one
    accessor, `cAtmosphereModel::pDynCeiling()`, read by both the solver that enforces it and
    the diagnostic that reports against it — the same defect class as a solver disagreeing
    with the RHS about the metric, which is what cost item 26 a factor of 20.

28. **The balance had only one component, and the other one drove the vertical wind to
    11 m/s (fixed; the 200-iteration A/B is running).**

    **Found by looking at a picture and not believing it.** The question asked was why the
    ParaView glyphs in the `u-v-Cell` field show no circulation at all when the `u` field
    plainly has one. That turned out to be a plotting defect (below), but checking it meant
    reading `u` out of the vtk files directly — and that is where the real thing was.

    **The u field in the unbalanced control is a correct Hadley cell.** At iteration 20,
    ascent over the equator out to ~10°, descent at 20–30°, surface `v` at 15° flowing
    *toward* the equator, and symmetric to the last digit (u(15°N) = u(15°S) = −4.395e-4,
    v(15°N) = −v(15°S) = 2.898), so invariant 1 holds. Between iterations 0 and 20 only the
    bottom ~5 km at the equator changes sign; the deep pattern keeps its sign and loses 8 %
    over 200 iterations.

    **The balanced run — the default since item 27 — does not.**

    | iteration | unbal max \|u\| | unbal u at EQ, 95 km | bal max \|u\| | bal u at EQ, 95 km |
    |---|---|---|---|---|
    | 0 | 1.14e-1 | +6.08e-2 | 1.14e-1 | +6.08e-2 |
    | 20 | 1.06e-1 | +6.00e-2 | 2.53e+0 | −7.91e-1 |
    | 100 | 8.63e-2 | +5.50e-2 | 8.08e+0 | −4.08e+0 |
    | 200 | 9.47e-2 | +4.45e-2 | **1.12e+1** | **−6.35e+0** |

    The vertical velocity grows ~100×, the tropics **reverse to sinking over the 1500 K
    equator**, and the poles rise. The aspect ratio says this atmosphere's vertical velocity
    should be ~(300 km/10 000 km)·v ≈ 0.09 m/s, which is what the control has. **And Ψ_max,
    the number item 27 flipped the default on, reports 157 555 against 153 301 throughout**:
    it is built from `v` alone. Item 26's lesson recurring one item later — an instrument
    that cannot see the component going wrong.

    **THE MECHANISM, AND WHY NOTHING OPPOSED IT.** `initBalancedState` integrates
    `dp_dyn/dθ` level by level and removes each level's mean, so the field's *shape* still
    follows w², which grows strongly upward. Measured at the equator: p_dyn goes 6.9e6 hPa
    at 9 km to 1.79e7 at 55 km, i.e. `−(1/ρ)dp/dz` of −630 to −1290 m/s², downward, where
    the flow is indeed now sinking. With the **default switches the radial equation contains
    nothing else at u = v = 0**: `coriolis_rad` carries the factor `nontrad` and
    `AtomUtils::coriolis_nontraditional()` is false by default; the `−(v²+w²)/r` curvature
    term sits inside `if(AtomUtils::metric_curvature())`, also false; and the buoyancy term
    is multiplied by `buoyancy_ramp`, which is **0 at iteration 0**. So `rhs_u` is
    `−dp_dyn/dr·exp_rm` and nothing else, and any radial gradient the balance writes is a
    pure unopposed vertical force.

    **The same switch touches the θ half, but this is the smaller half of the finding and
    was overstated when first written here.** `metric_curvature()` also governs the `w²cotθ`
    term the item-27 balance is built on, so with the defaults that term is not in `rhs_v`
    either — **the balance was derived from the RHS as written rather than as configured**,
    which is the family's constant trap one level up: not an Earth number this time, an
    Earth-tested code path read without checking whether it runs. **Measured, though, it is
    worth 0.3 %**: rms `F_θ` is 14.61 with the switch off and 14.66 with it on, because both
    balances are ~99.7 % Coriolis. So the defect that matters is the one-component
    balancing, not the missing term — and enabling the switch does **not** repair it, see
    below.

    **Enabling `ATOM_METRIC_CURVATURE` does not make the two components consistent** — the
    obvious hypothesis, tested and refuted at the initialisation:

    | | rms `F_r` | mode 2 leaves | rms `F_θ` | mode 2 leaves |
    |---|---|---|---|---|
    | curvature off | 0.000 | 0.712 | 14.61 | 5.658 |
    | curvature on | **0.0754** | 0.718 | 14.66 | 5.693 |

    The radial equation does gain a real force, exactly as the switch promises. It is
    **~4000× too small** to justify the radial gradient the θ-balance implies (0.075 against
    the 293 that mode 1 creates), so the solve returns essentially the same field either way.
    The θ-force's shape varies with height because `w` does, and **no radial force of
    comparable size exists to pay for that variation** — which is the real reason a
    one-component balance cannot be rescued by completing the metric.

    **THE FIX.** No `p(r,θ)` satisfies both components exactly — that needs
    `∂F_θ/∂r = ∂F_r/∂θ` and the curl of the force field is not zero; the rotational part is
    what buoyancy balances in a real atmosphere, and here buoyancy is ramped off at iteration
    0. So take the best compromise instead of pretending one exists, minimising the
    acceleration the model is left holding in its own units:

    ```
    J = SUM_ij [ (dp/dr*exp_rm - F_r)^2 + (dp/dtheta*inv_rm - F_theta)^2 ]
    ```

    whose Euler–Lagrange equation, `d/dr[exp_rm² dp/dr] + d/dθ[inv_rm² dp/dθ] =
    d/dr[exp_rm F_r] + d/dθ[inv_rm F_θ]`, is solved in conservative form by
    alternating-direction line relaxation. `F_r` and `F_θ` are read from the RHS **term by
    term and switch by switch**, so the balance follows the model's configuration instead of
    assuming one: with the defaults `F_r = 0` and the solve returns what the radial equation
    asks for, which is that p_dyn must not vary with r. The gauge is now a **single global
    constant** — a per-level mean is itself a radial gradient.

    **What each mode leaves behind, which is the diagnostic whose absence let this through:**

    | | radial residual (rms) | meridional residual (rms) | max \|dp_dyn\| |
    |---|---|---|---|
    | no balance at all | 0.000 | 14.61 | — |
    | mode 1 (item 27) | **293.4** (max 1827) | 0.062 | 455.07 |
    | mode 2 (new) | **0.712** (max 2.37) | 5.66 | 199.3 |

    Mode 1 buys a 0.4 % meridional balance by **creating** a radial imbalance of 293 rms
    where there was none. Mode 2 leaves 412× less of it and still removes 61 % of the
    meridional force. Both modes now print this pair, so neither can be judged again by a
    diagnostic that looks at one component. `ATM_BALANCED_MODE=1` restores item 27's field
    verbatim — it reproduces `max |dp_dyn| = 455.06686` exactly, so the legacy path is
    untouched.

    **The solver.** Serial on purpose: 41×181 = 7421 unknowns cost milliseconds, and a
    deterministic sweep order keeps the initial state bit-identical run to run, which a
    red-black OpenMP sweep would not (item 18 is about exactly that). Point SOR crawls on
    this operator — `exp_rm²/dr²` runs ~55× the `inv_rm²/dθ²` term — so it is line relaxation
    in both directions with over-relaxation. Sweeps to a 1e-10 relative residual against ω:
    **1.0 never gets there, 1.90 takes 4265, 1.98 takes 1247, 1.99 takes 2441, 1.995 never**.
    The converged field is identical to four digits across all of them, as it must be; ω
    buys wall clock, not an answer. 0.87 s at ω = 1.98.

    **Two things that did not work, recorded so they are not retried.** A finite-difference
    continuity check on the vtk output cannot distinguish the two runs: the residual scales
    with the field in both (relative residual 0.79 unbalanced against 1.0 balanced), because
    centred differences on a grid whose spacing runs 1.2 → 44 km are not accurate enough for
    the test. And the zonal term is not the missing partner either — `u` varies only ~2 %
    with longitude in both runs, so there is no `∂/∂φ` convergence to balance a radial term
    that is 36× the meridional one.

    **Still open, and deliberately not fixed here:** the `u-v-Cell` glyph vector
    (`Paraview_Atm.cpp:674`) is written in raw m/s onto points laid out in *index* units,
    `x = i*0.1` (level) and `y = j*0.05` (degree). One plot unit is 12–440 km of height
    against a fixed 2222 km of latitude, so an arrow whose true in-plane tilt is 27–70° is
    drawn at 0.2–3.8° and every arrow lies flat along the latitude axis. **No single glyph
    scale factor in ParaView can fix it**, because the distortion runs 150× at 4 km to 20× at
    133 km; the writer has to emit the vector in plot units. A second, smaller one found
    beside it: at the iteration-0 write `t_ref_level` has not been sized yet — only
    `computeLevelMeanTemperature()` fills it, from inside the RK4 step — so `ThermoAtm.h:521`
    takes its fallback `t_ref = 1.0`, and the iteration-0 `BuoyancyForce` field is
    `1e-3·ρ·g·(T/273.15 − 1)`, ~300× too large and shaped like the temperature rather than
    like an anomaly. It is not a force map until the first written iteration.

    **THE 200-ITERATION A/B** (moist physics from iteration 0, κ_H2O = 0.010, `im` = 41,
    `dt_visc` = 1e-4, identical but for the balance):

    | iteration | mode 2 Ψ_max | max \|u\| | mode 1 Ψ_max | max \|u\| | unbalanced Ψ_max | max \|u\| |
    |---|---|---|---|---|---|---|
    | 20 | 156 059 @15° | 0.119 | 157 584 @15° | 2.53 | 153 301 @15° | 0.106 |
    | 100 | 148 411 @15° | 0.173 | 157 054 @15° | 8.08 | 134 484 @15° | 0.086 |
    | 200 | 138 531 @15° | **0.215** | 157 555 @15° | **11.17** | 255 670 **@45°** | 0.095 |

    **Mode 1 reproduces the item-27 run to five digits** — 157 584 / 2.5291 at iteration 20,
    11.1670 against 11.1662 at 200, a 1e-5 difference that is item 18's reduction-order noise
    — so the legacy path came through the refactor untouched, and so did item 27's headline
    numbers (control cell −27 %, drift 255 670 @ 45°, mode-1 return branch 30 771, Ψ₀@45
    = −504).

    **What mode 2 fixes:** the vertical wind, 11.17 → 0.215 m/s, a factor of 52 and inside
    the range the aspect ratio predicts (~0.09 m/s); and the drift, which never takes over —
    Ψ_max stays at 15° for all 200 iterations and the 45° surface value is +1271 against the
    control's 255 670. Its mid-latitude return branch is stronger than mode 1's, 42 254
    against 30 771.

    **What it does not fix: the cell still decays**, 156 059 → 138 531, −11.2 %. That is the
    price of removing 61 % of the meridional force rather than all of it, and the two numbers
    match: mode 2 removes 61.2 % of the force and 60.6 % of the decay rate (−4.9 % against
    the control's −12.3 % at iteration 100). **The decay rate is proportional to the
    unbalanced θ-force**, which is a stronger statement than either measurement alone and
    says item 26 identified the whole mechanism. The tropical return branch erodes fastest of
    the three, −11 413 → −5 946 (−48 %) against mode 1's −27 % and the control's −18 %.

    **`ATOM_METRIC_CURVATURE` does not rescue it either — measured, not assumed.** Three more
    200-iteration runs, and the switch was demonstrably live (the balance's rms `F_r` reads
    7.539e-02 against exactly 0.000e+00 with it off):

    | at iteration 200 | Ψ_max | max \|u\| | Ψ₀ @45° |
    |---|---|---|---|
    | mode 2 | 138 530 @15° | 0.215 | 1271 |
    | curvature + mode 2 | 138 501 @15° | 0.217 | 1298 |
    | curvature + divergence + mode 2 | 138 263 @15° | 0.216 | 1261 |
    | unbalanced | 255 670 @45° | 0.095 | 255 670 |
    | curvature + unbalanced | 256 606 @45° | 0.091 | 256 606 |

    Every pairing agrees to **0.02–0.4 %**, including with `ATOM_METRIC_DIVERGENCE` on as
    well, which is how `Utils.h` says the two halves of the metric belong. **Including on the
    component the terms are supposed to act on**: `vw·cotθ/r` is described there as the
    meridian-convergence term that conserves angular momentum, "the one that makes jets", and
    `max|w|` at iteration 200 is 21.3974 with it off against 21.3929 with it on — 0.02 %. The
    terms are in the equation and are doing nothing at this amplitude.

    **So the decay is not a missing metric term and cannot be fixed by completing the
    metric.** Removing more of `F_θ` without reintroducing a radial gradient requires the
    rotational part to be balanced by buoyancy — a thermal-wind temperature perturbation —
    and that is blocked twice over: `buoyancy_ramp` is 0 at iteration 0, and `densities()`
    overwrites `t` with the adiabat every iteration unless `ATM_PROGNOSTIC_T=1`. **The cell's
    decay is therefore the same open task as invariant 3**, which is where item 29 arrives
    from the radiation side as well.

    **One diagnostic added while reading this output.** `[streamfn]` printed `Psi_max` and
    `Psi_min` and never compared them, though invariant 1 says `|Ψ_max|` = `|Ψ_min|` exactly —
    they are one cell reflected. It now prints the N–S asymmetry as a number and flags
    anything above 1e-3 (reduction order alone gives ~1e-5; the runs here sit at 4.5e-6).
    Item 13's 32 % Hadley asymmetry was found by reading a plot, after two printed extrema
    had disagreed for many runs with nobody subtracting them.

29. **The κ scan, run to 200 iterations: the converged OLR does not move with the opacity
    (done — and it is the answer item 25 feared).**

    Four runs at the shipped defaults, identical but for `kappa_H2O`, over a **64× span**:

    | κ_H2O | OLR @20 | OLR @200 | imbalance @200 | σT_lid⁴ | T_lid | albedo |
    |---|---|---|---|---|---|---|
    | 0.0025 | 334.08 | 272.44 | −1.32 | 271.11 | 262.96 | 0.4986 |
    | 0.010 | 333.07 | 272.42 | −1.31 | 271.11 | 262.96 | 0.4986 |
    | 0.040 | 329.18 | 272.37 | −1.25 | 271.11 | 262.96 | 0.4986 |
    | 0.160 | 315.98 | 272.17 | −1.05 | 271.11 | 262.96 | 0.4986 |

    **A 64× change in the water opacity moves the converged OLR by 0.10 %.** σT_lid⁴ is
    271.11 in all four, identical to five digits, and `t_skin` is 262.96 in all four. The
    5.4 % spread that exists at iteration 20 decays away: **κ changes how fast the column
    relaxes onto the skin temperature, not where it lands.** The single equatorial column
    behaves the same way, 273.7 → 273.1 W/m² across the same 64×, with `σT_surf⁴` identical
    at 275 680 because the surface temperature is prescribed.

    So item 25's discriminating test comes back the way it feared. **The grey scheme is
    reporting its top boundary condition. Item 11's rewrite postponed item 10's defect past
    iteration 20 rather than curing it**, and no claim about `kappa_H2O` is supportable in
    either direction — including the one this file used to make in the other direction. The
    mechanism is not in doubt: `updateSkinTemperature` targets `((absorbed SW + geothermal)/σ)^¼`,
    which contains no κ at all, and `ThermoAtm.h:1446` pins the profile's top at exactly
    `max(t_skin, T_ad)`. Once the effective radiating level sits inside that isothermal skin,
    OLR = absorbed is arithmetic.

    Two smaller things it settles. **Item 24's unexplained latitudinal sign was a transient**:
    at 200 iterations the equatorial column emits 273.7 against a global mean of 272.42, so
    the hottest surface is no longer under the least-emitting column and there is nothing left
    to explain. And **the claim was in the file before the run was**: CLAUDE.md carried "an 8×
    change in `kappa_H2O` moves it 0.9 % (item 25)" and item 27 carried "the κ scan showed the
    outgoing flux does not move with the opacity", both citing an item 25 whose own text says
    the test was "named, not yet run", with no scan in the repo record. The measurement now
    exists and is *stronger* than what was claimed — 0.10 % over 64× rather than 0.9 % over
    8× — which is luck, not vindication. **A cross-reference is not a check, and neither is a
    conclusion written ahead of its measurement.**

## Remaining work


- **Every 400-iteration number predates item 22, and every imbalance predates item 25.** The
  composition-ordering fix halved the OLR at 20 iterations, and the 200-iteration runs then
  showed that a 20-iteration imbalance is a transient. Item 18's table is a valid
  anelastic-versus-Boussinesq comparison and a stale radiation budget; **quote no radiation
  figure that was measured at 20 iterations, and re-measure anything from before item 22.**
  The longest post-item-22 runs are the two 200-iteration ones of items 24-25, and both are
  at `im = 41`.
- **The prescribed profile now has a switch and an adjustment behind it, and the flip has
  not been made.** `ThermoAtm::densities()` still re-imposes the adiabat + isothermal top on
  `t` by default, so what the dynamics and the radiation compute is overwritten before it
  can matter — invariant 3 in CLAUDE.md says radiation must *set* the profile, and the
  prescription still wins. Items 19–20 built the two things the flip needs:
  `ATM_PROGNOSTIC_T=1` and a grid- and cp-aware `ConvectiveAdjustment`. Item 23 measured the
  pair on the corrected baseline and **the case for flipping it has weakened, not
  strengthened**: prognostic now gives −226 W/m² against the prescription's −66, and the
  whole gap is the isothermal `t_skin` lid the prescription pins above 243 km. The next step
  is therefore not the flip — it is deciding what should set the top 60 km when nothing pins
  it, since the honest reading is that −66 W/m² is partly an assumption's doing.
- **The κ scan is done and the OLR is not a result** (item 29). 64× in `kappa_H2O` moves the
  converged OLR by 0.10 %, σT_lid⁴ is 271.11 in every run, and item 23's latitudinal sign
  turned out to be a transient too. **What replaces it is the top boundary condition itself.**
  The radiation cannot be made to respond to anything while `ThermoAtm.h:1446` pins the
  profile's top at `max(t_skin, T_ad)` and `updateSkinTemperature` solves `t_skin` from a
  budget with no κ in it. Breaking that loop — deciding what sets the top 60 km when nothing
  pins it — is now the single most informative thing to do, and it is the same task the
  prescribed-profile bullet above describes. Until it is done, **every OLR number in this
  file is a statement about `t_skin`**, and a grid-convergence check on the OLR would only be
  measuring how well the grid resolves a boundary condition.
- **The `c ≤ 1 − co2` ceiling has stopped firing** (item 23) — zero cells at 20 iterations,
  because it was reporting the misplaced cloud deck of item 22 and nothing else. With the
  transport exonerated (item 17), the column anchor fixed (item 19) and this gone, **no known
  mass sink remains**. The untested case is a run past `moist_phys_start_iter = 300`, where
  sedimentation runs for the first time; if the count comes back, measure *where* before
  proposing *why*, which is the one thing the two refuted mechanisms have in common.
- **The tropical cell decays 11.2 % over 200 iterations even when the balance is right**
  (item 28), and that is now the open question about the circulation. The decay rate is
  proportional to the unbalanced θ-force: mode 2 removes 61 % of the force and 61 % of the
  decay. Removing the rest without reintroducing a radial gradient needs the rotational part
  balanced by buoyancy — a thermal-wind temperature perturbation — which `buoyancy_ramp` = 0
  at iteration 0 and `densities()`'s overwrite of `t` both block. **Same task as the
  prescribed profile above**, approached from the dynamics instead of the radiation.
  `ATOM_METRIC_CURVATURE` was tested as the alternative and changes nothing (0.02–0.4 %).
- **The tropical return branch erodes fastest with the balance on**: −11 413 → −5 946
  (−48 %) against mode 1's −27 % and the control's −18 %. Item 27 first saw this and the
  two-component balance slightly worsens it. Unexplained.
- **The `p_dyn` ceiling is 2000 and is untested beyond 200 iterations** (item 27). It
  replaces a backstop that existed because unbounded `p_dyn` blew this solver up on Earth;
  that failure was orography-driven and should not arise over a featureless surface, but that
  is an argument, not a measurement. A long run is the check. Note item 28 halves what the
  balance asks of it — max |dp_dyn| 455 → 199.
- **The `u-v-Cell` glyph vector is drawn in the wrong metric** (item 28), so the meridional
  circulation is invisible in ParaView even when it is there. `Paraview_Atm.cpp:674` writes
  raw m/s onto points laid out in index units; the fix is to emit the vector in plot units
  (`u/(dz_i/0.1)`, `v/(dy/0.05)`). A ParaView-side scale factor cannot substitute, because
  the distortion runs 150× at 4 km to 20× at 133 km. The same applies to the `longal` writer.
- **The iteration-0 `BuoyancyForce` field is not a force** (item 28): `t_ref_level` is sized
  only inside `computeLevelMeanTemperature()`, called from the RK4 step, so the first write
  takes the `t_ref = 1.0` fallback at `ThermoAtm.h:521` and reports ~300× too much, shaped
  like the temperature. Either size the array at init or make the diagnostic say so.
- **The tropical return branch still erodes** — −11 969 → −8 689 over 200 iterations, and
  marginally faster with the balance on than without it. The drift that used to mask this is
  gone, so this is now the open question about the cell itself.
- **`dt_visc` no longer matches the committed grid.** The config ships 4e-5, validated for
  `im = 61`; item 24 validated 1e-4 at the `im = 41` that is now the default, worth 2.5×.
- **The run is not converged, and 400 iterations is not close** (item 18). The meridional
  wind is in free acceleration under an unopposed Coriolis torque — the pressure gradient
  that should balance it is at 1.8 % of it after 400 iterations and growing linearly, which
  puts geostrophic adjustment of order 10⁴ iterations away. Ψ_max grows linearly throughout.
  **This is the largest open question about the dynamics**, it is common to both continuity
  formulations, and item 18 ruled out the elliptic solver as its cause.
- **The OLR is not grid-converged either**: 519 W/m² at a 260 km shell against 581 at
  300 km, with `im` fixed at 61. Both figures predate item 22, so the check has to be redone
  as well as extended — refine vertically and repeat.
- **`moist_phys_start_iter = 300`** means a 400-iteration run is dry for three quarters of
  its length. Deliberate (it lets the circulation form before the stiff microphysics
  starts), but it must be stated whenever a run is quoted — and it is why the 20-iteration
  measurements above are all made with sedimentation and the ice schemes switched off.
- **Thread-count dependence at ~1e-8 remains** (item 18). The races are gone and a fixed
  thread count is bit-identical run to run, but OpenMP reduction order still moves the last
  digit, and it re-enters the physics through the global means the column is rebuilt from.
  Ordered reductions, `t_skin` first.
- **Deep convection is still inactive.** Its trigger thresholds (1000/970/900/800 hPa) are
  absolute Earth surface pressures and never fire at 250 bar. They need to become fractions
  of the local surface pressure — the same repair `initCloudIce` needs below.


- **`albedo_cloud = 0.50` IS the model's planetary albedo** and is an assumption. It is
  now the second-biggest lever after the opacities. A deep, cold, slowly sedimenting
  Hadean deck could plausibly be brighter.
- **`geothermal_flux` is the open number.** The ≥ ~195 W/m² figure was derived with the
  clear-sky albedo and the too-low insolation, so it has to be redone once the radiation
  can produce an OLR worth comparing to. Check against magma-ocean cooling estimates.
- **`initCloudIce`'s H_crit parabola is keyed to absolute pressure** (`p_crit = 1000` hPa,
  `p_mid = 550`), an Earth surface pressure. It should be a fraction of the local surface
  pressure, like the deep-convection triggers. Item 22 showed that what `initCloudIce`
  builds is not overwritten before it reaches the radiation, so this one has a path to the
  OLR.
- **A grey scheme cannot represent the window regions** that set the real runaway limit,
  whatever the three κ are set to.
- **The surface temperature is prescribed, not solved.** Everything above is conditional
  on that.
- **Boussinesq** remains untested against a column whose density spans two orders of
  magnitude (see CLAUDE.md).
