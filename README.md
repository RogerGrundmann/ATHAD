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

> **BEFORE QUOTING ANY OLR OR IMBALANCE FROM THIS FILE (2026-08-21, item 67).** `t_skin` was
> solving `sigma*T^4 = F` where the grey skin relation `sigma*T^4 = F/2` belongs — so the
> transparent lid was assigned the planet's entire energy input as its emission, and "the
> energy balance closes" was an identity rather than a result. The factor is in, and **on by
> default**. Consequence: the model now absorbs 271 W/m² and emits **147** at 200 iterations,
> imbalance **+123.76 and widening**, where before it closed to −1.41.
>
> **Every OLR and imbalance figure recorded below before that date was measured on the old
> branch.** None is retracted — they are correct measurements of a branch still reachable with
> `ATM_SKIN_GREY=0` — but none describes the shipped model. Item 67 carries the paired series
> that does.

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

30. **The radiation solver was never converging, and the two-stream system did not need
    iterating at all (done — solver written, default not yet flipped).**

    Chasing why the prognostic column grows a 1124 K inversion in its top four levels ended
    somewhere else. Those levels are not broken: the scheme's own comment says that at the
    top, where `dn → 0`, the update reduces to `σT⁴ = up/2` — the classical skin temperature
    *of the flux arriving from below* — and the measurement agrees, 82 000 W/m² giving 922 K
    predicted against 889 K measured. **They are slaved to the column's flux and cannot be
    fixed locally**, and pinning them to `t_skin` would re-impose exactly the lock item 29
    diagnosed. What was wrong was upstream.

    **`n_lambda = 4` is an Earth constant.** The Lambda iteration is a Jacobi relaxation on a
    41-link chain — information moves one layer per sweep, so it needs O(N²) — and 4 sweeps
    is nowhere near converged at 250 bar:

    | `n_lambda` | lid T | τ≈1 level | OLR | cost/call |
    |---|---|---|---|---|
    | 4 (inherited) | 889.2 K | 774.6 K | 65 982 | 0.370 s |
    | 40 | 687.3 K | 662.4 K | 50 656 | 0.443 s |
    | 400 | 477.9 K | 349.0 K | 15 248 | 4.551 s |
    | 2 000 | 478.0 K | 362.0 K | 15 754 | 34.5 s |
    | 10 000 | 478.0 K | 362.0 K | 15 769 | 116.3 s |

    A factor of 4 in the headline number, converged only by ~2000 sweeps at 500× the cost.

    **It does not need iterating.** With `a_i = 1 − ε_i/2` and `b_i = ε_i/2`, a layer in
    radiative equilibrium transfers `U_i = a_i U_{i−1} + b_i D_i` and
    `D_{i−1} = b_i U_{i−1} + a_i D_i`, and **`a_i + b_i = 1` exactly** — so `U − D` is the
    same at both faces and the net flux is constant. That is the definition of radiative
    equilibrium, and here it falls out of the discretisation rather than being imposed.
    Eliminating the fluxes leaves a march for the source function that is *linear in F*,
    `B_{i+1} − B_i = (F/2)[(1−ε_i)/a_i − 1/a_{i+1}]`, closed by `D_top = 0` into
    `F = σT_s⁴/(1/(2a_n) − G_n)`. **Two O(N) passes, exact, no linear solve** — Feautrier's
    tridiagonal is not needed either.

    | method | lid T | τ≈1 level | OLR | cost/call |
    |---|---|---|---|---|
    | `n_lambda` = 10 000 | 478.0 K | 362.0 K | 15 769.09 | 116.332 s |
    | **direct** | 478.0 K | 362.0 K | **15 769.11** | **0.036 s** |

    Eight significant figures against 10 000 sweeps, **3200× faster** than converging them
    and **10× faster than the wrong default**. Validated against the iteration run to
    convergence on synthetic transparent, thick and ATHAD-graded columns, and against both
    exact limits: `ε → 0` returns `F = σT_s⁴` and `T = T_s/2^¼`, the classical skin
    temperature of a freely radiating surface.

    **What it does not change is the useful part.** Over 200 iterations in the standard
    prescribed configuration the OLR goes **273.31 → 272.91 W/m²**, the imbalance −2.23 →
    −1.83: **0.15 %**. Every OLR figure this project has quoted survives the solver defect,
    because while `densities()` re-imposes the adiabat every iteration the reported flux is a
    property of the *prescribed profile* and does not depend on how well the radiation was
    solved.

    **And it confirms item 29 rather than overturning it.** Converged, prescribed, 200
    iterations: **272.91 at κ = 0.010 against 272.54 at 0.160 — 0.14 % over 16×.** The
    κ-insensitivity was never a solver artefact. The prescription is the lock, and that now
    rests on a converged calculation.

    **One thing it does change.** With the converged solver the prognostic column stops
    sitting on a hot branch and starts relaxing: OLR **15 764 → 8 542 → 6 211** over 200
    iterations, still falling, against ~66 000 and motionless at 4 sweeps. Still 23× the
    absorbed flux, so nothing is settled — and the earlier reading of this run, *"58× the
    absorbed flux, the column is too transparent"*, was a snapshot of a decaying quantity.
    **That is the trap item 25 exists to warn about, walked into twice in one day**: first
    predicting the hot branch was physical, then over-correcting to "it converges to the
    prescribed values" on a monotone trend that had not converged. The rule earns restating —
    *a monotone trend is not a limit* — and it applies to the 6 211 above as much as to
    anything else here.

31. **Earth's cell latitudes are wrong for a 5.5-hour day, and about a third of the item-28
    decay was the model rejecting them (done).**

    `VelocityInitializer` prescribes Hadley at 15°, Ferrel at 45°, polar at 75° — Earth's
    latitudes, which follow from Earth's thermal Rossby number. This atmosphere's is 12.5×
    smaller:

    | | scale height | Δθ/θ | Ω | Ro_T | Held–Hou edge |
    |---|---|---|---|---|---|
    | Earth | 8.4 km | 45/288 = 0.156 | 7.29e-5 | 0.0598 | 18.1° |
    | ATHAD | 59.3 km | 50/1500 = **0.033** | 3.17e-4 | **0.0048** | **5.1°** |

    Rotation is 4.35× faster (Ω² 18.9×) and the **fractional** contrast 4.7× weaker, only
    partly offset by a 7× deeper atmosphere. **A hot surface is not a strongly
    *differentially* heated one** — which is why the higher energy content narrows the
    circulation instead of widening it, the giant-planet direction rather than the Venus one.
    The Rhines scale agrees: ~2450 km against Earth's ~3490, so ~8 bands pole-to-pole
    against ~6.

    Four 200-iteration runs, mode 2 balance, moist physics from iteration 0, at what is now
    the `cell_lat_scale` config parameter (it was the `ATM_CELL_LAT_SCALE` environment knob
    when these were run — see item 32):

    | scale | Hadley anchor | Ψ @20 | Ψ @200 | decay |
    |---|---|---|---|---|
    | 1.00 | 15° (Earth) | 156 059 | 138 530 | **−11.2 %** |
    | 0.75 | 11° | 148 752 | 134 937 | −9.3 % |
    | 0.50 | 7° | 132 100 | 122 337 | −7.4 % |
    | 0.33 | 5° | 106 309 | 98 519 | **−7.3 %** |

    The decay falls monotonically as the cell moves to where the regime wants it and
    **saturates between 7° and 5°, where Held–Hou puts the edge**. So about a third of item
    28's decay was the imposed width, and the remaining ~7.3 % is the residual meridional
    force — consistent with mode 2 removing 61 % of it. Two separate causes, now separated.

    **The Ferrel and polar anchors are a different matter and are not transferable at all.**
    Earth's Ferrel cell is thermally indirect and eddy-driven; this flow is **axisymmetric to
    2 %** (measured: `u` varies 2.3 % across all 361 longitudes) and the column is neutrally
    stratified by construction (`cosmo_lapse_fraction = 1.0`, so N² ≈ 0), so there are no
    baroclinic eddies to drive one. Those two cells have no maintenance mechanism in **any**
    eon in this model, which is a model property rather than a Hadean one and should not be
    confused with the regime argument above.

    The 50 K contrast is **prescribed**, so this rests on an input — but robustly: recovering
    Earth's `Ro_T` at this rotation rate would need ΔT ≈ 630 K, and even 200 K still lands
    at ~10°.

32. **The cell latitudes become a config parameter and the default stops being Earth's
    (done).**

    Item 31 left the finding in an environment variable, which is the wrong place for it:
    an `ATM_*` knob is invisible in the committed configuration, absent from the diff, and
    silently defaults to the value the item argued against. `ATM_CELL_LAT_SCALE` is now
    **`cell_lat_scale`** in `param.py`, **default 0.33**, generated into both
    `config_athad.xml` files like every other parameter. The derivation and the four-run
    table live in the `param.py` comment, so the number and its justification are in the
    same place.

    `VelocityInitializer::latScale()` was a function-local `static` reading `getenv` once
    per process; it is now `m.cell_lat_scale`, which also removes a piece of hidden global
    state from a class that otherwise takes everything from the model. The anchor line now
    prints on **every** run rather than only when the scale was non-Earth — with a
    non-default default, a diagnostic that stays silent at 1.0 would report Earth's
    latitudes by saying nothing:

    ```
    cell_lat_scale = 0.330: Hadley anchor at 5 deg, Ferrel at 15 deg, polar at 25 deg
    cell_lat_scale = 1.000: Hadley anchor at 15 deg, Ferrel at 45 deg, polar at 75 deg   <- Earth's latitudes
    ```

    Verified at both values on a 1-iteration run: the anchors print as above and the fields
    genuinely differ (max `w_d` 0.415 against 0.969 m/s, the `v_d` extremum moving from 86°
    to 75°), so the parameter reaches the dynamics rather than only the printout.

    **Two things this does not settle.** The scaling is applied to all three anchors because
    it is one map on the latitude index, but only the Hadley one has a regime argument behind
    it — item 31's point that the Ferrel and polar cells have no maintenance mechanism in this
    model is unaffected, and scaling their anchors remains meaningless rather than
    justified. And **the default change rests on a single 200-iteration set**, in which the
    measurement cannot separate 0.33 from 0.50; the Held–Hou estimate is what chooses. Every
    number in this file measured before this item used 1.0.


33. **Scaling the cell latitudes without the amplitudes tripled the shear; a mode now
    scales both (done, default off).**

    Item 32 moved the anchors and nothing else, so the same velocity change was
    interpolated across a third of the latitude span and every meridional gradient in the
    initial state was 1/s times Earth's — 3× at the 0.33 default. `cell_amp_mode`:

    | mode | v × | w × | rationale |
    |---|---|---|---|
    | 0 | 1 | 1 | off — what item 32 committed |
    | 1 | s | s | continuity: preserves every prescribed gradient |
    | 2 | s | s² | continuity + angular momentum for the jet |

    Mode 2's s² is the small-angle form of the angular-momentum jet ratio
    `sin²φ/cos φ`; exact at s = 0.33 is 0.110 against s² = 0.109. The RADIAL amplitudes
    never scale in any mode, because continuity makes `du_r/dz` invariant when v and the
    latitude scale shrink together.

    Three 100-iteration runs, mode 2 balance, moist physics from iteration 0, measured on
    the *untiled* layout — see item 36, which invalidates the absolute numbers but not the
    comparison between modes:

    | mode | Ψ@20 | Ψ@100 | decay | max w | mean KE |
    |---|---|---|---|---|---|
    | 0 | 106 309 | 102 800 | −3.30 % | 22.60 | 19.50 |
    | 1 | 34 824 | 33 666 | −3.33 % | 7.456 | 2.124 |
    | 2 | 34 982 | 34 504 | **−1.37 %** | 2.459 | 0.244 |

    **The decay is scale-invariant but ratio-sensitive, and that isolates its cause.**
    Modes 0 and 1 differ by 3× in *every* amplitude and decay identically. Modes 1 and 2
    share v and Ψ and differ only in the jet — and the decay halves. Since w/v goes
    8.29 → 8.29 → 2.92, the residual tropical-cell decay tracks the prescribed **zonal
    jet**, not the overturning or the overall scale. That is a handle on item 28's open
    question and points where the README already suspects: an unopposed Coriolis torque,
    which a weaker jet exerts less of.

    Default stays 0. One run, and mode 2's absolute winds (max zonal 2.46 m/s) are far
    below the model's own emergent scale.

34. **`u_0` = 8 m/s is Earth's mean surface wind, and it cannot be calibrated in this
    configuration.**

    `u_0` non-dimensionalises the whole momentum equation. Term by term most of it is
    genuinely invariant — Coriolis `omega*L/u_0`, drag, eddy viscosity and Held-Suarez all
    scale so that the ratio to advection is fixed. **Buoyancy is the exception**: coded as
    `g*dt/u_0` where the advective-time-consistent form is `g*L/u_0²`, and its own comment
    admits it is "the calibrated, STABLE reference" rather than a derived one. `u_0` also
    sets the clock: physical time per iteration is `dt*L_atm/u_0` and `dt` is a constant
    independent of `u_0` (`cAtmosphereModel.cpp:1127`).

    Predicted therefore that `u16@100 ≈ u8@50`. **The two diagnostics disagree:**

    | | u_0 = 8 | u_0 = 16 | ratio | clock predicts |
    |---|---|---|---|---|
    | Ψ decay over 100 iter | −3.30 % | −1.83 % | **1.80** | 2.0 ✅ |
    | KE decay over 100 iter | −14.56 % | −14.20 % | **1.03** | 2.0 ❌ |

    The circulation behaves like a clock in `u_0`; the kinetic energy does not. Mean T is
    identical to six figures across every run, which is what `densities()` re-imposing the
    adiabat each iteration would produce — so the energetics do not experience elapsed
    time at all. **`u_0` cannot be calibrated from a prescribed-profile run**, which puts
    it behind the same blocker as items 25, 29 and 30. Not changed, and on present
    evidence should not be until `ATM_PROGNOSTIC_T` works.

    Chasing the buoyancy coefficient turned up a defect of item 30's shape. RK4 multiplies
    every RHS by `dt`. The Held-Suarez block carries an explicit note that its extra `*dt`
    was removed because "the extra `*dt` made HS enter as dt² = ~1e-4 too weak -> inert".
    **The buoyancy term (`RHS_Atm_Turb.cpp:1002`) and the Rayleigh drag (`:1019`) still
    carry it.** Both enter as dt², everything else as dt¹; for buoyancy that is a factor
    `L/(u_0*dt)` ≈ 2e7 below the consistent coefficient. Two caveats before this is called
    a one-line fix: the same comment records that a merely 336× larger coefficient "drove
    a polar vertical runaway", so the solver has already failed to carry a value four
    orders of magnitude short of consistent — that points at the Boussinesq open risk. And
    it makes both terms dt-dependent, so item 24's `dt_visc` 4e-5 → 1e-4 silently changed
    the effective buoyancy and drag by 6.25×.

35. **Earth's grid indices as `switch` labels: four of seven radial anchors were silently
    dropped at any scale ≠ 1 (fixed).**

    `compute()` passed `init_u` the SCALED index; `init_u` dispatched on `switch (j)` with
    Earth's UNSCALED indices as case labels. At `cell_lat_scale` = 1.0 they coincide, which
    is why this was invisible until item 32 flipped the default.

    `js()` clamps on its input, so `js(0)` and `js(jm-1)` return the poles unscaled: the
    poles and the equator always hit their case correctly. The **same four** anchors were
    dropped at every s < 1 — ±30° and ±60°, the ascent/descent branches at the
    Hadley/Ferrel and Ferrel/polar boundaries. `form_diagonals` then interpolated across
    anchors that had never been written. Ψ@20 moves 106 308.91 → 106 012.76 (−0.28 %) when
    fixed, so every item-31 row except `s = 1.00` shifts by about that.

    The anchors are now generated from `n_cells_hemisphere` rather than written as ~30
    latitude literals plus 24 hand-written `form_diagonals` calls, and the amplitudes are
    expressed by role. Verified bit-identical at n = 3.

    Two process notes, because the second is the interesting one. The commit message for
    the fix (`eaf6e8d`) states the poles were sign-flipped; they were not, for the `js()`
    clamping reason above. **The same oversight then produced the same bug in the fix**:
    `jN` initially inlined `round(90 - phi*s)`, moving the pole anchor to 30°N and leaving
    everything poleward unfilled. The tell was a 2e-7 residual in Ψ that had no business
    existing, on a change proved equivalent by enumeration. *An equivalence argument is
    not a check* — the same lesson as "a cross-reference is not a check".

36. **The cells did not tile the hemisphere, and item 31's scan was measured on a layout
    where narrowing meant deleting (fixed).**

    `cell_lat_scale` compressed every anchor and let `js()` pin the pole. At s = 0.33 that
    gives edges 0/10/20/90 — cells **10, 10 and 70 degrees wide**, with the entire
    extratropics a single linear ramp. Raising n made it worse: n = 4 packed four narrow
    cells into 0–22° and left a 68° cell whose prescribed core sat at 26°.

    **This was the committed behaviour from item 31 onward**, so item 31's four-run scan
    compared configurations that differed in more than cell width: narrowing mostly
    *removed* the outer cells. The Ψ_max figures stand — that peak is at 5°, inside the
    resolved part — but "the decay saturates at 0.33" was concluded from it.

    Held-Hou constrains the direct cell's edge and says nothing about the extratropical
    bands, which the Rhines scale sets independently. So s now scales the Hadley edge only
    and the remaining band is tiled by the other n-1 cells. This is the first layout on
    which the two estimates agree — extratropical bands 27° at n = 4 and **20° at n = 5**,
    against the Rhines ~22°:

    | n | s | edges | widths |
    |---|---|---|---|
    | 3 | 1.00 | 0/30/60/90 | 30/30/30 — **Earth exactly** |
    | 3 | 0.33 | 0/10/50/90 | 10/40/40 |
    | 5 | 0.33 | 0/10/30/50/70/90 | 10/20/20/20/20 |

    The residual 2:1 between Hadley and band is regime physics, not a tuning artefact: the
    Held-Hou edge goes as 1/Ω and the Rhines band only as 1/√Ω, so fast rotation separates
    two scales that nearly coincide on Earth (30° against 31°) — which is why Earth's
    three cells look uniform and this atmosphere's cannot.

    **Two 100-iteration runs, n = 3 against n = 5** (`cell_lat_scale` 0.33, moist physics
    from iteration 0). Every prescribed core is a Ψ maximum and every interior edge a Ψ
    minimum, so the structure is laid down as intended.

    - **Cell count does not reach the tropical cell.** Decay 20 → 100 is −3.48 % at n = 3
      against −3.41 % at n = 5. The Hadley edge is 10° either way.
    - **The extratropical cells GROW rather than decay** — +11 to +25 % at the mid-latitude
      cores over 20 → 100 while the tropical cell loses 3.5 %. This was predicted to go the
      other way, on the grounds that nothing maintains an indirect cell here.
    - **n = 5 keeps its boundaries; n = 3 does not.** Core-to-edge contrast at n = 3's single
      interior edge falls 24:1 → 10:1 as the two 40° cells merge, while n = 5 *deepens* two
      of its three interior boundaries (4.2:1 → 5.3:1 and 4.5:1 → 6.4:1).

    **Mean KE is 38–40 here against 19.50 on the untiled layout.** Restoring the
    extratropics roughly doubled the model's kinetic energy, so every KE figure measured
    before this item — item 34's `u_0` comparison included — was made with half the
    circulation missing.

    n = 5 is now better supported than n = 3 on three independent grounds, but the default
    is unchanged at 3: 100 iterations is a trend, and the extratropical growth is exactly
    the kind of monotone trend this file has twice been caught extrapolating.


37. **400 iterations: item 36's boundary finding was a transient and reverses, but n = 5
    survives on a different argument.**

    Two 400-iteration runs, n = 3 against n = 5, otherwise identical to the 100-iteration
    pair of item 36 (`cell_lat_scale` 0.33, `cell_amp_mode` 0, moist physics from
    iteration 0, `checkpoint` 20). Ψ@20…Ψ@100 reproduce item 36 exactly.

    **The tropical cell is linear and n-independent.** Ψ_max falls ~470 per 20 iterations
    with no turnover anywhere in the trace:

    | | Ψ@20 | Ψ@200 | Ψ@400 | 20→400 |
    |---|---|---|---|---|
    | n = 3 | 106 788 | 98 457 | 88 930 | −16.7 % |
    | n = 5 | 106 511 | 98 393 | 89 245 | −16.2 % |

    The two converge on each other (0.35 % apart at 400). Cell count does not reach the
    Hadley cell at 20, 100 or 400 iterations.

    **Item 36's boundary claim does not survive.** It reported that n = 5 *deepens* two of
    three interior boundaries while n = 3 erodes its only one. At 400 every boundary in
    both layouts has eroded to the same ~3:1 core-to-edge contrast:

    | | edge | 100 iter | 400 iter |
    |---|---|---|---|
    | n = 3 | 50° | 10.2 : 1 | **3.2 : 1** |
    | n = 5 | 30° | 5.3 : 1 | 3.4 : 1 |
    | n = 5 | 50° | 6.4 : 1 | 3.0 : 1 |
    | n = 5 | 70° | 6.9 : 1 | 3.2 : 1 |

    The edges fill in faster than the cores grow — +104 to +223 % at the edges against
    +31 to +48 % at the cores over 100→400. **The deepening was an adjustment and it was
    extrapolated**, which is the third time this file has recorded that error and the
    second time in one day. Item 36's own closing sentence named the risk and the write-up
    made the claim anyway.

    **What replaces it is stronger, because it does not rest on a trend.** Peak |Ψ| by
    latitude at iteration 400, cores marked:

    ```
    n=5   20° 56071*  30° 16502   40° 51590*  50° 17470   60° 36193*  70° 11384   80°  6019*
    n=3   20° 61119   30° 57052*  40° 34952   50° 17634   60° 16243   70° 13563*  80°  4233
    ```

    n = 5's maxima sit **exactly on its prescribed cores** at 20/40/60°. n = 3's do not:
    its 30° core is no longer a peak, the maximum having migrated equatorward past 20°. The
    model keeps the layout n = 5 gives it and rearranges the one n = 3 gives it. That is a
    statement about position rather than about a rate, so unlike the boundary contrast it
    is not a quantity still in motion.

    A second signal, common to both: everything equatorward of ~60° intensifies while
    70–90° weakens by 25–31 %. The circulation concentrates into mid-latitudes whatever the
    cell count.

    **Controls.** OLR 272.84 and 272.75 — both have relaxed onto σT_lid⁴ = 271.11, the lock
    of items 25/29/30, and are identical across cell counts. Mean T 809.41 in both. Mean KE
    29.5 / 29.3, down from ~39 at iteration 100. N–S asymmetry −3.7e-06 and −3.5e-06: inside
    the 1e-5 tolerance, but ~60× the 20-iteration value, so it **grows with run length** and
    should be watched rather than assumed.

    **Nothing here is converged and the default stays at 3.** Ψ is still linear at 400, the
    boundary contrast is still falling, and the mid-latitude intensification is still
    monotone — all consistent with the standing estimate that geostrophic adjustment is
    O(10⁴) iterations away. n = 5 is the better-supported layout on three independent
    grounds now (Rhines band width, the geometric artefact, and core position), none of
    which is a converged measurement.

38. **`n_cells_hemisphere` defaults to 5 (done).**

    Items 36–37 built the case; this changes the default from Earth's 3. The committed
    configuration now lays down

    ```
    n_cells_hemisphere = 5, cell_lat_scale = 0.330
    cell edges at 0 / 10 / 30 / 50 / 70 / 90 deg N; cores at 5 / 20 / 40 / 60 / 80 deg;
    widths 10 / 20 / 20 / 20 / 20 deg
    ```

    — a 10° direct cell at the Held–Hou edge plus four 20° extratropical bands against the
    Rhines scale's ~22°, which is the first layout in this model on which the two
    independent estimates agree.

    **What the default rests on, and what it does not.** Three arguments, none converged:
    band width, the removal of the 40°-wide-cell artefact, and — the only one that is not a
    rate — n = 5's Ψ maxima sitting exactly on its prescribed cores at 400 iterations while
    n = 3's migrate off theirs. It does **not** rest on item 36's interior-boundary claim,
    withdrawn by item 37. And it changes nothing about the tropical cell, which decays
    −16.2 % at n = 5 against −16.7 % at n = 3 and does not respond to cell count at any run
    length tested.

    **This is a structural default, not a claim that the model sustains five cells.** The
    caveat item 31 raised about the Ferrel and polar anchors applies to all four
    extratropical cells: the flow is axisymmetric to 2 % and neutrally stratified by
    construction, so no baroclinic eddies exist to maintain an indirect cell. n sets the
    structure the model is handed. Everything in this file measured before item 38 used
    n = 3.

39. **`exp_rm` is not the Jacobian it is documented to be, and `im` is the wrong question.**

    The question this started from was which vertical resolution is physically right, 41
    levels like the rest of the family or 61. The answer is that both are being asked about
    a grid the solver and the radiation do not agree exists.

    **The radiative case first, because it is what `im` actually decides.** Rebuilding the
    column from the model's own formulas — Shomate `cp_of`, the exact `tau_gas` of
    `MultiLayerRadiation.h`, hydrostatic on the layer-mean T anchored at `p_0` — and sorting
    the layers by optical depth measured down from the top:

    ```
                                          im = 41    im = 61
    layers with tau_above in [0.1, 10]          1          3
    dtau of the layer where tau_above = 1     54.8        7.6
    height of that crossing                 201 km     218 km
    max dz/H                                  2.19       1.48
    ```

    At `im = 41` the **entire photosphere is one grid cell**, and that cell has `dtau = 55`:
    level 35 sits at `tau_above = 10.8` and level 36 at 0.82, so the column goes from opaque
    to transparent inside one layer. The two-stream sweep is first order in `dtau`, so the
    escaping flux is computed across a layer whose source function is unresolved. That is
    very likely what the un-converged OLR of the shell-depth check (519 W/m² at 260 km
    against 581 at 300 km) is made of. `im = 61` is not converged either — you want
    `dtau` ≲ 1 at the crossing — but it is poor rather than broken. **ATHAD_COND is the
    kinder case and its shipped 61 is defensible**: `dtau@tau=1` = 2.9 and 4 photosphere
    layers, because its 120 km shell is ~8 surface scale heights against ATHAD's 13.7
    pressure e-folds, so the stretch has far less room to run away.

    **Then the metric, which outranks the count.** `TurbulenceAtm.h` states that
    `exp_rm = 1/(rm+1)` "is the Jacobian of the radial coordinate transformation", and every
    radial derivative in the core is an index difference times `inv_2dr` times `exp_rm`.
    `PressureSolverAtm.h` says the same in equation form: `dp/dr_physical = exp_rm *
    dp/d(rad.z)`. It is not the Jacobian. `exp_rm = 1/(rm+1)` implies `dz/d(rm) = rm+1`,
    a **quadratic** stretch `z ~ (rm+1)²/2`, while `init_layer_heights` builds an
    **exponential** one, `z = (exp(zeta·(rm − rad.z[0])) − 1)·L_atm`. `zeta` appears nowhere
    in the momentum, continuity or transport differences — only in `init_layer_heights`,
    `InitValues_Atm.cpp`, `metricShellLength()` and one turbulence `log()`.

    The test is **unit-free**, which is what makes it conclusive rather than a units
    argument. For `(index difference)·inv_2dr·exp_rm` to be `d/dz_physical` under *one*
    constant length normalisation — whatever that normalisation is — the quantity
    `ratio(i) = [inv_2dr·exp_rm(i)] / [1/(z[i+1] − z[i−1])]` must be constant in `i`.
    Measured at the shipped configuration (`im = 41`, `zeta = 3.0`):

    ```
      i    rad.z    exp_rm    z[km]   dz_true[km]   J_true[m]   J_code[-]   ratio[m]
      1   1.0250   0.49383      1.2        1.320        50830       2.025    25124.7
      9   1.2250   0.44944     15.2        2.404        92618       2.225    41665.1
     17   1.4250   0.41237     40.5        4.381       168761       2.425    69657.3
     25   1.6250   0.38095     86.8        7.983       307502       2.625   117253.6
     33   1.8250   0.35398    171.0       14.546       560306       2.825   198524.2
     39   1.9750   0.33613    277.2       22.813       878734       2.975   295649.9
    ```

    `ratio` rises monotonically by **11.77×**. The true `dz/d(rad.z)` spans 17.3× across the
    shell; `1/exp_rm` spans 1.47×. Physically: **the core's radial length unit is 25 km near
    the surface and 296 km at the top**, where it should be one number. No choice of `L_unit`
    can absorb that — the metric has the wrong *shape*, not the wrong scale.

    **It is not a function of `im`, and refining marginally hurts** — 11.77× at 41 levels,
    12.29× at 61. The spread is set by `zeta` and the `rad.z` span, neither of which the
    level count touches. So the original question has no good answer in the form it was
    asked.

    ```
    spread of ratio(i):   zeta 3.0   2.0    1.5    1.0    0.5   0.405    0.3    0.1
                 im=41       11.77  4.55   2.83   1.76   1.09    1.02   1.10   1.34
    ```

    The minimum at `zeta ≈ 0.405` is not a coincidence: it is `ln(1.5)`, where the quadratic
    form's range over `rad.z ∈ [1,2]` matches the exponential's `exp(zeta)`.

    **Three independent arguments now point the same way — cut `zeta`, do not raise `im`.**
    Photosphere resolution (`dtau@tau=1` 54.8 → 3.4 at `zeta = 1.0`), the diffusive timestep
    (the CFL goes as the physical surface spacing squared, so `zeta = 1.0` at `im = 41` buys
    a **13×** larger `dt` — about 30× less wall clock per unit physical time than `im = 61`
    at `zeta = 3.0`), and the metric (11.8× → 1.0). The reason they agree is the same in each
    case: a surface-refining stretch is an **Earth** choice, made to resolve a boundary layer.
    Here the bottom cells carry `dtau ≈ 9e4` each — radiatively dead, in perfect LTE — and the
    surface temperature is a *prescribed boundary condition*. `zeta = 3.0` spends the grid
    where the model has nothing to learn and starves the one region that sets its only output.

    **What is verified and what is not.** The shape is measured and is not in doubt. The
    *absolute* consequence is not: against the `L_atm = 15719 m` that `RHS_Atm_Turb`'s
    coefficients divide by, `ratio` is 1.6× at the surface and 18.8× at the top, suggesting
    radial derivatives are overweighted everywhere and increasingly so with height — but the
    full non-dimensionalisation has not been traced, so treat that reading as indicative.
    And the two repairs are **not** equivalent: replacing `exp_rm` with the true Jacobian is
    the principled fix and makes any `zeta` valid, but it touches ~91 sites across 10 files
    including the Poisson operator and moves every number this model has produced; cutting
    `zeta` makes the existing metric accidentally correct and is far cheaper. That decision
    wants its own measurement.

    **This is inherited, and it is worse upstream.** `ATOM_Precipitation` has the identical
    `exp_rm = 1/(rm+1)`, the identical exponential `init_layer_heights`, and `zeta = 3.715`
    hard-coded at `cAtmosphereModel.h:273` — a spread of **23.2×**, about twice ATHAD's.
    `ATHAD_COND` is at `zeta = 3.0`, so ~12× like ATHAD. ATJUP/ATSAT/ATURAN/ATNEPT have no
    `exp_rm` at all. Checked, not assumed.

    `cAtmosphereModel::checkRadialMetric()` now prints the spread at every startup and the
    per-level table under `ATM_METRIC_CHECK=1`. It is print-only: a 1-iteration run against
    the pre-change binary differs in timestamps, thread-banner order and wall-clock, and in
    nothing else. The diagnostic exists because this defect was **documented in a comment**
    and survived anyway — the same lesson as the pressure-solver race, and the reason
    CLAUDE.md says a cross-reference is not a check.

40. **Where the 5.5 h day comes from, and what it is worth (documentation).**

    `omega = 3.17e-4 rad/s` has been in `param.py` since item 8 with a two-line note calling it
    an assumption inside a 4–6 h literature range. The argument is worth more than that note
    carries, and weaker in one respect the note omits. Nothing here changes a parameter; this
    records the derivation so the number can be argued with.

    **The argument is angular-momentum conservation of the Earth–Moon system.** Tidal friction
    moves angular momentum from Earth's spin into the Moon's orbit — the Moon recedes 3.8 cm/yr
    by laser ranging, the day lengthens ~2 ms/century. The total is very nearly conserved, so
    the relation inverts: a closer Moon is a faster Earth.

    ```
    I_earth      = 8.016e+37 kg m^2      (0.3307 M R^2)
    L_spin(now)  = 5.845e+33   (17 %)
    L_orbit(now) = 2.856e+34   (83 %)
    L_total      = 3.441e+34 kg m^2/s
    ```

    Placing the Moon at a given distance and giving Earth the remainder:

    ```
      a [R_E]      day [h]
          3.0         4.99
          5.0         5.34
          5.9         5.49
         10.0         6.14
         20.0         7.79
         60.3        23.90     <- today, a check on the arithmetic
    ```

    **A 5.5 h day is the Moon at 5.95 R_E**, giving ω = 3.173e-4 rad/s = 4.35× modern, which is
    the configured value.

    **Why the number survives the Hadean being unconstrained.** The useful property is that
    table's flatness. Moving the Moon from 3 to 10 R_E — a factor 3.3 in distance, and an
    enormous span of tidal-evolution time — moves the day only from 5.0 to 6.1 h, because once
    the Moon is close Earth's spin holds most of the angular momentum and the split changes
    slowly. Dating *when* the Moon was where is therefore not required to bracket the day
    length; it is enough that the Moon was within a few tens of Earth radii, which it was for
    the whole Hadean. That is why 4–6 h is robust and 5.5 h a defensible midpoint — a better
    argument than "estimates range 4–6 h" was carrying.

    **What would break it, in order of seriousness.**

    - **The high-angular-momentum impact scenarios.** Ćuk & Stewart (2012) and Canup (2012)
      proposed Moon-forming impacts leaving Earth spinning at ~2–3 h with substantially more
      angular momentum than the system now has, the excess later shed through the evection
      resonance with the Sun. If that is right `L_total` was *not* conserved across the Hadean
      and the inversion above is invalid — in the direction of a **shorter** day.
    - **Nothing measures this.** Tidal rhythmites and cyclostratigraphy reach ~2.5 Ga
      (~18–19 h); at 4.4 Ga the record is essentially Jack Hills zircons. The 4–6 h range is
      theory, not observation.
    - **The timescale problem cuts the other way.** Extrapolating today's recession rate
      linearly backwards puts the Moon at Earth's surface ~1.5 Ga, far too recent, so present
      dissipation is anomalously high (ocean tides near a basin resonance) and past averages
      were lower. That *supports* a still-short day at 4.4 Ga.

    **Two torques specific to this model that the bookkeeping omits — the part the generic
    4–6 h range hides.** A 250 bar atmosphere is ~250× Earth's air mass, so thermally driven
    atmospheric tides are not a negligible torque; on Venus they drive the rotation the *other*
    way. And tidal dissipation in a molten or quenching surface is not that of solid Earth or
    of an ocean. Both act during precisely ATHAD's epoch and neither is in the angular-momentum
    inversion, so `omega`'s justification is weaker **here** than the literature range implies.
    Neither is quantified — flagging them is not bounding them.

    **What the uncertainty costs the model.** Rotation enters as `f = 2Ω sinφ`, so it sets the
    cell geometry. Item 31's table is confirmed independently (Earth `Ro_T` = 0.0596, Held–Hou
    edge 18.1°; ATHAD 0.00474 and 5.09°). Propagating the range:

    ```
     day[h]      Ro_T   HH edge   Rhines band   n_cells/hemi
        4.0   0.00251     3.70d         17.0d            5.3
        4.5   0.00317     4.17d         18.1d            5.0
        5.0   0.00392     4.63d         19.1d            4.7
        5.5   0.00474     5.09d         20.0d            4.5
        6.0   0.00564     5.56d         20.9d            4.3
    ```

    The two estimates degrade differently, as item 36 noted: the Held–Hou edge goes as 1/Ω and
    the Rhines band only as 1/√Ω. Over 4–6 h the Hadley edge spans 3.7–5.6° and the
    extratropical bands 17–21°.

    **This is a caveat on item 38 that item 38 does not carry.** Its Rhines argument — 20° bands
    against ~22° — is one of the three grounds for `n_cells_hemisphere = 5`, and that ground
    moves with `omega`. At the nominal 5.5 h the Rhines estimate is **4.5 cells per hemisphere**,
    which selects 4 or 5 equally; it takes the short end of the rotation range to make 5 the
    clear answer. So n = 5 is *inside* the rotation uncertainty but is not *selected* by it, and
    the layout is structurally reasonable rather than derived. Item 38's other two grounds — the
    removal of the 40°-wide-cell artefact, and Ψ maxima still on their prescribed cores at 400
    iterations — are untouched by this and remain the stronger pair.

41. **The `zeta` scan: the OLR does not respond to the grid either, and the −1.26 W/m²
    residual is a `zeta = 3.0` artefact.**

    Item 39 gave three converging arguments for cutting the stretch. This tested them. Five
    runs at `im = 41`, shell held at 300 km, `dt_visc` fixed at 1e-4, 200 iterations,
    `zeta` ∈ {3.0, 2.0, 1.5, 1.0, 0.5}. All five completed at 200, rc = 0, no NaN.

    **Two design constraints, both forced and both worth stating.** Holding the shell at
    300 km forces `L_atm = 300000/(exp(zeta)−1)`, which runs 15.7 km → 462 km, a **29×**
    swing. `L_atm` is the length normalisation for the Held–Suarez relaxation, the Rayleigh
    drag, `coeff_MC_*`, `coeff_S`, `coeff_L`, `nue_max` and `coeff_u_p` — `RHS_Atm_Turb.cpp:441`
    already says so and calls them "STILL ON L_atm, and therefore still 40x too weak".
    `force_nd` is exempt, having been moved to `metricShellLength()`, which returns 300 km for
    every run here. So **the dynamical half of the scan is confounded and the radiative half is
    not**: under the default prescribed `t`, `densities()` re-imposes the adiabat over the grid
    every iteration, so the OLR depends on the grid and not on those coefficients — provided
    moist physics is off, since `coeff_MC_*` and `coeff_L` feed the clouds that feed the
    radiation. Hence `moist_phys_start_iter > nm`. This is item 39's defect in another place:
    `L_atm`, an exponential-stretch amplitude, used as though it were a grid length.

    ```
    run    zeta  spread   albedo   absSW      OLR   imbalance   t_skin   lid eps
    z3      3.0  11.77x   0.4986  271.11   272.37     -1.25     262.96    0.0069
    z2      2.0   4.55x   0.4975  271.37   271.64     -0.27     263.02    0.0434
    z1p5    1.5   2.83x   0.4974  271.38   271.38     -0.00     263.02    0.0000
    z1      1.0   1.76x   0.4961  271.74   271.74      0.00     263.11    0.0000
    z0p5    0.5   1.09x   0.4945  272.13   272.13      0.00     263.20    0.0002
    ```

    The measured metric spreads reproduce item 39's analytic table digit for digit, which is an
    independent check on `checkRadialMetric()` and on the offline reconstruction both.

    **The scan cannot answer its primary question.** `OLR = absorbed SW + geothermal =
    σT_lid⁴` exactly, to the last printed digit, in the three low-`zeta` runs. Total OLR
    spread over a 6× change in stretch is **0.99 W/m², 0.36 %** — and that is not a radiative
    response to resolution. The chain runs the other way: albedo falls 0.4986 → 0.4945,
    absorbed SW rises 271.11 → 272.13, `t_skin` rises 262.96 → 263.20, and the OLR tracks
    absorbed SW penny for penny. **`zeta` → cloud deck → albedo → absorbed SW → `t_skin` → OLR**,
    with the photosphere playing no part. This is item 25's fixed point, not the opaque-lid
    failure: lid emissivities are ~0 in every run and the diagnostic's opaque-lid flag never
    fired. **So item 39's photosphere argument is not confirmed, and could not have been** —
    nothing moves the OLR while `t_skin` pins it. That `dtau@tau=1` falls 54.8 → 2.8 across
    this range remains a computed property of the grid, not a measured improvement in output.

    **This is the fifth non-response.** κ 64× → 0.10 % (item 29), circulation 500× → 0.03 %
    (item 27), `im` analytically (item 39), and now stretch 6× → 0.36 %, albedo-driven.

    **The −1.26 W/m² residual is a grid artefact.** CLAUDE.md calls it the least trustworthy
    number in the file; it is worse than that. At `zeta = 3.0` it reproduces exactly (−1.25,
    and still falling at 200 iterations). At `zeta ≤ 1.5` the imbalance closes to **0.00
    exactly**. It is not physics — it is the shipped grid failing to reach its own fixed point
    within 200 iterations while the finer-topped grids land on it. The convergence behaviour
    agrees: z3's OLR is still falling monotonically at 200 (276.50 → 272.37, −0.52 per 20
    iterations), z2 falling slowly, z1p5 flat, z1 and z0p5 slowly *rising*. **This is an
    argument for cutting `zeta` that does not depend on the photosphere at all**, and it
    confirms item 25's diagnosis by an independent route.

    **`initCloudIce`'s grid dependence is now measured, and it is the only live path to the
    OLR in this model.** Max cloud water jumps 37.3 → 49.9 g/kg between `zeta` 3.0 and 2.0.
    Moist physics was off, so this is purely the initialisation deck, laid by a parabola keyed
    to `p_crit = 1000` hPa and `p_mid = 550` — Earth surface pressures — landing on different
    levels for different grids. This README has listed that as an open item with the note that
    it "has a path to the OLR". It does, and in this scan it is the only thing that used it.

    **What this does not test.** The 13×-timestep claim: `dt_visc` was held fixed, and the
    clean completion of `zeta = 0.5` is weak evidence at best, since that run also carried 29×
    shifted drag, relaxation and viscosity. The diffusive CFL remains untested and remains
    blocked behind moving those coefficients onto `metricShellLength()`.

    **The headline is not "cut `zeta`".** It is that **`t_skin` has to be broken before any
    grid or opacity question can be measured at all** — it now blocks item 29's κ question,
    item 39's photosphere question and this one. Three lines of work behind one fixed point.


42. **Instruments for the things the model could not see — and the first one corrected a
    number this file had been quoting from me rather than from the model.**

    Items 39 and 41 both rested on the photosphere: the level where the cumulative optical
    depth reaches 1, which is where the emission to space originates. **The model computed it
    nowhere.** `epsilon` is a 3D array of *per-layer* emissivity, `tau` was a local scalar
    overwritten every cell, and the cumulative value existed only in an offline Python column
    of mine. Nor could it be recovered after the fact: `epsilon = 1 − exp(−tau)` saturates to
    exactly 1.0 in double precision at `tau ≳ 37` and the deep column carries `tau ~ 9e4`
    **per layer**, so the information is destroyed at write time and has to be accumulated
    where `tau` is still in scope.

    Five diagnostics, all print-only:

    | field | what it is for |
    |---|---|
    | `tau_above` | cumulative LW optical depth from the lid down; the photosphere is where it crosses 1 |
    | `tau_layer` | per-layer `dtau` — the model's own measure of how well it resolves its photosphere (item 39) |
    | `ubud_*` (6) | the **radial** momentum budget, which did not exist while θ and φ both had one |
    | `N2` | Brunt–Väisälä frequency squared — measures "neutrally stratified by construction" |
    | `Psi` | the meridional streamfunction as a field, not just a CSV and a scalar |

    All are in the three ParaView writers and in Results min/max. `printPlanetaryBalance` now
    prints the `tau_above = 1` height, the temperature there, and the fraction of columns
    radiating from within 1 K of `t_skin`, flagged past 50 %.

    **The measurement that matters: the photosphere is not where this file said it was.**

    ```
                          this file (my offline column)     the model
    photosphere height    218 km  (level 36)                239.4 km
    temperature there     402.6 K                           325.2 K
    ```

    **21 km higher and 78 K colder.** The offline column was dry and cloud-free; the runs have
    moist physics on, and cloud liquid and ice add `k_liq·LWP + k_ice·IWP` to the layer optical
    depth, pushing the crossing upward. Every photosphere figure quoted before this item came
    from the reconstruction, not the model, and should be read with that 21 km in mind — item
    39's "the entire photosphere is one grid cell of `dtau = 55`" included.

    **`initCloudIce`'s H_crit is now on a fraction of surface pressure, and the fix is a
    correctness cleanup rather than the lever it was billed as.** The inherited parabola used
    `x = p/1000` hPa with its minimum at 550 hPa — Earth's surface pressure and Earth's
    mid-troposphere. At 250 bar levels 0–35 have `p > 1000` hPa, so `x > 1`, the parabola goes
    negative and only the clamp saves `H_crit` → 1.0; then at level 36 it hits its minimum of
    0.8010. Rescaled, the minimum moves to ~36 km — deep in the supercritical column where
    nothing can condense — and `H_crit ≥ 0.9926` everywhere condensation is possible, so the
    Earth profile goes inert here **by derivation** rather than by fiat. Exactly equivalent on
    Earth, so it ports upstream. `ATM_HCRIT_ABS=1` restores the old form.

    Two 100-iteration arms, moist physics from iteration 0, 24 threads, one binary:

    ```
                       fixed      inherited      delta
    albedo             0.4987     0.4987         0
    absorbed SW        121.07     121.07         0
    OLR                299.04     298.49         +0.55 W/m2  (0.18 %)
    max cloud water    13.556493  13.557265      -0.006 %
    T at photosphere   325.20 K   322.56 K       +2.64 K
    ```

    **The predicted mechanism failed, twice over, and both failures were caught by
    instruments added after the prediction.** The chain proposed was: higher threshold at the
    photosphere → less cloud → lower albedo → more absorbed SW → higher `t_skin` → higher OLR.
    The albedo did not move **at all**, because reflectivity saturates on the *presence* of
    condensate and not its amount — so that chain is structurally impossible in the shipped
    configuration, whatever the cloud amount does. And the perturbation was mis-sized because
    of the 21 km above: at the real photosphere (~239 km, levels 37–38) the inherited `H_crit`
    was **0.917–0.985**, a 2–8 % effect, not the 20 % the parabola's minimum at 218 km implied.
    An order of magnitude of the ~200× shortfall is accounted for by my own misplacement of
    the level. The sign held; nothing else did.

    This also leaves item 41 owing an explanation: if cloud *amount* cannot move the albedo,
    then whatever moved it across the `zeta` scan (0.4986 → 0.4945) was cloud **coverage** —
    which levels hold any condensate at all — and that attribution is unverified.

    **An unexpected result, and the first testable prediction this session has produced that
    is not blocked behind `t_skin` — because it is a measurement of `t_skin`.** At iteration
    100 only **17.1 %** of columns radiate from within 1 K of `t_skin`; the flag needs 50 % and
    did not fire. So the pinning is not yet dominant at 100 iterations. But these runs are
    28 W/m² from balance, and item 41 showed that convergence drives the OLR down onto
    `σT_lid⁴` exactly. If item 25's mechanism is right, the skin fraction must **rise toward
    100 % as the imbalance closes** — the radiating level migrating into the isothermal top is
    precisely what item 25 claims. That is now falsifiable rather than inferred. See item 43.

    **The fields validate against independent references, which is why they were checked before
    being trusted.** At the first populated checkpoint: `tau_above` max **2 196 033** at the
    surface against the offline column's 2.28e6 (3.6 %); `tau_layer` max **111 181** at 12.9 km,
    peaking above the surface as the offline `dtau` profile also does; `ubud_pgf` **+1.19 / −1.26**
    nd at 18°S/94.8 km and 7°N/72.5 km; and `Psi` max **1.0651e14 kg/s** against the existing CSV
    writer's own `Psi_max = 106510.84 (1e9 kg/s)` — an **exact** match, which is the check that the
    replication across `k` is right.

    **`N2` is the one that earns its place immediately: it confirms invariant 4 by measurement
    for the first time.** It is 0.000000 through the column and reaches **2.74e-4 s⁻² only at
    277 km** — neutral exactly where the adiabat is imposed, stably stratified only in the
    isothermal skin. An isothermal layer gives `N² = g²/(cp·T) ≈ 2e-4` at 263 K, so the magnitude
    is right as well as the shape. "Neutrally stratified by construction" is no longer an
    assertion.

    **`ubud_*`'s first result measures two things this file had only asserted.** All six terms
    from a single capture, max |value| in nondimensional tendency:

    ```
    ubud_pgf    1.153554     <- dominates the next term by 136x
    ubud_cor    0.000000     <- identically zero
    ubud_advv   0.000840
    ubud_advh   0.008502
    ubud_diff   0.000071
    ubud_buoy   0.000001     <- 1.15e6 times smaller than the pressure gradient
    ```

    `ubud_cor` is **exactly** zero in the same capture in which `ubud_pgf` is 1.15, which confirms
    `coriolis_nontraditional()` is genuinely off rather than merely documented as off — the check
    CLAUDE.md asks for ("read the switches, not just the RHS") now runs itself.

    And `ubud_buoy` at 1e-6 against a 1.15 pressure gradient is **item 34 measured**: the
    buoyancy body force carries an extra `*dt`, so with `dt = 1e-4` it enters ~1e-4 too weak and
    is effectively absent from the radial balance. That makes item 28's claim quantitative — with
    the default switches `rhs_u` really is the radial pressure gradient *and nothing else*: no
    non-traditional Coriolis, no curvature, and a buoyancy term six orders of magnitude down.

    **Two ordering caveats in the new fields, from the run loop rather than the wiring** — both
    now confirmed empirically: the first checkpoint reports `ubud_*` and `Psi` as identically
    zero and the second reports the values above.
    `ubud_*` in the **vtk lags one checkpoint interval**: the vtk is written at
    `cAtmosphereModel.cpp:1632` and `ubudget_capture` is set at 1671, after it, so iteration
    20's file carries the split captured during iteration 0's RK4. `vbud_*`/`wbud_*` escape this
    because their CSV is written at 1748, after capture. And `Psi`'s **Results min/max is one
    checkpoint stale** (zero on the first), because `print_min_max_atm()` runs at 1630, one line
    before the fill at 1631 — the vtk value itself is current. Both are fixed by moving the
    capture flags ahead of the checkpoint block; not done, to avoid reordering the run loop
    while a measurement was being taken on it.

43. **The prediction of item 42 is confirmed in direction and unsettled in magnitude — and the
    imbalance turns out not to be an independent quantity at all.**

    Item 42 predicted that if item 25's mechanism is right, the fraction of columns radiating
    from within 1 K of `t_skin` must rise as the imbalance closes. Tested on the shipped
    configuration, moist physics from iteration 0, 24 threads, **stopped at iteration 200 once
    the trend could be judged** rather than run to 400:

    ```
    iter   z_ph    T_ph  skin%    imbal      OLR   t_skin
      20   237.4  368.41    1.1    -67.78   338.85   262.95
      40   237.4  368.00    1.1    -52.70   323.77   262.95
      60   237.5  358.01    1.1    -42.54   313.61   262.95
      80   238.2  339.10    7.2    -33.98   305.05   262.95
     100   239.4  324.17   17.1    -27.63   298.70   262.95
     120   240.9  305.03   20.4    -22.98   294.06   262.95
     140   242.2  294.99   20.4    -18.85   289.92   262.95
     160   243.4  282.04   30.3    -15.49   286.56   262.95
     180   244.7  273.24   37.0    -12.79   283.86   262.95
     200   245.9  270.36   39.2    -10.51   281.59   262.95
    ```

    **Confirmed, in direction.** The skin fraction rises 1.1 → 39.2 % while the imbalance closes
    −67.78 → −10.51 and the photosphere descends from 368 K to 270 K against a `t_skin` of
    262.95. The radiating level migrating into the isothermal top is exactly item 25's claim, and
    it is now observed rather than inferred. The stair-stepping (plateaus at 1.1 and at 20.4) is
    what a threshold diagnostic does on a hemispherically symmetric field: the surface
    temperature is a smooth parabola in latitude, so whole bands cross "within 1 K" together.

    **Not settled, and this file's own rule applies.** 39.2 % is not near 100 %, and both rates
    are decelerating — the photosphere cooled 8.8 K over iterations 160–180 and only 2.9 K over
    180–200, and the imbalance decrements decay 15.1, 10.2, 8.6, 6.4, 4.7, 4.1, 3.4, 2.7, 2.3.
    **Where the skin fraction saturates is unknown and is not extrapolated here.** A monotone
    trend is not a limit; this file has been caught on that twice on radiation numbers already.

    **The sharper result, which was not the thing being tested.** `t_skin` is **262.95 K at
    every one of the ten diagnostics** — constant to five significant figures across a run in
    which the OLR fell by 57 W/m². Since the absorbed flux is the `t_skin` fixed point,
    σT_skin⁴ = 271.2 W/m², and the imbalance is therefore *identically* the OLR's distance from
    that constant: 281.59 − 271.2 = 10.4 against a reported −10.51. **"The energy balance
    closes" and "the OLR arrives at a number that never moved" are the same event, not two
    agreeing measurements.** Item 25 said the converged imbalance was the least trustworthy
    number in the file; this shows it is not really a second number at all. Any future run
    reporting a small imbalance is reporting how close the OLR has drifted to σT_skin⁴, and
    nothing else, until the fixed point is broken.

    **One incidental observation.** The photosphere *rises* 237.4 → 245.9 km, 8.5 km, while
    cooling 98 K. It is not descending through a fixed profile; the profile is cooling underneath
    it and the optical depth is redistributing. Worth remembering when reading item 39's
    resolution argument, which treated the crossing height as a property of the grid.

    **Method note.** The run was stopped at 200 deliberately. It also produced no usable vtk —
    `checkpoint` was set to 400 to cut I/O, so only the iteration-0 file was written, where
    `tau_above` is identically zero because radiation has not yet run. A separate short run with
    `checkpoint = 20` was needed to get the new fields into ParaView at a meaningful state.


44. **The surface drag is Earth's, twice over: its strength was fitted to a jet off Chile and
    its depth is a cell count, which makes it 30× deeper here.**

    **CLOSED by item 49: both defects are real and neither matters.** The scan ran, and a
    correctly-scaled drag accounts for 0.034 % of the Ψ decay (0.31 % at ten times the rate).
    The paragraph below arguing this is "a first-order suspect" was right that it had never
    been varied and wrong that varying it would explain anything — and the scan as first
    written could not have shown either way, because the drag entered as dt². Read item 49
    before using anything here.

    Neither is a config parameter. Both are `constexpr` in `RHS_Atm_Turb.cpp`, and both carry
    comments justifying them by Earth's geography — the pattern CLAUDE.md catalogues, in a file
    nobody has re-read since the fork.

    **The strength.** `rayleigh_kf = 1.0/86400.0`, with the comment: *"baseline 1/day gave ~34 m/s
    eastward w off W-coast S-America; 10× cut the surface to ~28 m/s and 20× gained nothing (the
    jet max sits at ~1 km, above the drag layer), so 10× is the settled strength."* ATHAD has no
    South America. Invariant 1 guarantees it: `h ≡ 0`, `is_land()` false everywhere.

    **The depth, and this is the part that bites.** `drag_n_layers = 5.0` is a count of **air
    cells**, not a length, so its physical depth is whatever the grid makes it:

    ```
    ATOM_Precipitation (Earth, L_atm = 400 m, zeta = 3.715)   5 cells =   236 m
    ATHAD              (L_atm = 15719 m, zeta = 3.0)          5 cells =  7152 m
                                                              ratio   =  30.3x
    ```

    The comment calls 5 "the physical 5" and defends it against a 10-cell test by noting the
    coastal jet max at ~1 km sits *above* the drag layer — true when the layer is 236 m. Here the
    same constant applies a momentum sink through **7.2 km** of the column. It is the same defect
    shape as `init_tropopause_layers`' `round(h / L_atm)`: a grid index treated as a physical
    length, which is only a length on the grid it was written for.

    **It is live, and it carries two other known defects.**
    `surf_drag = (rayleigh_kf * L_atm / u_0 * dt) * drag_profile` multiplies `rhs_v` and `rhs_w`
    (not `rhs_u` — deliberately). The `* dt` inside it is item 34's extra factor, so the drag
    enters as dt² like the buoyancy term; and it divides by `L_atm`, so it is on item 41's list of
    coefficients that a shell-preserving `zeta` change moves by 29×. Three separate recorded
    problems meet in one expression.

    **Why this matters now rather than as bookkeeping.** The open question about the dynamics is
    the tropical cell decaying with no turnover (items 28, 31, 37), and drag is the obvious sink.
    A momentum sink 30× deeper than intended, entering at dt² and riding on `L_atm`, is a
    first-order suspect that has never been varied. **Unlike the radiative questions this one is
    not blocked behind `t_skin`** — Ψ responds to drag directly, and `ubud_*`/`vbud_*` can now
    attribute it. Nothing is changed here; both constants should become parameters first, so the
    scan is possible at all.


45. **The `t_skin` fixed point does break under prognostic temperature — and what escapes it
    radiates away eleven times what it absorbs.**

    Items 41–43 established that every radiative measurement was pinned: `OLR = σT_skin⁴ =
    absorbed`. This is the run that tests whether the pinning is a property of the *model* or of
    the *prescribed profile*. `ATM_PROGNOSTIC_T=1 ATM_RAD_DIRECT=1`, shipped config otherwise,
    moist physics from iteration 0, 24 threads, 400 iterations, `checkpoint = 20`. Completed
    rc = 0, no NaN, 21 vtk checkpoints.

    ```
    iter         OLR       imbal     z_ph     T_ph   skin%   t_skin
      20    21239.70   -20968.72    240.4   459.83     0.0   262.92
     100     8700.29    -8429.29    249.6   397.08     0.0   262.93
     200     6254.32    -5983.31    252.4   385.70     0.0   262.93
     300     3930.43    -3659.41    254.0   374.55     0.0   262.93
     400     2926.65    -2655.63    254.8   368.05     0.0   262.93
    ```

    **The pinning is a property of the prescribed profile, not of the model.** `skin%` is
    **0.0 at every one of twenty diagnostics** — not one column radiates from within 1 K of
    `t_skin`, and the photosphere sits at 368 K against a `t_skin` of 263. With the adiabat no
    longer re-imposed there is no isothermal top for the radiating level to migrate into, so the
    mechanism item 43 measured cannot operate. That answers the question items 29/39/41 were
    blocked on: the fixed point *can* be escaped.

    **But what escapes it is not a steady state.** The OLR is **2927 W/m², 10.8× the 271 W/m²
    absorbed**, and still falling at iteration 400, with the first non-monotonicity in the series
    at 360 → 380 (3217 → 3244).

    **This extends item 30 and vindicates its own warning.** Item 30 measured 15 764 → 8 542 →
    6 211 over 200 iterations and said explicitly that a monotone trend is not a limit. This run
    reproduces it closely (8 700 at 100, 6 254 at 200) and then falls a further **53 %** by 400.
    So 6 211 was nowhere near a floor — and **2 927 is not claimed to be one either.**

    **Item 43's result survives into a regime ten times different.** `t_skin` is 262.93 K at all
    twenty diagnostics here too, because it is a fixed point of σT⁴ = absorbed SW + geothermal and
    depends on the albedo and insolation, not on the column. So the imbalance is *again* exactly
    the OLR's distance from a constant: 271.02 − 2926.65 = −2655.63, matching the printed value to
    the last digit. **The imbalance is never an independent measurement, in either configuration.**

    **The open question this creates.** The effective emission temperature at 400 iterations is
    (2927/σ)^¼ = **476 K**, well above the **368 K** at `tau_above = 1`. Emission is therefore
    coming from deeper and hotter than the photosphere, which should not happen in a grey
    two-stream column unless the opacity is letting deep hot layers leak out. That wants
    understanding before more integration is spent — `tau_layer` is now the instrument for it.

    **Method notes.** A restart dump `atm_restart_0Ma_300.bin` was written **despite
    `restart_stride = 0`**, which is supposed to disable it; useful here, but the knob does not
    do what it says. The run is being continued from that checkpoint rather than recomputed.


46. **The prognostic OLR does stop falling — at 5.2× the absorbed flux, because the prescribed
    surface is an infinite reservoir.**

    Item 45 left the endpoint unknown and refused to extrapolate. Continued from the
    iteration-300 restart to **901** with the same knobs (`ATM_PROGNOSTIC_T=1
    ATM_RAD_DIRECT=1`, moist from 0, 24 threads, `checkpoint = 20`). Completed rc = 0, no NaN.

    ```
    iter      OLR    T_ph    z_ph        iter      OLR    T_ph    z_ph
     401  2900.29  367.92   254.6         661  1332.21  357.89   255.5   <- minimum
     501  2190.19  364.16   255.1         701  1422.79  358.49   255.5
     601  1673.70  361.07   255.3         801  1360.07  359.42   255.7
     641  1477.96  359.47   255.3         901  1429.75  360.23   255.8
    ```

    **There is a genuine turning point near iteration 660–680.** The OLR bottoms out around
    1300 and then flattens at **~1400 W/m², scattering ±6 % with no further trend** — mean 1363
    over 661–781 against 1425 over 801–901, i.e. flat to slightly *rising*. `T_ph` turns at the
    same place: a minimum of **357.89 K at 661**, back to 360.23 by 901. `z_ph` saturates at
    255.7–255.8 km. So the answer to item 45's question is **no** — it does not approach
    σT_skin⁴ = 271, and it does not keep falling either.

    **Why it lands there, and this is the part that matters.** 1400 W/m² is **5.2× the 271 W/m²
    absorbed**, with a persistent imbalance of ≈ −1130 W/m². That is not an energy-conserving
    steady state, and it cannot be: **`t_surf_equator` = 1500 K is prescribed**, so the surface
    is an infinite reservoir. The column drains heat from that boundary and radiates it away, and
    the steady state is set by how much flux the prescribed surface can drive through the column,
    not by a top-of-atmosphere balance. The TOA imbalance does not close because the boundary
    condition supplies it indefinitely.

    So this run converged to a real fixed point *of the model as posed*. CLAUDE.md's standing
    caveat — "the surface temperature is prescribed, not solved; every result is conditional on
    it" — is not a footnote to prognostic mode, it **is** prognostic mode's answer. Freeing the
    column from the adiabat (item 45) removes one prescribed profile and leaves the model pinned
    to another, one layer down. **A genuine radiative equilibrium needs the surface temperature
    solved, not just the column freed.**

    **Item 45's "first non-monotonicity" was over-read, and this corrects it.** I flagged
    360 → 380 (3217 → 3244) as if a single reversal might signal an approaching floor. It did
    not. The series carries ±5–6 % scatter throughout, with reversals at 361, 421, 481, 581, 701
    and after. The real floor arrived **280 iterations later** and announced itself as a turning
    point in `T_ph`, not as a reversal in a noisy OLR. One reversal in a scattered series is not
    a signal — the same lesson as "a monotone trend is not a limit", in the opposite direction.

    **Item 45's open question is shrinking on its own.** The gap between the effective emission
    temperature and the photosphere narrows from **476 K vs 368 K (108 K)** at iteration 400 to
    **398 K vs 360 K (38 K)** at 901, consistent with the column approaching a self-consistent
    profile rather than with a leak. Worth re-checking with `tau_layer` but no longer alarming.

    **Restart fidelity: attempted, INCONCLUSIVE, and the design is why.** Recomputed iterations
    against the original differ by ±5–6 %, but the series' own scatter between diagnostics is
    ±5 %, so the comparison cannot discriminate a faithful restart from a mildly unfaithful one.
    Two further confounds: the restart shifts the diagnostic phase by one iteration (321 rather
    than 320), so an equality check matches nothing and any comparison needs interpolation across
    a 20-iteration gap in a noisy series; and the binary differed between the two runs (12:20
    against 12:25 build — items 2 and 4, physics-neutral by construction but not proven).
    What *is* established: `t_skin` is 262.93 in both, and the 29 serialized arrays are the full
    prognostic set including the RK4 `n` copies. Evidence favours fidelity; it is not verified.
    The decisive test is a from-scratch rerun to 401 on the current binary, which would settle
    the rebuild question at the same time. **Sequence runs after builds, not before.**


47. **The reservoir is 250 bar of gas, not a boundary condition — and items 45–46 are corrected
    by that. The model cannot reach radiative equilibrium by integration.**

    Item 46 concluded that the prognostic OLR plateaus at 5.2× the absorbed flux "because
    `t_surf_equator` = 1500 K is prescribed, so the surface is an infinite reservoir." **That
    mechanism is wrong.** The surface is not prescribed during the run.

    **The surface already floats.** `MultiLayerRadiation.h:452` solves a linearised surface
    balance — absorbed SW + geothermal + downwelling LW − σT_s⁴ + sensible exchange with the
    first air level — and writes `m.t.x[i_mount][j][k]` directly, with `i_mount` = 0. Under
    `ATM_PROGNOSTIC_T=1` nothing overwrites it, so it was floating throughout items 45 and 46.
    It settles at **1498.7 K** (min 1449.0, mean 1481.9) because the column immediately above
    radiates ~280 kW/m² downward onto it: `σT_surf⁴ = 280 258 W/m²` against an OLR of ~1190, a
    suppression of **×236**. Surface and first air layer are radiatively locked, so floating the
    surface changes nothing while the air above it is hot.

    **The decisive number, which no item before this one computed.**

    ```
    column heat capacity   p/g·cp            = 5.20e9 J/m²/K
    at the measured -1130 W/m² net           = 6.85 K per YEAR
                                             -> 100 K of cooling takes 14.6 years

    simulated time, 901 iterations at dt=1e-4:
        time unit L_atm/u_0     (15.7 km)    =  177 s   (0.05 h)
        time unit shell/u_0     (300 km)     = 3375 s   (0.94 h)
    ```

    **Under an hour of physical time, against a multi-decade radiative timescale — 5 to 7 orders
    of magnitude short.** The surface cooled 1500 → 1498.7 K, which is exactly what 6.85 K/yr
    predicts. There was never a paradox in item 46's persistent −1130 W/m²: the column *is*
    draining at precisely that rate, and will be for decades.

    **So item 46's plateau is not a converged state.** 1400 W/m² is the quasi-steady flux from a
    reservoir that has cooled by ~1 K. The turning point near iteration 660 is a fast local
    adjustment of the thin upper column against a deep reservoir that has barely moved, not
    thermal equilibrium. **Item 45's refusal to extrapolate was right; item 46's claim to have
    found the endpoint was not.**

    **A new defect: no run in this README has a stated physical duration.** The two candidate
    time units differ by **19×** because the non-dimensionalisation divides by "L/u_0" and it is
    ambiguous which L — `L_atm` (the stretch amplitude) or `metricShellLength` (the shell). That
    is item 41's exact defect appearing in the *time* coordinate, where nobody had looked. Until
    it is settled, every "N iterations" in this file is a count and not a duration.

    **This reframes invariant 3.** That invariant treats the prescribed adiabat plus the `t_skin`
    fixed point as a deficiency to be removed — "radiation must *set* the profile, not nudge it
    toward a prescribed one." But integration cannot set it: reaching equilibrium needs of order
    10⁸ iterations. **The prescribed profile is not laziness, it is the only equilibrium this
    model can currently obtain.** The honest alternatives are to operator-split the thermal
    equation and relax temperature implicitly at a far larger effective step, to reduce the
    column heat capacity during spin-up and restore it afterwards, or to solve the
    radiative-convective column directly — which is what the prescribed profile already is.

    **And a Neumann condition on `p_stat` at the surface would not help either.** The surface
    pressure of a hydrostatic atmosphere is not a free boundary value; it *is* the weight of the
    air above, `p_s = g·M/A`. A zero-gradient condition asserts `dp/dz = 0` exactly where the
    true gradient is largest in the domain (`−ρg` = −421 Pa/m), contradicting the hydrostatic
    integration `densities()` performs immediately afterwards. There is no pressure reservoir to
    release: the thermal inertia is `p_s/g·cp` = column MASS × cp, so letting `p_s` drift would
    create or destroy air rather than reduce the heat capacity. Item 19 already established this
    and pinned `p_prev = m.p_0` (`ThermoAtm.h:1467`), cutting the mass drift from −0.2128 % to
    −0.0011 %. A cooling column contracts at constant mass: `p_s` stays 250 bar and the shell
    gets shallower, which is already the behaviour.

48. **ATHAD has no functioning surface boundary-condition layer. Every surface BC is an
    identity.**

    Found while checking whether a different surface condition could be substituted. Three
    separate sites impose surface values, and all three are exact no-ops:

    | site | writes | reduces to |
    |---|---|---|
    | `BC_Atm.h:433` (`bcSolidGround`) | `t`, `radiation` | `t.x[0] = t.x[0]` |
    | `BC_Atm.h:1108` (`bcScalarSurfSur`) | `t`, `radiation`, `c`, `cloud`, `ice` | identity |
    | `ThermoAtm.h:1540` (`densities`) | `p_stat`, `r_dry`, `r_humid` | identity |

    Every one is written as "copy from the topography level" — `x.x[0][j][k] =
    x.x[i_topography[j][k]][j][k]` — and invariant 1 fixes `i_topography ≡ 0`, so each is
    `x[0] = x[0]`.

    **So the surface values are simply whatever the last physics routine happened to write**:
    `MultiLayerRadiation` for `t` (its surface energy balance, item 47), `densities()` for
    `p_stat`, `r_dry` and `r_humid`. Nothing *chooses* the surface condition; it is a side effect
    of call order. That is why "make the surface float" turned out to be already in effect
    without anyone having decided it.

    **This is the dead-land-branch pattern at structural scale.** CLAUDE.md catalogues it as
    individual Earth-only fallbacks that never run at home and are never tested — `clampAndFade`,
    the dilute Newton loop, `M_nonwater`. Here it is not one constant or one branch but an entire
    boundary-condition layer, dead for the same reason: it is keyed to topography, and on Earth
    it is only non-trivial *over* topography, so it is barely exercised there either.

    **Consequence for any future work on the surface.** There is currently no place to impose a
    surface condition even if one were wanted — a surface energy balance, a prescribed flux, a
    fixed temperature — because the sites that look like they do it do nothing. Any such
    condition has to be *written*, not enabled. Not repaired here, because which condition is
    correct depends on item 47's open question, and repairing the plumbing before deciding the
    physics would just make the wrong answer look deliberate.


49. **The drag scan ran, and it measured the reproducibility noise floor. Drag is eliminated
    as the cell-decay driver — but only after a defect was found that made the scan incapable
    of measuring anything.**

    **What the four arms returned.** 200 iterations, `rayleigh_kf` × 0.1 / 1 / 10 and
    `drag_n_layers` 5 → 1, one thread count throughout:

    | iter | `base` | `kf01` | `kf10` | `n1` |
    |---|---|---|---|---|
    | 20 | 106510.84 | 106510.84 | 106510.84 | 106510.84 |
    | 100 | 102879.15 | 102879.15 | 102879.15 | 102879.15 |
    | 200 | 98393.19 | 98393.19 | 98393.18 | 98393.18 |

    **Eight significant figures across a 100× spread in drag.** The obvious reading — drag is
    irrelevant — is not available, because the arms agree to within the ~1e-8 thread-order
    noise this file already documents. **A null at the noise floor is not a measurement.**

    **Why: the drag entered as dt², so the scan varied a term ~7 orders of magnitude too weak.**
    `surf_drag` carried its own `* dt` and RK4 multiplies every RHS by `dt` again
    (`RungeKutta_Atm_Turb.cpp:169,201,228,255`). This is the second half of item 34's pair, and
    the Held–Suarez block thirty lines above had made exactly this repair on 2026-07-04, with
    the note that the extra `*dt` made it *"enter as dt² = ~1e-4 too weak → inert"*. The two
    have been inconsistent inside one file ever since.

    ```
    surf_drag                = 2.2742e-06   (already carried *dt)
    RK4's own *dt -> per step = 2.2742e-10
    over 200 iterations       = 4.5e-08     <- AT the ~1e-8 reproducibility floor
    ```

    100× of 1e-10 is still below the floor, so **the four arms could not have differed.**

    **The term is connected, which was checked rather than assumed.** Forcing `kf` = 1.0
    (86400×) moved `Psi_max` 106510.84 → 106507.38 at 20 iterations. Structurally, `surf_drag`
    sits at top-level indentation with no enclosing conditional and no early `return` before
    it, and `turb_model = k_omega_SST` with a single RHS file, so the laminar branch at
    `RHS_Atm_Turb.cpp:740` is not in play. *This check exists because the same day's
    `waterVapourEvaporation` claim was retracted in ATHAD_COND for exactly the opposite
    error — liveness inferred from a grep instead of read from the control flow.*

    **The answer, on a real measurement.** Two 200-iteration arms emulating the fix by
    pre-multiplying `kf` by 1/dt, against the shipped baseline:

    | | `Psi_max` @ 200 | drag costs | share of the decay |
    |---|---|---|---|
    | shipped (dt², drag inert) | 98393.19 | — | — |
    | correctly scaled, shipped `kf` | 98390.43 | 2.76 (0.0028 %) | **0.034 %** |
    | correctly scaled, 10× `kf` | 98367.67 | 25.52 (0.0259 %) | **0.314 %** |

    The two effects stand in a ratio of 9.25, against 10 expected — linear in `kf`, as a
    Rayleigh damping must be, which is the check that the emulation is doing what it claims.
    Ψ falls 7.62 % over those 200 iterations, and **the correctly-scaled drag accounts for
    0.034 % of it. Even at ten times the shipped rate it accounts for 0.31 %.**

    **So drag is eliminated — and item 44's Earth calibration does not matter.** Both halves
    of it, the strength fitted to a jet off Chile and the depth that is 30× too deep because
    it is a cell count, are real defects and neither is capable of driving this decay.

    **The reason is a timescale, and it generalises past drag.** `rayleigh_kf` = 1/86400 s⁻¹,
    while 200 iterations is **39 s to 12.5 min** of physical time depending on which `L` the
    time unit uses (item 47). A friction with a one-day e-folding time cannot do anything to a
    circulation observed for twelve minutes. **Anything with a timescale longer than minutes is
    disqualified as an explanation of a decay that completes inside 200 iterations** — which
    removes drag, radiative relaxation and surface exchange together, and points the remaining
    search at the initialisation transient and geostrophic adjustment (items 26–28), which act
    on the advective time the model actually integrates.

    **Fixed in `3e2d78f`, and ported to ATHAD_COND in its `a609d2b`.** The repair is safe
    against item 34's caveats for a reason specific to this half of the pair: drag is a
    *damping* term, so a larger coefficient is stabilising, and the "336× drove a polar
    vertical runaway" caveat recorded there is about buoyancy, a body force, which is **not**
    touched. It also removes a `dt` pathology rather than adding one — under dt² the effective
    drag scaled as dt², so item 24's `dt_visc` 4e-5 → 1e-4 had silently changed it by 6.25×,
    and ATHAD and ATHAD_COND differed in real drag strength by that same 6.25× purely through
    timesteps neither chose as a physics setting.

    **Verified in ATHAD_COND to the printed precision** (its numbers, 20 iterations): pre-fix
    42212.38, pre-fix at `kf` = 1.0 42212.33, post-fix at the shipped `kf` **42212.37** against
    a linear-scaling prediction of 42212.366, with the coefficient rising by exactly
    25000 = 1/`dt`. `make test` passes there, 0 failures.

    **THIS CHANGES RESULTS. Every Ψ, KE and wind number in this file predates the fix.**

    **What the scan does not settle.** It measured Ψ, and Ψ is built from `v` alone — item 28's
    warning that a radial failure is invisible to it still stands. `vbud_*`/`wbud_*` are written
    to CSV and were not analysed here; the attribution of *which* term does drive the decay is
    still open, and is now the question worth spending runs on.


50. **The buoyancy body force is item 34's other half, and it is wrong in three ways at once.
    A knob corrects all three, and it is shipped OFF — the arithmetic says so.**

    Item 49 fixed the Rayleigh drag's extra `*dt`. The buoyancy body force carries the same
    defect and was left alone, deliberately, and the reason is worth stating because it is the
    opposite of the drag's:

    ```
    current      g*dt/u_0   = 1.2263e-04
    consistent   g*L/u_0^2  = 2409.4          -> 1.96e+07x the current coefficient
    ```

    The drag repair multiplied a **damping** term by 1/dt = 1e4, where a larger coefficient is
    stabilising. This is 2e7 on a **body force**. The largest value ever run on this term is the
    intermediate 336 recorded in the comment above it, which drove a polar vertical runaway; the
    consistent value is **7.2× larger than that**. Shipping it on would not be a fix, it would be
    an experiment, and one with a recorded failure mode.

    *That "336" is itself Earth's*: `g/(omega*L_atm)` is 336 at `ATOM_Precipitation`'s
    `omega` = 7.29e-5 and `L_atm` = 400, and **1.97 here**. Item 34 described it as "a merely 336×
    larger coefficient", which understates the gap — against today's coefficient, 336 is already
    2.7e6×.

    **The knob corrects three things and they are inseparable.** Rescaling an incorrect force is
    not a repair, so `ATM_BUOY_CONSISTENT=1` changes all of:

    - **(a) the non-dimensionalisation** (item 34): `g*dt/u_0` → `g*L/u_0²`, RK4 supplying the `dt`.
    - **(b) the Boussinesq reference temperature — which is new, and is not a `dt` question at
      all.** The anomaly `(t − t_ref_level[i])` is `(T − T_ref)/t_0`, because `t` is `T/t_0`, so
      the shipped term divides the buoyancy by **`t_0` = 273.15 K**. Boussinesq buoyancy is
      `g·(T − T_ref)/T_ref`; non-dimensionally that is `(t − t_ref)/t_ref`, the `t_0` cancelling.
      The shipped form therefore overstates the force by `T_ref/t_0`: **5.49× at a 1500 K surface,
      0.96× at a 263 K top.** That is not a rescaling, it is a height-dependent distortion of a
      body force, strongest exactly where the convection is. `t_ref_level[i]` is the horizontal
      mean at each level and is the right base state; only the divisor was wrong. **This is the
      Earth-constant pattern in its purest form — `t_0` = 273.15 K is a non-dimensionalisation
      constant standing where a physical reference temperature belongs, and on a 288 K planet the
      two are within 5 % of each other.**
    - **(c) `L_atm`, not `L_coeff`.** `ATM_COEFF_SHELL` is an experiment on the damping and
      forcing *coefficients*, worth ~29×. Letting a body force join it would mean flipping a
      diagnostic switch changes the dynamics by 29× — precisely the confound that switch exists
      to remove.

    **Measured, 20 iterations, with (a) alone** (`run_buoy_off` / `run_buoy_on`):

    | | off | on (a only) | |
    |---|---|---|---|
    | max `ubud_buoy` | 0.000000 | 6.250198 | absent → present |
    | min `ubud_buoy` | −0.000002 | −34.569365 | ~4× the pressure gradient |
    | max `ubud_pgf` | 1.198185 | 10.575467 | grows 7× to oppose it |
    | max radial wind | 0.297520 m/s | 0.355395 m/s | **+19.5 %** |
    | `Psi_max` @ 20 | 106510.42 | 106510.62 | **+0.0002 %** |

    So the term goes from absent to dominant in the radial budget, the pressure gradient absorbs
    most of it, and the meridional overturning does not notice. **No runaway at that length — and
    that is not evidence of stability**: `buoyancy_ramp` was still 0.067 at iteration 20, so the
    force was at 6.7 % of its final strength. The measurement also predates (b), which cuts the
    near-surface term ~5.5×, so **it should not be quoted as the consistent term's behaviour.**

    **If the consistent form does blow up, that is a Boussinesq result, not a bug in this line.**
    CLAUDE.md's open risk stands: density spans ~2 orders of magnitude across this column and the
    solver rests on Boussinesq. A body force 2e7 larger is the first thing that would find out.

    **Bit-identity when off: NOT achieved, and stated rather than rounded away.** The off branch
    is the original expression verbatim — the first version of the knob hoisted a shared
    `buoy_coeff`, which is *not* bit-identical, because `*` and `/` are left-associative and
    `(A*g)*dt` ≠ `A*(g*dt)` in floating point unless `-ffast-math` chooses to reassociate. Even
    with the verbatim off branch, hoisting the result into a named `buoy_term` changes which
    `a*b+c` pairs may fuse into an FMA. Measured against the pre-knob binary at matched thread
    count over 200 iterations: **6 of 10 checkpoints identical, 4 differing in the last printed
    digit of Ψ (1e-7 relative)**, `residuum_atm` identical to 8 digits and the wind extrema to 6.
    For calibration, that is **smaller than the 1.2e-6 this file measures between 12 and 24
    threads.** True bit-identity is available by duplicating the whole `rhs_u` assignment under
    `if/else`; not done, because the duplication costs more than the 1e-7 buys.

    **Open**: a long run with the knob on, and the decision of what `buoyancy_ramp` should reach
    when the coefficient is 2e7 larger. Committed in `daaddf2`, default off, so every existing
    result stands.

51. **`cp_l` = 2040 J/(kg·K) is the mixture cp at ~1101 K, and the moist convection ran on it
    everywhere. Correcting it moves rain production by −18 % and nothing downstream at all.**

    `MoistConvection.h` used the constant `cp_l` throughout. Measured against the model's own
    `AtmMixture::cp_of` at the configured composition:

    | where | T | `cp_of` | `cp_l` is |
    |---|---|---|---|
    | isothermal skin — **where all condensation happens** | 263 K | 1562.8 | **1.305×** |
    | 218 km | 378.7 K | 1597.8 | 1.277× |
    | 186 km | 569.0 K | 1699.4 | 1.200× |
    | surface | 1500 K | 2341.2 | 0.871× |

    So the constant is **+30.5 % wrong exactly where the moist physics runs** and −12.9 % at the
    surface, and no single value can be right across a column that spans 1563–2341. (`cp_of` is
    clamped below 298 K, so the +30.5 % is a lower bound.)

    **Converted: the five sites where cp converts energy to temperature and does not cancel** —
    the surface buoyancy flux and `delta_T_sfp`, both θ_e exponents, the updraft moist adjustment
    (`G` and the latent heating, with `cp_u` refreshed each pass as `T_u` moves), and the `MC_t`
    tendency. The θ_e case was **already self-inconsistent**: `kappa` two lines above it used
    `AtmMixture::cp_of` while the exponent used the constant, in the same expression.

    **Not converted, deliberately: the nine `s ↔ T` sites.** `s = cp_l·T/s_0` and
    `T = s·s_0/cp_l` are an exact inverse pair, so `cp_l` **cancels** and never reaches a result.
    Making one side local would break the round trip and turn `s` into an implicit function of
    `T`. Documented in place, so the file has one rule and a stated exception rather than a mix.

    **Measured, 200 iterations, moist physics from iteration 0, 24 threads**, against a baseline
    differing *only* in this change (`run_drag_fixed` → `run_cplocal`):

    | | baseline | local cp | |
    |---|---|---|---|
    | max `S_r` (rain) | 2.692884 | 2.201288 | **×0.817** |
    | max `S_s` (snow) | 0.898367 | 0.734446 | **×0.818** |
    | max `S_v` | 0.157471 | 0.158541 | ×1.007 |
    | max `Q_Latent` | 0.218706 | 0.218709 | ×1.0000 |
    | cloud water | 15.198408 | 15.227884 | ×1.002 |
    | albedo | 0.4987 | 0.4987 | unchanged |
    | OLR | 281.55 | 281.61 W/m² | **+0.02 %** |
    | photosphere | 245.9 km / 270.36 K | identical | |
    | `Psi_max` | 98390.55 | identical | |

    **A −18 % correction to rain and snow production changes nothing integrated.** The albedo
    absorbs it, exactly as it absorbed `initCloudIce`'s H_crit repair (0.18 % on OLR, item 42):
    **the reflectivity saturates on the *presence* of condensate**, so how fast condensate is
    produced cannot reach the radiation. That is now three separate microphysical corrections
    that die at the same wall, and it is a statement about the albedo parameterisation, not about
    the microphysics.

    **A note on method, because the first answer was wrong.** The first comparison used
    `run_drag_base` and was confounded: that run predates `3e2d78f` and used 6 threads, and it
    made the cp change look like **+3 % on OLR** (273.33 → 281.61). Almost all of that was the
    **drag fix** (273.33 → 281.55). Same trap as items 45–46 — *a binary that changed
    mid-experiment*. What identified the correct baseline was `Psi_max` being bit-identical
    between `run_cplocal` and `run_drag_fixed`: a quantity the change cannot touch is the label
    that tells you which run you are actually comparing against.

    **Also found, documented not fixed**: `s_0` = 274515.75 = 1005 × `t_0` carries **Earth's
    dry-air cp**, not ATHAD's 2040 (which would give 557226), while `param.py` calls it
    "`cp_l * t_0`" so it reads as derived. ~~It changes no result precisely because `cp_l`
    cancels.~~ **RETRACTED THE NEXT DAY BY ITEM 52**: the cancellation holds for the explicit
    `s ↔ T` pair, and `rhsForcing` contains a third, *implicit* conversion — it multiplies the
    `s`-flux divergence by `t_0`, which inverts `s = cp_l·T/s_0` only if `s_0 = cp_l·t_0`. The
    convective transport heating is therefore **2.03× too large**. The claim rested on the one
    place the scale was thought to matter (`s` = 1.0 at `BC_Atm.h:219,368`) being dead under
    invariant 1, which it is — but it was not the only place.
    And two of the five converted sites currently multiply zero, because `Q_sensible_2D` is
    declared, reset, read and output but **never written**; converted anyway so the file has one
    rule. Committed in `e7bd455`.

52. **`s`, `s_u`, `s_d` are neither entropy nor dry static energy, and the scaling they do use
    is wrong in three places. One of them pins the convective heating at its safety cap for
    the whole spin-up, and it retracts item 51's "changes no result".**

    Audited on request: the units, dimensions and environment relations of the three `s`
    fields. What follows separates *what the quantity is* from *where its arithmetic is
    wrong*, because they are independent problems.

    **THE FIELD IS A NORMALISED TEMPERATURE, AND THE TREE CALLS IT TWO OTHER THINGS.**
    `MoistConvection.h` defines it once and inverts it once:

    ```
    s = cp_l * T / s_0          T = s * s_0 / cp_l
    ```

    That is `T` divided by a constant — dimensionless, correctly printed as `./.`. It is not
    entropy (which would be `cp ln T − R ln p`, in J/(kg·K)), and it is not dry static energy.
    Three names are in use for one array:

    | site | name used |
    |---|---|
    | `param.py:636` | "`s_0` — entropy at 0 °C, `cp_l * t_0` in J/kg" |
    | `cAtmosphereModel.h:717-719` | "dry static energy" / "…in the updraft" / "…in the downdraft" |
    | `Results_Atm.cpp:190` | "entropies in the up-and downdraft" |

    **THE MISSING TERM IS `g·z`, AND ON THIS SHELL IT IS THE LARGER HALF.** Dry static energy
    is `cp·T + g·z`; the geopotential is absent. That matters because a mass-flux scheme uses
    DSE precisely so that a non-entraining parcel cools at the dry adiabatic rate for free —
    conserving `cp·T` alone, a rising parcel keeps its temperature. On ATHAD's grid:

    | level | z | `g·z` [J/kg] | `cp·T` [J/kg] | `g·Δz/cp` per layer |
    |---|---|---|---|---|
    | i = 0 | 0 km | 0 | 3.5e6 (1500 K) | 5.1 K |
    | i = 20 | 54.7 km | 5.4e5 | 2.6e6 | 33.6 K |
    | i = 35 | 201.3 km | 2.0e6 | 8.0e5 | 103.6 K |
    | i = 39 | 277.2 km | 2.7e6 | 4.2e5 (263 K skin) | 139.9 K |

    `g·z` overtakes `cp·T` at about **157 km** and is **6.5×** it at the lid — and the skin
    above ~240 km is exactly where item 51 measured all the condensation. The adiabatic cooling
    the scheme omits is **5 K per layer at the surface rising to 140 K per layer at the top**.
    **The model already knows how to do this correctly elsewhere**: `findCloudBaseLFS` lifts a
    genuine θ_e-conserving parcel with local `cp` and pressure, and that is what sets the
    trigger, the LNB and CAPE. So the scheme carries two parcels — one that cools and decides
    *whether* to convect, one that does not and decides *how much*. **Not fixed here**: adding
    `g·z` changes what every `s` number means and needs its own measurement.

    **THREE PLACES WHERE THE SCALING IS WRONG.** All three are behind `ATM_MC_S_CONSISTENT`,
    default off, and the off branch is bit-identical (verified below).

    - **(1) The updraft recurrence divides by `s_0` when nothing needs it.**
      `dummy_s_u = (M·s_u + step·(E·s − D·s_u)) / s_0` — every term is already in normalised
      `s`, so the division makes `s_u` ~2.7e5× too small. **The tell is in the same block**:
      the four sibling recurrences beside it (`q_v_u`, `q_c_u`, `v_u`, `w_u`) have no such
      division.
    - **(2) The downdraft recurrence divides the whole bracket when one term needs it.**
      There, `L_latent·r_h·e_d` genuinely is J/(m³s) against its neighbours' kg/(m²s)·`s`, so
      **that term alone** wants the `/s_0`. It was applied to all three.
    - **(3) `MC_t` converts `s` back to kelvin with `t_0`.** `rhsForcing` multiplies the
      `s`-flux divergence by `t_0`, which inverts `s = cp_l·T/s_0` **only if `s_0 = cp_l·t_0`**.

    **AND (3) IS WHY ITEM 51 IS WRONG.** `s_0` = 274515.75 is `1005 × 273.15` — Earth's dry-air
    `cp` times `t_0`, exact to the digit, and correct upstream where `cp_l` **is** 1005. ATHAD
    raised `cp_l` to 2040 and left `s_0`, so the true inverse factor is `s_0/cp_l` = **134.57 K**
    and the shipped one is `t_0` = 273.15 K. Item 51 recorded this as harmless "because `cp_l`
    cancels". **The cancellation is real for the explicit `s ↔ T` pair and there is a third,
    implicit conversion that assumes the derived value** — so the convective transport heating
    is `cp_l·t_0/s_0` = **2.03× too large**. The in-code comment asserting harmlessness has been
    corrected in place. *The lesson is narrow and repeatable: to claim a constant cancels, find
    every conversion, including the ones written as a bare multiplication by something else.*

    **MEASURED — 40 iterations, 24 threads, moist physics from iteration 0, one binary, the
    knob the only difference** (`python/run_s_base` vs `run_s_fix`):

    | | shipped | `ATM_MC_S_CONSISTENT=1` |
    |---|---|---|
    | max `s_u` @ 10 | **342.33** → parcel T = **46 066 K** | 3.512 → **472.6 K** |
    | min `s_u` @ 10 | **−980.29** → **−131 915 K** | −3.863 (lid artefact, see below) |
    | max `c_u` @ 10 | **118.84 g/kg/s** | 0.0671 → **×1771 smaller** |
    | max `MC_t` @ 10 | **0.010000 K/s — AT THE CAP** | 0.000083 K/s (120× below it) |
    | max `MC_t` @ 20 / 40 | 0.003076 / 0.001492 | 0.001048 / 0.000425 (**×2.9 / ×3.5**) |
    | max `S_s` @ 40 | 0.435406 | 0.409573 (**−5.9 %**) |
    | max `S_r` @ 40 | 0.889990 | 0.888215 (−0.20 %) |
    | OLR @ 40 | 323.73 W/m² | 323.67 W/m² (**−0.02 %**) |
    | albedo, photosphere, `Psi_max` @ 40 | 0.4987, 237.4 km/368.00 K, 105590.83 | **identical** |

    **So the shipped scheme spends its spin-up with updraft parcels at tens of thousands of
    kelvin and `MC_t` pinned to `MCt_max`.** A safety cap was absorbing an arithmetic error, which
    is the worst way for one to hide: the field looks bounded and physical in every plot. The
    ×2.9–3.5 at iterations 20–40 is close to (3) acting alone — 2.03 — with the recurrences
    supplying the rest.

    **AND IT REACHES NOTHING INTEGRATED. FOURTH TIME.** OLR −0.02 %, albedo and photosphere
    identical, `Psi_max` identical to all eight digits. Same wall as item 51's local `cp`
    (−18 % on rain, +0.02 % on OLR) and item 42's `H_crit` (0.18 %): **the reflectivity
    saturates on the presence of condensate, so nothing about its amount or rate can reach the
    radiation.** ~~Anything that only changes microphysical rates is unmeasurable in this model
    until `albedo_cloud` responds to something.~~

    **THE SENTENCE ABOVE IS REFUTED BY ITEM 75 (2026-08-23), AND THE FOUR NULLS THAT PRODUCED IT
    WERE MEASURING AN ANNIHILATION.** The albedo clause is right — reflectivity does saturate on
    presence — but the conclusion drawn from it is wrong: condensate reaches the OLR through the
    LONG-WAVE cloud opacity `k_liq*LWP + k_ice*IWP`, which responds to AMOUNT. That path read as
    dead because `SaturationAdjustment` was manufacturing condensate into superheated cells and
    `ThreeCatIceScheme` deleted all of it every iteration (item 74), so there was no condensate
    left downstream for any rate change to act on. With it alive, one repair moved the OLR
    **−49.4 %** and the photosphere **44 km**. *Four consecutive nulls should have prompted a
    check that the field being varied still existed where it was supposed to matter.*

    **BIT-IDENTITY WHEN OFF: achieved and checked against a foreign run.** `run_s_base`
    reproduces `run_cplocal` (item 51's baseline, built before this knob existed) at iteration
    20 to every printed digit including `Psi_max` = 106510.37 — so the off branch is the shipped
    code and `nm` does not enter the physics.

    **WHAT IS CORRECT, stated because an audit that only lists faults is not one.** The flux
    form `M_u(s_u − s) + M_d(s_d − s)` is the right relation to the environment; the
    `|M| ≤ coeff_recurr` fallbacks that set `s_u = s`, `s_d = s` are right and are the reason
    the recurrence defect stayed survivable; `precompute` rebuilding `s` from the current `t`
    before anything reads it is right and documented; `downdraftEntrainment` seeding
    `s_d(i_LFS) = s(i_LFS)` is right; and `step[]` is built from `get_layer_height()`, so the
    module is on the radiation's true metric rather than the core's `exp_rm` (item 39).

    **FOUR NEIGHBOURING DEFECTS FOUND WHILE LOOKING, none fixed, all real:**

    - **`MC_t`'s second term carries an extra `t_0`.** `(L/cp)·conv_src` is already K/s before
      the `* m.t_0`, so that term is ~273× its own scale — and it is the term that holds
      `MC_t` at its **negative** cap in both arms of the A/B above. Inherited verbatim from
      `ATOM_Precipitation`, where it is equally wrong. Not folded into the knob because it is
      not an `s` question and would have confounded the measurement.
    - **The lid values of `s`, `s_u`, `s_d` are extrapolation artefacts.** `bcRadius`'s
      `both_cubic` list applies `x[iml] = x[iml-3] − 3x[iml-2] + 3x[iml-1]` to them at both
      ends. That is the exact overshoot the same function's comment says it removed *for
      `u`, `v`, `w`* — and it produces **negative `s`, i.e. a negative absolute temperature**:
      `s_u` = −3.86, `s_d` = −5.85, `q_v_u` = −1459 g/kg, every one of them at exactly
      300 005 m. `min s` = 0.575 → 77 K where the skin is 263 K. The same list is applied at
      the poles (`fields_cubic`).
    - **`q_v_u` reaches 3565 g/kg — 3.5 kg/kg of vapour in a mass fraction.** Interior, not a
      lid artefact (218 km). Its recurrence has a `≤ 0` clamp and no upper bound, unlike
      `q_c_u`'s `q_c_u_max`. Repairing `s` does not touch it: **both arms print the same
      3564.52**.
    - **`CAPE` divides a physical `step[i]` by `exp_rm`.** `computeCAPE` uses
      `g·step[i]/exp_rm·Δθ_v/θ_v`, mixing a true thickness with the core's quadratic-stretch
      Jacobian (item 39). It also uses environment T + a fixed `t_add_u` = 0.2 K as "the
      parcel" at every level, so this second CAPE never consults `s_u` either.

    **Also:** `s = s_u = s_d = 1.0` at initialisation (`UtilsAtm.h:442`) and in the dead
    `is_land` branches means **134.57 K**, not the intended 273.15 K. Harmless today — the
    cold-shutoff branches that also write 1.0 need `T ≤ 236.15 K` and the column's minimum is
    263.0 K — but it is the same `s_0` error, and it becomes live the moment anything cools
    the skin.

53. **The repair pass: item 52's four neighbouring defects fixed, the `s` scaling made the
    default, and a fifth defect found while repairing the third — which turned out to be the
    one that actually produced the impossible number.**

    Five builds, all 40 iterations, 24 threads, moist physics from iteration 0, one config,
    each adding to the one above it:

    | run | what it adds | |
    |---|---|---|
    | `run_s_base` | — | the shipped model |
    | `run_s_fix` | the three `s` scalings (item 52) | |
    | `run_s_repair` | `MC_t` term 2, the lid BC, the bounded recurrence, `computeCAPE` deleted | |
    | `run_s_all` | the `cc_factor` reference | |
    | `run_s_all2` | the recurrence narrowed to the shipped activation set | |

    **(1) `MC_t`'s second term lost its `*t_0`.** `(L/cp)·conv_src` is already K/s — `L/cp` is
    kelvin per unit mass fraction, `conv_src` is (kg/kg)/s — so the factor made that term 273×
    its own scale. Two terms in one sum cannot both be right when one is a temperature and the
    other a temperature over `t_0`. Fixed unconditionally, not behind the `s` knob, because it
    is not an `s` question. **Inherited verbatim from `ATOM_Precipitation`, where it is equally
    wrong.**

    **(2) The cubic extrapolation came off BOTH radial boundaries** for the convection scalars,
    microphysical rates, mass fluxes and diagnostic forces (`bcRadius`'s Pattern A list),
    replaced by a zero-gradient copy. This is the **third** time this stencil has had to come
    off a boundary in that one function — `v`/`w` "overshoot THROUGH zero", the turbulence
    scalars "amplify the concavity", `p_stat` reached −36 hPa — and the fourth argument was
    already written in `bcTheta` ("condition number ~7 … NaN at the pole"). Item 52 measured
    what it was doing at the lid; **ATHAD_COND supplied the surface half**: `s_u` = −5.012 at
    0 m, an updraft parcel at **−1020 K** over a 513 K sea. There was even a symptom-level
    patch already sitting under the i=0 loop — *"updraft moisture at the surface must be
    non-negative"*, clamping `q_v_u` and `q_c_u` — put there by whoever met the same negative
    values and treated them at the output instead of at the source.

    **(3) The updraft recurrence is bounded by construction, not by a clamp.** A clamp is what
    hid the `s` defect, so the repair is algebraic instead. `d(M·φ_u)/dz = E·φ − D·φ_u` with
    `dM/dz = E − D` expands to `M·dφ_u/dz = E·(φ − φ_u)`: **detrainment cancels exactly**, and
    φ_u can never leave the interval spanned by its previous value and the environment.
    Discretely that identity holds only if the denominator is the mass the numerator was built
    from, and three things break it here — `clamp_M` and the `is_land`/`t_00` kills rewrite
    `M(i)` after the fact; `E_u` is a moisture *convergence* and goes negative where the flow
    diverges moisture, which an entrainment rate cannot; and `step·D` can exceed `M(i-1)` on
    the deep layers, where `step` reaches 23 km. So the two weights are formed explicitly,
    floored at zero, and the sum is the denominator. `q_c_u` keeps the flux form — its bracket
    carries real sources, so it is not a convex combination of anything — but gets the
    consistent denominator.

    **(4) `computeCAPE()` deleted, not fixed.** It was a *second* CAPE that nothing read:
    `m.CAPE` was written here and never read anywhere in either tree. It divided a physical
    `step[i]` by `exp_rm`; its "parcel" was the environment temperature plus a fixed 0.2 K at
    every level, so it never lifted anything and never consulted `s_u`; and `m.CAPE` was a 1-D
    array indexed by LEVEL written inside a (j,k) loop, so the last column overwrote every
    other one. `cape_col[j][k]` in `findCloudBaseLFS` — a θ_e-conserving ascent with local `cp`
    and true thicknesses — is the real one and already correct. **One correct CAPE is better
    than one correct and one wrong.**

    **(5) THE FIFTH DEFECT, and it is the one that produced item 52's impossible number.**
    Repairing (3) did not move `q_v_u` at all: 3564.52 g/kg before and after, to the digit. The
    seed was the source, not the recurrence. `cc_factor` returns `q_sat(T)/q_sat(T_ref_cc)`
    with **`T_ref_cc` = 288.15 — Earth's mean surface temperature, a bare literal in a physics
    kernel.** On Earth that makes the factor O(1) by construction, which is the entire design
    ("warmer columns get proportionally more seed moisture, ~7 %/K"). Here the numerator is
    `q_sat(1500 K)` = 1, because the surface is supercritical (invariant 2), and the
    denominator is `q_sat(288 K)` at 250 bar ≈ 6.8e-5 — so the "ratio" was ~3.5e4 and the seed
    cap `q_v_u_add · cc_factor` came out at **3.47 kg/kg**. Referenced to the model's own
    surface temperature scale instead, the factor is O(1) on any planet and still carries the
    equator-to-pole contrast the design wants.

    **MEASURED:**

    | | `base` | `fix` | `repair` | `all` |
    |---|---|---|---|---|
    | max `q_v_u` @40 | 3564.52 g/kg | 3564.52 | 3564.52 | **724.67** |
    | min `q_v_u` @40 | −742.79 g/kg | −742.79 | **0.000000** | 0.000000 |
    | min `s_u` @40 | −3.865 (−520 K) | −3.865 | **0.000000** | 0.000000 |
    | min `s_d` @40 | −4.723 (−636 K) | −4.723 | **0.000000** | 0.000000 |
    | max `s_u` @10 | 342.33 (46 066 K) | 3.512 | 3.512 | 3.512 |
    | max `S_r` @40 | 0.889990 | 0.888215 | **0.000611** | 0.000611 |
    | max `S_s` @40 | 0.435406 | 0.409573 | 0.333249 | 0.400351 |
    | max `MC_t` @40 | 0.001492 | 0.000425 | 0.000425 | 0.000420 |
    | OLR @40 | 323.73 W/m² | 323.67 | 323.67 | **323.76** |
    | albedo / photosphere / `Psi_max` @40 | — | — | — | **identical in all four** |

    **`max S_r` = 0.889990 → 0.000611 is worth its own line**: that extremum was *at the lid*,
    300 005 m, i.e. the largest rain-production rate in the model was an extrapolation artefact.
    After the BC repair the maximum sits at 236 km, in the cloud.

    **And the OLR moves 0.03 % across all of it.** Fifth correction to die at the albedo
    (items 42, 51, 52 and both halves of this one): `Psi_max` is identical to eight digits,
    the albedo and photosphere unchanged. **Everything repaired here was unphysical and none of
    it was reaching the radiation**, which is a statement about `albedo_cloud` saturating on the
    presence of condensate, not about the repairs.

    **ONE THING GOT WORSE AND IS REPORTED AS SUCH.** At iteration 10 `max MC_t` sits at
    **0.010000 — its `MCt_max` cap** — in every build containing the non-`s` repairs, where
    `run_s_fix` had 8.3e-5. Iterations 20 and 40 agree to 0.3 % across all builds, so it is a
    spin-up transient, but it saturates a cap, and this file has just spent an item on what
    caps hide.

    **Attributed, and it is the bounded denominator.** Two hypotheses, both tested. *Refuted:*
    narrowing the recurrence to the shipped activation set (`run_s_all2`) leaves it at 0.010000,
    and `run_s_all2` is otherwise identical to `run_s_all` in every printed number at iteration
    40 — so the conservative form costs nothing and explains nothing. *Confirmed:*
    `ATM_MC_UNBOUNDED_UPDRAFT=1` with every other repair in place (`run_s_attrib`) returns
    **0.000083 at iteration 10 — `run_s_fix`'s value to the digit** — so it is the consistent
    denominator, not the boundary condition and not `MC_t`'s `t_0`.

    **The square closes** (`max MC_t` at iteration 10, all other repairs present in both
    right-hand runs):

    | | legacy denominator | consistent denominator |
    |---|---|---|
    | i=0 cubic | 0.000083 (`run_s_fix`) | **0.010000** (`run_s_all`) |
    | i=0 copy | 0.000083 (`run_s_attrib`) | **0.010000** (`run_s_final`) |

    — and outside that transient the two denominators agree exactly: `max MC_t` is 0.001050 at
    iteration 20 and 0.000420 at iteration 40 in *both* `run_s_attrib` and `run_s_all`. So the
    difference is confined to the spin-up, measured across two builds rather than inferred from
    one.

    **The i=0 half of the BC repair is confirmed and costs one number.** `run_s_final` — every
    repair, including the surface copy — is identical to `run_s_all` at iteration 40 in the OLR
    (323.76), `Psi_max`, albedo, photosphere, `S_r`, `S_s`, `q_v_u` and `s_u`, and differs in
    exactly one place: **`max s_d` 11.099406 → 11.085196 at 0 m**, i.e. the surface downdraft
    value goes from a cubic extrapolation (1493.7 K) to a copy of the level above it (1491.8 K),
    a 0.13 % change in the one quantity the surface boundary condition actually sets.

    **And that reframes it: nothing got worse, an Earth-calibrated ceiling started binding.**
    With the denominator consistent, `s_u` is a bounded, physically meaningful parcel value that
    genuinely differs from its environment, so the flux divergence it produces is a real
    convective heating instead of the near-zero one a rewritten mass was manufacturing. Where
    that heating exceeds `MCt_max` = 0.01 K/s the cap does what a cap does — but **`MCt_max` is
    36 K/hr, annotated in the code as "still 3× realistic", and next to it `inv_step_rh` floors
    the density at 0.01 kg/m³, "≈ air at 50 hPa"**. Both are Earth numbers, and at 256 km in a
    250 bar column both are being asked to bound something they were never sized for. That is
    the next thing to size, and it is a *new* open item, not damage from this repair.

    **The `s` repair is now the default** and `ATM_MC_S_LEGACY=1` restores the three defects
    together; `ATM_MC_UNBOUNDED_UPDRAFT=1` restores the old recurrence. Item 52's knob was
    default-off pending measurement; it has been measured.

    **NOT REPAIRED, deliberately, and both are stated in item 52:** the missing `g·z` — `s` is
    `cp·T`, dry static energy is `cp·T + g·z`, and `g·z` is the larger term above ~157 km here,
    so the updraft still does not cool as it rises. That is a redefinition of the field, not a
    units fix, and it needs its own measurement. `MC_t`'s **negative** cap, driven by `e_d`, is a
    mixed result: in ATHAD it binds at iteration 10 and is off the cap by 20 (−1.1e-5) both
    before and after, while **in ATHAD_COND it was at −0.010000 at every diagnostic before the
    repair and is −0.000716 after** — there, removing the 273× was enough to unbind it.

    **`bcTheta`'s `fields_cubic` was checked and was NOT a defect** — both its loops are plain
    copies; only the name was left over from the stencil the pole argument had already removed.
    Renamed `fields_pole_copy`, because a name is how the next reader decides where to look.


    **PORTED TO ATHAD_COND, and checked by reading that file rather than assuming the fork
    shares the code.** All of it applies there: same `/s_0` in both recurrences, same `*t_0` in
    both `MC_t` terms, same Pattern A cubic, same `computeCAPE`, same `cc_factor`. The
    magnitudes differ where `cp_l` does — `s_0` is Earth's `1005 × t_0` there too, against
    `cp_l` = 1349, so its `MC_t` transport ran **1.342×** hot rather than ATHAD's 2.03×, and its
    `s`-to-kelvin factor is `s_0/cp_l` = 203.50 K.

    **Measured there, 40 iterations, 24 threads, its own config, against a baseline run on the
    pre-repair binary** (`run_s_cond_base` vs `run_s_cond`), at iteration 10:

    | | baseline | repaired |
    |---|---|---|
    | max `s_u` | 5.545 → parcel **1128 K** | 2.502 → **509.1 K**, over a 513 K sea |
    | min `s_u` | −5.012 → **−1020 K**, at 0 m | **0.000000** |
    | min `MC_t` | **−0.010000 — its cap, at every diagnostic** | **−0.000716** |
    | max `MC_t` | 0.000000 | 0.000807 |
    | max `q_v_u`, `S_r`, cloud water | 333.18 / 0.0778 / 47.08 | **identical** |
    | max `c_u` | 0.067791 | **0.000000** |

    and at iteration 40, where the integrated quantities can be compared:

    | | baseline | repaired |
    |---|---|---|
    | max `s_u` | 5.737 → **1167.5 K** | 2.5014 → **509.0 K** |
    | min `s_u` | −5.012 | **0.000000** |
    | min `MC_t` | −0.010000 (cap) | −0.000457 |
    | max `c_u` | 0.067955 | 0.000000 |
    | OLR | 283.72 W/m² | **283.72** |
    | `Psi_max` | 42186.73 | 42186.38 (−8e-6) |
    | albedo, `S_r`, `S_s`, cloud water | 0.5000 / 0.0777 / 0.4628 / 47.01 | **unchanged to 5 digits** |

    **`c_u` → 0 is the most informative number in this table, and it is not a repair failing.**
    The baseline's updraft condensation was computed from a parcel at 1128 K — a temperature the
    `/s_0` defect invented. With the parcel at its correct 509 K it is sub-saturated and nothing
    condenses. **And it never will**, because the parcel does not cool as it rises: that is the
    missing `g·z` of item 52, and with the arithmetic now right it is no longer masked. So
    repairing the scaling has promoted the design defect from "documented" to "the only thing
    left between this scheme and a working updraft". `c_u` is one of the two terms in
    `conv_src`; in ATHAD it was already identically zero at every diagnostic in every build.

54. **`p_dyn` has no radial structure because the model has no radial force. Switching the two
    off-by-default metric terms on gives it one, and the radial pressure gradient appears
    immediately — 192× at mid-latitudes, at iteration 0.**

    Reported observation: `p_dyn` shows no variation in i, when the velocity field should give
    it some. Confirmed, diagnosed and then tested by turning the terms on.

    **THE FIELD IS 99.97 % A PRESCRIBED BALANCE.** Read out of the zonal VTK
    (non-dimensionalised, k = 87), the shipped model at iteration 0:

    | | latitudinal range | radial range j=45 | j=90 |
    |---|---|---|---|
    | iteration 0 | 444.61 | 0.247 (**0.06 %**) | 0.840 (0.19 %) |
    | iteration 20 | 444.18 | 0.166 | 0.769 |

    The radial structure is already negligible *before any dynamics have run*, and over 20
    iterations it **shrinks**. With `ATM_BALANCED_INIT=0` the field the velocity builds on its
    own is amplitude −0.077…+0.159 at iteration 20 — genuinely two-dimensional (radial range
    0.0179 at j=45, **16 % of its own amplitude**, non-monotonic) but **4000× smaller than the
    balance**. The flow's pressure response is not missing; it is swamped.

    **THREE COMPOUNDING CAUSES.**

    - **`initBalancedState` has no radial force to balance.** It builds
      `F_r = force_nd·nontrad·2sinθ·wbar (+ curvature)`, and `coriolis_nontraditional()` and
      `metric_curvature()` are both default-false (`lib/Utils.h:52,90`), so **`F_r ≡ 0`
      identically**. `balancedStateSolve` then minimises
      `(dp/dr·exp_rm − 0)² + (dp/dθ·inv_rm − F_the)²`, which *penalises every radial gradient*,
      and its operator is anisotropic against them: `cE = exp_rm²/dr²` = 400 at the surface
      against `cN = inv_rm²/dθ²` = 7.3, i.e. **55× stiffer radially**.
    - **The startup projection throws its own pressure away.** `project_initial_velocity` runs
      200 Jacobi sweeps, applies `v ← v − ∇p`, then **step 4 clears `p_dyn` to zero** so RK4
      does not double-correct. Measured: with the balance off, `p_dyn` at iteration 0 is
      identically 0 everywhere.
    - **One Gauss–Seidel sweep per iteration is a smoother, not a solve** — the file says so
      itself. With `num1` ≈ 400 against `num2` ≈ `num3` ≈ 155 it erodes radial structure
      preferentially, which is the 0.247 → 0.166 above.

    **THE TEST: `ATOM_CORIOLIS_NONTRAD=1 ATOM_METRIC_CURVATURE=1`, 40 iterations, 24 threads.**
    (`ATM_METRIC_RADIUS` defaults to 6370 km — metric r0 = 21.23 — so the documented
    prerequisite for the curvature terms is already met.)

    | | lat range | radial j=45 | radial j=90 |
    |---|---|---|---|
    | shipped, iter 0 | 444.61 | 0.247 | 0.840 |
    | **nontrad+curv, iter 0** | 443.66 | **47.42 (10.7 % of lat)** | **20.47 (4.6 %)** |
    | nontrad+curv, iter 20 | 443.18 | 44.55 | 19.88 |
    | shipped, iter 40 | 444.14 | 0.137 | 0.716 |
    | nontrad+curv, iter 40 | 443.05 | **43.37** | 19.88 |

    **192× at j=45 and 24× at j=90, present at iteration 0**, with the latitudinal range
    unchanged (444.61 → 443.66) — an addition of radial structure, not a rescaling. The profile
    is monotonic and physical: at j=45 `p_dyn` climbs 6.72 → 51.27 from surface to lid, while at
    the equator it *falls* 206.5 → 186.1. It decays only 6 % over 20 iterations against the
    shipped field's 33 %.

    **AND THE RADIAL MOMENTUM BUDGET BECOMES A BUDGET.** rms over the zonal slice at
    iteration 20:

    | term | shipped | nontrad+curv |
    |---|---|---|
    | `ubud_pgf` | 0.3025 | **18.257** |
    | `ubud_cor` | **0** | **18.901** |
    | `ubud_advv` | 0.00017 | 0.00013 |
    | `ubud_advh` | 0.00095 | 0.00086 |
    | `ubud_diff` | 0.000014 | 0.00034 |
    | `ubud_buoy` | 0.0000005 | 0.0000005 |
    | **NET** | **0.3027** | **2.183** |

    Shipped, the "budget" is one unopposed term: the net *is* the pressure gradient, because
    there is nothing for it to balance (item 42's finding, re-measured). With the terms on,
    a Coriolis term of rms 18.9 stands against a pressure gradient of rms 18.3 and they cancel
    to a residual of 2.2 — **88 % cancellation where there was 0 %**. Both grew ~60× from the
    shipped pgf; `p_dyn`'s radial gradient is now doing real work.

    **THE COST, AND IT IS GROWING — WITHDRAWN BY ITEM 58.** At 400 iterations, against a
    matched control, the metric arm ends **1.5 % ABOVE** the control after a −44 % excursion and
    a recovery. Everything in this paragraph is a transient. Original text follows.
    `Psi_max` 106510.37 → 101801.92 at iteration 20
    (**−4.4 %**) and 105590.83 → 96019.85 at iteration 40 (**−9.1 %**). The radial wind goes
    0.2975 → 0.2394 m/s (−19.5 %) — *weaker*, which is what a pressure gradient that is now
    opposed rather than acting alone should give. Everything else is unmoved at 20: `v`
    −0.02 %, `w` +0.01 %, OLR identical at 338.84, max T to 7 digits. Invariant 1 holds — N–S
    asymmetry −3.75e-07 at iteration 40 against the 1e-5 tolerance, though that is 2.4× the
    shipped run's −1.56e-07 and worth watching. **No runaway at this length, and 40 iterations
    is not a stability result** — item 28 records a polar vertical runaway from a radial force
    that was 300× too large, and the Ψ trend here is *doubling* between iteration 20 and 40, so
    where it lands is unknown. This file has been caught extrapolating a monotone trend twice;
    this is not a third.

    **A DIAGNOSTIC BLIND SPOT FOUND IN THE ACT OF USING IT.** `ubud_advh` is built from the raw
    advective pieces, `-(v_invrm·dudthe_adv + w_invrs·dudphi_adv)`, while
    `metric_curvature()` adds `−(v² + w²)·inv_rm` to `transport_u`. **So with the curvature
    terms on, one live term of the radial budget is captured by no `ubud_*` field at all**
    (≈0.5 by estimate, against pgf 18 — small, but the budget cannot be closed without it).
    The instrument was built when the term was guaranteed inert, and it inherits that
    assumption. Fix `ubud_advh` before using the budget to judge the curvature terms.

    **NOT FLIPPED.** These stay default-off: what is measured here is 40 iterations of an
    initialisation-dominated transient, `lib/Utils.h` calls both terms "small corrections" on
    the evidence of a 20-iteration Earth test, and the one number that moved — Ψ, −4.4 % — is
    the number this file has spent items 26-37 trying to interpret. The honest next step is a
    long run, and a repaired `ubud_advh` to read it with.

55. **`N²` gets its shape from the prescribed temperature profile and from the truncation error
    of the integration that writes it. It contains no dynamical information at all — measured
    across three runs whose `p_dyn` differs by 2500×.**

    Asked where a diagnostic that no prognostic equation produces gets its structure from. The
    chain is short and entirely upstream of the dynamics.

    **IT IS COMPUTED FROM WHAT THE SAME FUNCTION JUST WROTE.** The `N²` block sits at the end of
    `ThermoAtm::densities()` (`ThermoAtm.h:1592`), and it reads `m.t` and `m.p_stat` — both of
    which `densities()` assigned a hundred lines earlier in the same call:

    ```
    t_surf_equator / t_surf_pole   (prescribed)
        -> T_ad = T_prev - (g/cp_of(T_prev)) * dz        the dry adiabat, forward Euler
        -> T_i  = max(t_skin, T_ad)                      isothermal above ~240 km
        -> p_i  = p_prev * exp(-g*dz/(R*T_mean))         hydrostatic on that same T
        -> theta = T*(p_0/p)^(R/cp),  N^2 = (g/theta) dtheta/dz
    ```

    So `N²` is a re-reading of the prescribed profile. Nothing in the momentum or pressure
    solution enters it, and `ATM_PROGNOSTIC_T=0` (the default) guarantees the profile the
    dynamics and the radiation computed is overwritten before `N²` is formed.

    **MEASURED: IT IS IDENTICAL ACROSS RUNS THAT DIFFER ENORMOUSLY IN THE FLOW.** Iteration 20,
    same zonal slice:

    | run | `p_dyn` radial range | `N²` at i=8 / 32 / 39 | max abs difference vs shipped |
    |---|---|---|---|
    | shipped | 0.166 | 6.00e-07 / 1.515e-05 / 2.737e-04 | — |
    | nontrad+curvature | **44.55** | 6.00e-07 / 1.515e-05 / 2.737e-04 | **8e-08** |
    | no balanced init | 0.0179 | 6.00e-07 / 1.515e-05 / 2.737e-04 | 2e-08 |

    **2500× in the radial pressure structure and 4 % in Ψ move `N²` by 0.5 % at its smallest
    value and 0.03 % at its largest.** Whatever `N²` is reporting, it is not the circulation.

    **AND BELOW THE SKIN IT IS NOT STRATIFICATION, IT IS THE GRID.** An exact adiabat has
    `N² ≡ 0`. The model's is not exact: it steps `T_ad = T_prev − (g/cp(T_prev))·dz` with `cp`
    taken at the *bottom* of each layer, and `cp` falls as `T` falls (2341 → 1563 across this
    column), so every layer comes out slightly too warm and θ drifts upward. The signature is
    unambiguous:

    | i | z [km] | dz [km] | `N²` | `N²/dz` | `N²/dz²` |
    |---|---|---|---|---|---|
    | 2 | 2.5 | 1.37 | 1.3e-07 | 9.5e-11 | 6.9e-14 |
    | 8 | 12.9 | 2.15 | 6.0e-07 | 2.8e-10 | 1.30e-13 |
    | 20 | 54.7 | 5.29 | 2.96e-06 | 5.6e-10 | 1.06e-13 |
    | 29 | 122.6 | 10.39 | 9.50e-06 | 9.1e-10 | 8.8e-14 |
    | 35 | 201.3 | 16.29 | 2.50e-05 | 1.53e-09 | 9.4e-14 |

    Over a **12× range in layer thickness**, `N²/dz` varies by 16× and **`N²/dz²` is constant to
    ±30 %**. That is a second-order truncation error, not a physical stratification: the
    monotonic rise of `N²` up the column is the exponential grid stretch, nothing else. On a
    uniform grid it would be flat, and refining `zeta` would shrink it quadratically.

    **THE SKIN VALUE IS THE ONE GENUINE NUMBER, AND IT IS A CHECK.** In the isothermal layer
    `N² = g²/(cp·T)` exactly. Measured at i = 38: `N²` = 2.372e-04 with T = 263.0 K implies
    **cp = 1543 J/(kg·K)**, against `AtmMixture::cp_of`'s clamped value of ~1563 at 298 K —
    **1.3 %**. So the diagnostic is arithmetically correct; it is the profile it reads that is
    an input. (At the lid i = 40 the one-sided difference drops it to 1.93e-04; edge artefact,
    not physics.)

    **WHAT THIS COSTS THE ARGUMENT IT WAS BUILT FOR.** Item 42 added `N²` to stop the model
    *asserting* that it is neutrally stratified by construction, because that assertion carries
    the claim that no baroclinic eddy can maintain an indirect cell. The instrument now answers:
    the column is neutral **because the adiabat is imposed**, and the only departure from
    neutrality is the integrator's own error. That is not evidence about the atmosphere. It
    becomes a measurement the moment `ATM_PROGNOSTIC_T=1` — and CLAUDE.md's "N2 confirms
    invariant 4" should be read as "N2 confirms `densities()` integrates its own adiabat to
    O(dz²)", which is a real and useful check, but a different one.

56. **CO₂ shows no convective influence for three independent reasons, any one of which would
    be sufficient — and the third one is a bookkeeping choice that makes CO₂ inert against
    everything, not just convection.**

    Measured first: `max co2` = `min co2` = **0.205300 kg/kg** at every diagnostic of every run
    in this file. Not approximately uniform — uniform to every printed digit.

    **(1) THERE IS NO CONVECTIVE TERM IN ITS EQUATION.** Side by side:

    ```
    rhs_c.x   = -transport_c   + diffusion_c + coeff_trans*S_v*r_humid + coeff_MC_q*MC_q
    rhs_co2.x = -transport_co2 + diffusion_co2
    ```

    There is **no `MC_co2` field anywhere in the tree**. `MoistConvection` carries five updraft
    scalars — `s_u`, `q_v_u`, `q_c_u`, `v_u`, `w_u` — and their downdraft partners, and it has
    no tracer slot at all. CO₂ appears in that file only as a thermodynamic argument to
    `cp_of`/`R_of`, never as something transported. **So a CO₂ field WITH gradients would still
    be invisible to the convection scheme**; convection could reach it only indirectly, through
    the resolved velocity field it modifies.

    **(2) THERE IS NOTHING TO TRANSPORT, AND `rhs_co2` IS IDENTICALLY ZERO.** `co2Atmosphere`
    fills the field uniformly (initial condition only since item 12) and no source or sink for
    CO₂ exists anywhere in the model. For a uniform field `transport_co2 = u·∂co2/∂r + … ≡ 0`
    and `diffusion_co2 = ∇²co2 ≡ 0`, so `rhs_co2 ≡ 0` in every cell for all time. CO₂ is not a
    slowly-evolving tracer that happens to stay flat; it is **frozen by construction**, which is
    why the printed min and max agree to six decimals rather than to noise.

    *The `co2 column` diagnostic does show structure — 523353.42 against 523352.74 kg/m², a
    7e-7 spread — but that is the AIR column varying with the surface temperature, not CO₂.*

    **(3) AND IT CANNOT EVEN DILUTE, BECAUSE THE COMPOSITION CLOSES ON THE BACKGROUND.**
    `AtmMixture::split` sets `q_b = 1 − q_v − q_c`: the background is the residual. Water is 67 %
    of this atmosphere's mass and it is **not** uniform — 536.9 to 739.0 g/kg across the domain —
    so the total mass of a parcel changes as water moves in and out of it, and the mass fractions
    of everything else must respond. In this model **all of that response lands on the
    background** and none on CO₂:

    | | q_v | q_CO₂ | q_bg | R_mix |
    |---|---|---|---|---|
    | wettest cell | 0.7390 | 0.2053 (pinned) | 0.0557 | 398.9 |
    | reference | 0.6724 | 0.2053 | 0.1223 | 387.9 |
    | driest cell | 0.5369 | 0.2053 (pinned) | 0.2578 | 368.4 |

    The background fraction spans **4.6×** while CO₂ does not move at all. Concentrating CO₂ and
    background *together* — which is what removing water from a parcel physically does — puts
    `R_mix` **1.4 % lower at the wettest point and 2.6 % higher at the driest**, and that error
    is correlated with the water field, so it is a systematic, water-shaped bias in `R_mix`,
    `cp_mix` and hence in the density and the hydrostatic column, not noise.

    **So "CO₂ is prognostic" (item 12) is true and empty at the same time.** The transport
    equation is integrated, its result is exactly the initial condition, and nothing in the
    model — convection, microphysics, or dilution by the water it shares a mass budget with —
    can change it. That is a correct outcome for reasons (1) and (2), which are consequences of
    there being no CO₂ source at 4.4 Ga in this model; reason (3) is a modelling choice worth
    revisiting, because it is the one that would otherwise have given CO₂ a real spatial
    structure to transport.

    **What would make CO₂ move**, in increasing order of work: renormalise the non-water
    carrier so CO₂ and background concentrate together (fixes (3), and fixes the `R_mix` bias
    whether or not CO₂ is interesting); add a CO₂ slot to the mass-flux scheme (fixes (1)); give
    CO₂ a source — magma-ocean degassing or dissolution — which is the only thing that makes
    (2) false and the only one that is a Hadean science question rather than a plumbing one.

57. **The non-water carrier is renormalised, so CO₂ and background now concentrate together as
    water leaves a parcel — and the CO₂ initial condition becomes `initCO2()`, where it
    belongs. Item 56's third reason is fixed; the stored field keeps its numbers.**

    **THE CONVENTION.** The transported `co2` field is read as the CO₂ mass fraction **at the
    reference water content**, and `AtmMixture::q_CO2_of` scales it to the local carrier:

    ```
    q_CO2 = co2_stored * (1 - q_v) / (1 - q_H2O_ref)
    ```

    At `q_v = q_H2O_ref` this is the identity, so `co2_0`, the initial condition and every
    printed value keep their meaning and their numbers. Away from it, `q_CO2` and `q_bg` scale
    **together** and the composition of the carrier is invariant — which is the physical
    statement item 56 said was missing. CO₂ stays prognostic: a transported perturbation in the
    stored field still shows. `carrierRef()` is set once inside `resolve()`, the one function
    that knows the configured composition and runs before any physics; left at 0 it selects the
    legacy path, so nothing can depend on call order. `ATM_CO2_DILUTE=0` restores the old
    behaviour.

    **PUTTING IT IN `split()` FIXES SIX FUNCTIONS AT ONCE** — `R_of`, `cp_of`, `M_of`,
    `M_nonwater`, `x_CO2_of`, `x_H2O_of`. Four sites read `co2` raw and were converted
    individually:

    - **`MultiLayerRadiation`** — the CO₂ optical depth, and the one with teeth. Pinning `q_c`
      handed the whole water-driven variation to `q_b`, i.e. it moved mass between two species
      whose opacities differ by **1000×** (`kappa_CO2` = 1e-3 against `kappa_bg` = 1e-6).
    - **`ThermoAtm::co2Column`** — so the integral is a mass, not a mass fraction integrated as
      though it were one.
    - **`IceSchemeCommon` ×2** — the `c ≤ 1 − co2` ceiling. Under the new convention
      `q_v + q_CO2 + q_bg = 1` for *any* `q_v`, because the carrier shrinks as water grows, so
      that ceiling has no meaning and 1.0 is the real bound on a mass fraction. (It had already
      stopped firing — item 23.)

    A side effect worth noting: `M_nonwater` now returns a value **independent of the water
    content**, which is what its own docstring always claimed it computed ("renormalised to
    exclude H₂O"). Under the old split it varied with water.

    **MEASURED — 40 iterations, 24 threads, one binary, the knob the only difference:**

    | | `ATM_CO2_DILUTE=0` | dilution on |
    |---|---|---|
    | co2 column, max − min | **6.8e-04 kg/m²** (1.3e-6 relative) | **220.0 kg/m²** (4.2e-4) |
    | co2 column average | 523352.862 | 523314.417 → 523314.708 |
    | max water vapour | 744.965 g/kg | 777.730 (+4.4 %) |
    | OLR | 323.83 W/m² | 323.87 (+0.012 %) |
    | photosphere T | 368.00 K | 368.30 K |
    | mean albedo | 0.4987 | 0.4988 |
    | `Psi_max` | 105590.82 | 105590.66 |
    | max temperature | 1219.345378 | 1219.345314 |

    **The CO₂ column gains 320× more spatial structure** — that is the whole point, and it is
    the direct answer to item 56: CO₂ mass now follows the inverse of the water field instead of
    being flat to 1e-6. The global mean falls 0.0073 %, because the domain-mean water sits
    slightly above `c_0` so the mean carrier is slightly smaller, and it then **holds to
    5.6e-07 across the run's diagnostics** — the mass is conserved, which was the number most
    likely to misbehave under a convention change.

    **BIT-IDENTITY OF THE LEGACY PATH, checked against a foreign binary.** `ATM_CO2_DILUTE=0`
    reproduces `run_s_final` — built before any of this existed — digit for digit at iteration
    10: water vapour 696.916724/611.243205, co2 column 523435.587906/505987.734976, column
    average 520130.791.

    **A misattribution caught in the act.** The `co2 column` field shows a 3.4 % latitudinal
    spread at iteration 10 in *both* arms, and it would have been easy to publish that as the
    dilution working. It is a pole artefact of that diagnostic, present identically with the
    knob off. The real signature is the max−min at iteration 40 and the mean.

    **`initCO2()` — the initial condition moves to `InitValues_Atm.cpp`.** `ThermoAtm::
    co2Atmosphere()` is retired. It never belonged in a thermodynamics class: it is an initial
    condition, and keeping it there is what made it easy to call from inside the time loop,
    which is exactly the defect item 12 had to remove. The new routine states the two contracts
    that were previously implicit — that it must precede `initTemperatureData` and `densities()`
    (item 22: with `co2` still zero, R collapses to the background's 317.3 against 387.9, an
    18 % short scale height, and the lid snapshot the run is pinned to is built from it), and
    what the stored value now means.

    **And it can make the tracer non-uniform, which is what item 56 said was untested.**
    `ATM_CO2_INIT_PERTURB` (default 0.0, bit-identical) lays down a surface-rich, top-poor
    gradient, linear in **true height** rather than grid index — on this stretch those are very
    different things (item 39). Verified at amplitude 0.25: 0.2566 kg/kg at the surface against
    0.1618 at 277 km. **It is a measurement tool and not a Hadean claim**: item 56's second
    reason — a uniform tracer with no source has nothing to transport — means the CO₂ transport
    and its convective redistribution are *untested*, not tested and found working. Anything
    measured with this knob set describes the numerics, and no number from such a run belongs in
    a statement about the atmosphere.

    **NOT FIXED, and it is second order.** `c` is water **vapour**; condensate lives in separate
    arrays, so when vapour condenses the parcel keeps that mass as liquid or ice and the
    gas-phase carrier grows less than `1 − c` implies. The strictly correct normalisation is by
    the gas mass, `1 − (cloud + ice + graupel)`. Condensate here runs 12–47 g/kg, so this is a
    ≤5 % residual on a correction that is itself ±15 %.

58. **The 400-iteration pair: the two metric terms do NOT cost the circulation. Ψ takes a deep
    excursion, restructures, and comes back level with the control — and every intermediate
    number said the opposite.**

    Item 54 measured `ATOM_CORIOLIS_NONTRAD=1 ATOM_METRIC_CURVATURE=1` at 40 iterations, found
    Ψ_max down 4.4 % at 20 and 9.1 % at 40, and left the endpoint open on the grounds that a
    doubling gap over one interval is not a trend. It was right to. Two matched 400-iteration
    runs, 24 threads, same binary, same config, the control under `ATM_CO2_DILUTE=0` so its
    physics matches the metric arm's:

    | iter | control (shipped) | nontrad + curvature |
    |---|---|---|
    | 20 | 106510 @ 5° | 101802 @ 5° |
    | 100 | 102877 @ 5° | 79492 @ 5° |
    | 180 | 99291 @ 5° | **59370 @ 5°** |
    | 260 | 95674 @ 5° | 67330 @ **17°** |
    | 340 | 92014 @ 5° | 78189 @ 16° |
    | 380 | 90168 @ 5° | 85995 @ 9° |
    | 400 | **89241 @ 5°** | **90606 @ 9°** |

    **The control decays monotonically and almost linearly, −16.2 % over 380 iterations**, its
    maximum pinned at lat 5° throughout. That reproduces item 37's −16.7 % / −16.2 % on the
    repaired binary, which is worth having on its own: items 52-53 did not disturb the
    circulation.

    **The metric arm does something structurally different.** It falls to a minimum near
    iteration 180-200 — **−44 %**, far worse than the control — and then *recovers*, ending
    **1.5 % ABOVE** the control. The cell core migrates out to 17° and back to 9° while it
    happens. That is a reorganisation, not a decay: the terms take the prescribed cell apart and
    something else assembles.

    **SO ITEM 54'S "COST, AND IT IS GROWING" IS WITHDRAWN.** The sequence of readings was
    −4.4 % at 20, −9.1 % at 40, −25 % at 100, −40 % at 180, and **+1.5 % at 400**. Every
    intermediate value pointed the wrong way, and the two that looked most like a trend —
    the doubling between 20 and 40, and the near-linear fall to 180 — were the least
    informative. **This is the third time this file has been caught by a monotone trend**, and
    the only reason it did not become a claim is that the run was paired with a control of the
    same length on the same binary from the start. *The lesson is not "extrapolate more
    carefully". It is that an unpaired long run cannot answer a question of this shape at all.*

    **What is NOT settled.** 400 iterations is 78 s to 25 min of physical time (item 47), so
    neither arm is near equilibrium and "they converge" is a statement about this window. The
    two curves are still moving in opposite directions at 400 — the control down 927 per 20
    iterations, the metric arm up 4611 — so the crossing is recent and the next hundred
    iterations could put them anywhere. **Do not read the +1.5 % as a result either.** What is
    solid is the negative: the 4-9 % losses item 54 reported are transient, and the terms
    cannot be rejected on the evidence that rejected them.

59. **The gas-mass normalisation moves into `split()`, and `densities()`'s `water_factor` is
    deleted in the same edit — because the model already had this correction, in one place,
    partial, and applying both would have counted the condensate twice.**

    Item 57 left this as the known second-order residual: `c` is water **vapour**, so
    normalising the carrier by `1 − q_v` hands the suspended condensate's share of the mass to
    the gas. A droplet is still in the parcel and still weighs something.

    **THE MODEL ALREADY KNEW.** `ThermoAtm::densities()` carried

    ```cpp
    const double water_factor = std::max(0.5, 1.0 - m.cloud.x[i][j][k] - m.ice.x[i][j][k]);
    m.r_humid.x[i][j][k] = 1e2 * p_i / (R_loc * T_i * water_factor);
    ```

    — the gas-mass normalisation, applied to `r_humid` and to nothing else, spelled without
    graupel and floored at 0.5. **One concept with two implementations is how it stayed
    partial**: every other consumer of the composition — the radiation's optical depth, the
    saturation carrier, the mixture `R` and `cp`, the hydrostatic step three lines above — used
    a carrier that silently included the condensate.

    **THE FIX.** `split()` takes `q_cond` (cloud + ice + **graupel**), defaulted to 0 so every
    existing call compiles unchanged, and returns fractions that **sum to `1 − q_cond`, not to
    1**. That is the point rather than a wart: they are per unit TOTAL parcel mass, so `R_of`
    returns `(1 − q_cond)·R_gas` and `p/(R_of·T)` is the total density directly — which is
    exactly what `water_factor` was constructing by hand. With the factor inside `R_loc`, the
    divisor had to go, and the 0.5 floor moves onto `R_loc` where it belongs (a parcel that is
    all water leaves no carrier; the guard is against dividing by zero, not against the
    physics).

    **AND IT FIXES THE HYDROSTATIC STEP AS A BY-PRODUCT.** `p_i = p_prev·exp(−g·dz/(R_loc·T))`
    sits three lines above the deleted divisor and used the un-normalised `R_loc`, i.e. it
    integrated `dp/dz = −ρ_gas·g` where the correct weight is `ρ_total`. `water_factor` never
    reached it. Now it does.

    **WHERE `q_cond` IS AND IS NOT PASSED, as a rule rather than a list**: pass it where the
    result multiplies or produces a **total density**, because there `(1 − q_cond)` is exactly
    the gas-to-total conversion. Do **not** pass it to a per-mass gas property. The adiabat's
    `cp_of` is the case that matters: `(1 − q_cond)·cp_gas` would make the lapse rate *steeper*,
    which is backwards — suspended condensate adds its own ~4200 J/(kg·K) and makes a parcel
    harder to cool. Doing that properly is a moist adiabat, which is physics, not
    normalisation. The site is commented to that effect.

    **MEASURED — 40 iterations, 24 threads, one binary, against the dilution-only run:**

    | | dilution only | + gas mass |
    |---|---|---|
    | OLR | 323.86 W/m² | **323.52** (−0.11 %) |
    | imbalance | −52.79 W/m² | −52.45 |
    | co2 column average | 523314.529 | 523318.560 (+7.7e-06) |
    | max water vapour | 777.932 g/kg | 783.164 (+0.67 %) |
    | photosphere T | 368.30 K | 368.28 K |
    | mean albedo | 0.4988 | 0.4988 |
    | `Psi_max` | 105590.66 | 105588.74 (−1.8e-05) |

    **The OLR moves nine times further than item 57's dilution did** (−0.11 % against
    +0.012 %), which is the opposite of what "second-order residual" suggests. The reason is
    where the condensate is: it peaks at 37-47 g/kg in a thin band near 236 km, which is also
    where the photosphere sits, so a correction that is negligible in the deep column lands
    exactly on the levels the outgoing flux comes from. **Both numbers are 40-iteration
    transients and neither is a converged result** — what they establish is that the correction
    is live and its size is set by the condensate's location, not its column mean.

    At iteration 10 the two runs agree to every printed digit and differ only in `residuum_atm`
    at 2e-06; the separation grows through the run. A coarse extremum agreeing is not evidence
    that nothing changed.

60. **The background was six gases wearing one opacity, and two of them are not inert. Splitting
    it raises the background's optical depth by 1867× — and the converged OLR moves 0.08 %.**

    `kappa_bg` = 1e-6 m²/kg is the value of a radiatively inert diatomic, and it was applied to
    all six background gases at once. Two of them are **NH₃ and CH₄ at 1.4 mole-% each**, among
    the strongest infrared absorbers in the set, and **SO₂** is another 1.4 % with strong bands
    at 7.3 and 8.7 µm.

    **THE SPLIT COLLAPSES TO ONE NUMBER, WHICH IS WHY IT IS CHEAP.** All six are well mixed and
    source-free, so their mass ratios are fixed by the composition and never evolve. The
    per-species sum `Σ f_i·κ_i` is therefore a constant, and the radiative transfer does not
    change shape — only its coefficient does. Printed every diagnostic by `Results_Atm`:

    | species | fraction of bg | κ [m²/kg] | share of κ_bg_eff |
    |---|---|---|---|
    | N₂ | 0.3207 | 1.0e-06 | 0.0 % |
    | CH₄ | 0.0857 | 3.0e-03 | 13.8 % |
    | **NH₃** | 0.0910 | 1.0e-02 | **48.7 %** |
    | H₂ | 0.0108 | 1.0e-05 | 0.0 % |
    | CO | 0.1496 | 1.0e-04 | 0.8 % |
    | **SO₂** | 0.3422 | 2.0e-03 | **36.7 %** |
    | **TOTAL** | 1.0000 | **1.867e-03** | **1867× the lumped value** |

    **κ_bg_eff = 1.87e-3 is nearly 2× κ_CO2**, and 32 % of the background by mass — N₂, H₂ and
    CO together — contributes under 1 % of it. The lumped value was not an average of the six;
    it was the value of the least absorbing one.

    **AND THE OLR BARELY NOTICES, WHICH IS THE POINT.**

    | iteration | lumped | split |
    |---|---|---|
    | 20 | 338.84 W/m² | 337.54 (**−0.38 %**) |
    | 40 | 323.76 W/m² | 323.49 (**−0.083 %**) |

    Photosphere 368.41 → 368.69 K at iteration 20 and 368.05 → 368.30 at 40, height unchanged at
    237.4 km; `Psi_max` 105590.83 → 105588.89; albedo unchanged. **The effect SHRINKS with
    iteration** — a 1867× change in one opacity term is worth 0.38 % of the OLR at iteration 20
    and 0.08 % by 40.

    **That is items 29 and 41 measured a third time, and the cleanest instance yet.** Those found
    the converged OLR moving 0.10 % over 64× in `kappa_H2O` and 0.36 % over a 6× grid stretch,
    and concluded the OLR is `σT_skin⁴` read back out. This is the same conclusion reached by
    changing a *different* opacity, in a *different* species, by a *different* mechanism — and
    the OLR still returns to the same place. **Do not read the 1867× as a licence to expect a
    radiative response; read it as one more measurement of `t_skin`.**

    Arithmetically the smallness is unsurprising once written down: the background is 12 % of the
    mass, so its τ contribution goes from 1.2e-07 (negligible) to **2.3e-04, now comparable with
    CO₂'s 2.1e-04** — but both are ~3 % of water's 6.7e-03. Total τ rises 3.3 %. **The structural
    correction is large and the radiative consequence is small, and both statements are true.**

    **THE SIX KAPPAS ARE ASSUMPTIONS OF THE SAME STANDING AS `kappa_H2O` AND `kappa_CO2`** — grey
    band-averages chosen for the right order and the right *ordering*, not measured, and set
    relative to this scheme's own two anchors. A grey scheme cannot represent the window regions
    these gases fill differently, so **the split is better founded than any individual number in
    it**: that NH₃ and CH₄ are not N₂ is certain, and that NH₃'s κ is 1e-2 rather than 3e-3 is not.
    `ATM_BG_LUMPED=1` restores the single `kappa_bg` and every number that predates this item.

    **In ParaView and Results.** The six mass fractions `q_N2 … q_SO2` are written to all three
    VTK slices, computed on the fly rather than stored: `q_i = q_bg·f_bg[i]` exactly, so six 3-D
    arrays would carry one array's worth of information at ~128 MB. `q_bg` comes from
    `AtmMixture::split`, the routine the thermodynamics and the radiation use, so the plotted
    field cannot drift from the integrated one. Verified against the table: NH₃ 0.01478/0.16245
    = 0.0910, SO₂ 0.3422, N₂ 0.3207. `q_bg` itself now spans 0.081–0.162 — item 57's carrier
    dilution, visible for the first time.

61. **`g·z` goes into the static energy, and the updraft cannot use it — because the updraft is
    ONE GRID LEVEL DEEP. Item 53's "the only thing between this scheme and an updraft that
    condenses" is retracted: the blocker is geometry, not thermodynamics.**

    **(a) The geopotential is written.** Item 52 found `s`, `s_u`, `s_d` holding `cp_l·T/s_0` —
    a normalised temperature — while `cAtmosphereModel.h` called them dry static energy, which is
    `cp·T + g·z`. `ATM_MC_GEOPOTENTIAL=1` adds `g·z(i)/s_0` at every write of the three arrays,
    subtracts it at the one place a parcel temperature is read back out
    (`T_u = (s_u − φ)·s_0/cp_l`), and moves the two `<= 0` floors onto `φ(i)` so they stay floors
    on temperature rather than on DSE. Default OFF.

    It is the larger term here, as item 52 said: `g·z` = 2.94e6 J/kg at the lid against
    `cp·T` = 4.2e5. `s` stops looking like a temperature profile and starts looking like a static
    energy — **1.567–11.108 becomes 10.229–14.268**, nearly uniform, which is what a column on
    its own adiabat must have.

    **MEASURED — 40 iterations, 24 threads, ONE binary, moist physics from iteration 0:**

    | | off | on |
    |---|---|---|
    | `s` range (iter 0) | 1.567 – 11.108 | 10.229 – 14.268 |
    | `max s_u` (iter 40) | 3.6578 | 11.2577 |
    | `max MC_t` (iter 40) | 1.078e-03 K/s | **2.06e-04 (−81 %)** |
    | `min MC_t` (iter 40) | −2.5e-05 K/s | **−8.31e-04 (33×)** |
    | **`max c_u`** | **0.000000** | **0.000000** |
    | `max q_v_u` | 748.345 g/kg | 748.413 |
    | `max M_u` | 135.026 g/m²s | 134.886 |
    | OLR (iter 20 / 40) | 337.54 / 323.49 W/m² | 337.54 / 323.59 |
    | photosphere | 237.4 km, 368.30 K | identical |
    | mean albedo | 0.4988 | 0.4988 |
    | `Psi_max` (iter 40) | 105588.89 | 105588.91 |

    **(b) `c_u` is still identically zero, and the reason is geometry.** The parcel temperature at
    every active updraft cell is **302.7 K in BOTH arms, to the digit** — the geopotential cannot
    cool a parcel that never rises, and this one never does:

    - the dominant column has **cloud base at level 37 (236.4 km) and the LFS at level 38
      (256.0 km)** — 99 of the 173 convecting columns in the zonal slice, with base 38 / LFS 39 in
      another 13;
    - the updraft recurrence is `for(i = i_base+1; i <= i_LFS-1)`, which is **an empty loop when
      the two are adjacent**. The whole updraft is then whatever `initUpdraft` seeded at one level;
    - `|M_u| > coeff_recurr` (0.1 kg/m²s) holds in exactly **one level (i = 37) and 7 of 181
      columns** of that slice.

    So the parcel sits at 302.7 K with `q_v_u` = 508.6 g/kg at p = 0.083 bar, in an environment at
    265 K whose own `q_sat` is 0.186 — **the environment is condensing while the parcel is 38 K
    too warm to**. Item 53 predicted that adding `g·z` would let the updraft cool and condense.
    It would — over a hundred-kilometre ascent. There is no ascent.

    **Where the geopotential DOES act is the downdraft**, whose recurrence runs
    `for(i = i_lfs-1; i >= 0)` — from 256 km to the surface — so its DSE excess picks up a real
    `g·Δz` that does not cancel against the environment. That is the entire source of the `MC_t`
    change, and it is a change of PATTERN, not of scale: peak convective heating −81 %, peak
    convective cooling ×33. **A scheme whose updraft is one cell deep and whose downdraft spans
    the column is not a convection scheme yet**, and CLAUDE.md's "deep convection is inactive"
    (absolute-hPa triggers at 250 bar) is now measured rather than inferred.

    **The knob stays off.** Not because the arithmetic is in doubt — `cp·T` with no geopotential
    is simply not dry static energy — but because turning it on changes `MC_t` by 5× in a scheme
    that cannot presently produce an updraft, i.e. it improves a term nothing consumes. The
    ordering is: fix the trigger levels, then flip this.

    **(c) One name per species, and the `N2` collision is gone.** The VTK files carried `N2` for
    the **Brunt-Väisälä frequency squared** and `q_N2` for **nitrogen**, side by side in the same
    file. The array is now `brunt_N2` and the field `BruntVaisala_N2`. The eight species are
    written under bare names — `H2O CO2 N2 CH4 NH3 H2 CO SO2` — all mass fractions, all from
    `AtmMixture::split()`, the routine the thermodynamics and the radiation use. Before this,
    three conventions coexisted and one of them was wrong in substance: the field plotted as
    `CO2` was the **raw tracer array**, which carries no dilution, so it read 0.2053 everywhere
    while what the physics uses spans **0.1359–0.2805** in the zonal slice at iteration 40.
    Item 57's 320× carrier structure existed in the model and was invisible on screen for that
    reason. The transported array is still written, as **`CO2_tracer`** — that is where CO₂
    *transport* is visible, and the two fields answer different questions. The eight sum to
    1 − condensate (0.9807–1.0000, zonal slice). Physics-neutral, checked: the post-rename run
    reproduces OLR 323.49 W/m², `Psi_max` 105588.89 and photosphere 368.30 K exactly.

    **(d) A NULL A/B, because two runs of identical physics turned out not to agree.** Checking
    that the knob's off-branch reproduced the old binary produced the most useful measurement in
    this item. Four pairs, all *the same physics*, 40 iterations, moist physics from iteration 0:

    | pair | seed | first divergence | OLR @40 | water vapour @40 | `max M_u` | `Psi_max` |
    |---|---|---|---|---|---|---|
    | same binary, 24 threads, twice | run-to-run only | residuum #25, **7e-7** | identical | 1.2e-10 of scale | 0 | identical |
    | rebuild, no physics code touched (the ParaView rename) | none measurable | residuum #25 | **identical** | — | — | identical |
    | rebuild that re-arranged physics arithmetic (this item's edit, **knob OFF**) | codegen | residuum #3, **5e-8** | 323.85 → **323.49 (0.11 %)** | **795 of 7421 points, up to 64 g/kg** | **57 % of scale** | identical to 8 digits |
    | same binary, **24 vs 23 threads** | reduction order | residuum #2, **3.5e-7** | 323.49 → **325.14 (0.51 %)** | **1087 points, 12.2 % of scale** | **56 % of scale** | 1e-6 |

    **Three things follow, and the third is a warning about this file's own numbers.**

    **The moist and convective fields are the amplifier; the dynamics are not.** The same 5e-8
    seed leaves `Psi_max` identical to eight digits and `p_dyn` at 1e-7, while moving `M_u` by
    57 % of its own maximum and shifting 795 grid points of water vapour by up to 64 g/kg — all
    of them at levels 38–39 (256–277 km), the condensing band, where a trigger flipping on or
    off relocates a whole cell. The global max and min of water vapour still agree to seven
    digits: what changes is *where* it is, not how much.

    **A rebuild is not automatically neutral, and this is checkable rather than assumable.** The
    ParaView-only rebuild is indistinguishable from run-to-run noise; the header edit is not.
    The difference is whether the optimiser re-arranged arithmetic inside a physics loop under
    `-ffast-math -march=native`.

    **The 40-iteration OLR is reproducible to 0.11 % across a rebuild and only 0.51 % across a
    THREAD-COUNT CHANGE — and that is the number to hold against this file's OLR results.**
    Items 51, 52, 53, 59 and 60 reported their OLR effects as 0.02–0.11 %, i.e. **between one
    fifth and one twenty-fifth of the model's own 40-iteration reproducibility envelope**. Each
    was a same-binary, same-thread-count A/B, so each is valid *as a controlled comparison* and
    none is retracted; but none of them is a robust property of the model either, and none
    should be quoted as a physical effect without saying at what thread count it was measured.
    The pattern those items describe — "every microphysics correction dies at the albedo" — is
    unaffected and is if anything strengthened: the corrections are not merely small, they are
    small compared with re-running the same physics on 23 cores instead of 24. **`Psi_max`,
    albedo and the photosphere are unmoved in every pair** (Ψ to 1e-6 or better) and remain
    trustworthy to the digits printed; the OLR at 40 iterations is not.

    **And the honest note about this knob**: `ATM_MC_GEOPOTENTIAL`'s off-branch is
    arithmetically the old code — `+ phi_s(i)` with `phi_s ≡ 0` is exact, and the iteration-0
    min/max block confirms `s`, `s_u`, `MC_t` and `c_u` come out bit-identical after the first
    `MoistConvection` call — but **it is not bit-identical to the old *binary***, because the
    rebuild moved unrelated loops. This file's usual "off-branch bit-identical" standard is a
    statement about the source, and for the moist fields at 40 iterations that is not the same
    thing.

    **Reproduced at**: `python/run_gz_base` (off) and `python/run_gz_on` (on) for (a)–(b);
    `python/run_rename40` for (c); `python/run_gz_ctrl`, `python/run_null_same` and
    `python/run_null_23t` for (d).

62. **The startup energy-balance check was lit by 62 % of the model's own insolation, and the
    albedo in it was a literal.** Four lines of a printout, both defects in
    `cAtmosphereModel.cpp`:

    **(1) The planetary mean of a latitude parabola is not the average of its endpoints.**
    `short_wave_radiation[]` is `P(φ) = pole + (equator − pole)·(1 − (2φ/π)²)`, and a flux's
    planetary mean is its cos(latitude)-weighted mean, which integrates in closed form to
    `equator·8/π² + pole·(1 − 8/π²)`. At the shipped 298/0 that is **241.55 W/m² = S/4** for
    S = 0.71·1361 — exactly what `param.py` says the pair was *fitted* to deliver. The shipped
    `0.5·(equator + pole)` gave **149.0**, 61.7 % of it.

    **(2) The albedo was a bare `0.08` commented "molten surface"**, duplicating the
    `albedo_surface` parameter that exists for it — and inherited unchanged into ATHAD_COND,
    whose surface is a 240 °C **ocean**. `param.py`'s own complaint about `albedo_pole`/
    `albedo_equator` being "configuration theatre" applies in reverse: a literal beside a
    parameter of the same meaning is a second source of truth, and it is the one nobody edits.

    | | shipped | repaired |
    |---|---|---|
    | TOA mean | 149.00 W/m² | **241.55** |
    | absorbed + geothermal | 287.08 W/m² | **372.23** |
    | implied `t_skin` | 266.75 K | **284.64 K** |
    | reported gap to the configured 254 K | 12.75 K | **30.64 K** |

    **Diagnostic-only, and that is the whole of the good news.** `planetaryShortWave()` has
    always done the cos-weighting correctly, so the `t_skin` fixed point, the OLR and every
    result in this file are untouched — verified: the runtime line still reads `absorbed +
    geothermal 271.63 W/m² from albedo 0.4965`, i.e. 241.6 W/m² of insolation, as it always
    did. But this estimate is *the* line that judges whether the configured `t_skin` is
    consistent with the budget, and it was **17.9 K out** in the reassuring direction. The
    print now also states the insolation and the albedo it used, so the next reader can check
    it without reading the source.

    Found by asking whether `rad_equator_short`/`rad_pole_short` could be wrong "because the
    surface is water" in ATHAD_COND. They cannot — they are TOA fluxes, and the surface-absorbed
    error they are suspected of is one this file already fixed once (the 116/71 pair). **The
    suspicion was right about the neighbourhood and wrong about the parameter**, which is how
    the two defects above turned up. `albedo_surface` = 0.08 is still justified in `param.py`
    by "measured basaltic-melt albedos are 0.05–0.10" **in both trees**, and in ATHAD_COND that
    justification does not apply — open, and a physics decision rather than a repair.

63. **The convective triggers are fractions of surface pressure now — it changes nothing in
    ATHAD, and it is the reason ATHAD_COND's updraft condenses for the first time.**

    **The defect.** `p_stat_beg/end/Cloud_Base/base/diff/midlevel` were
    1000/970/900/800/200/700 **hPa**, absolute Earth surface pressures, inside every scan that
    locates the convective column. At 250 bar the whole column sits above 1000 hPa except its
    top ~70 km, so each scan ran in the **top** few levels instead of the bottom few. They are
    now fractions of the column's own surface pressure (`f_stat_* = Earth value / 1013.25`,
    times `p_stat.x[0][j][k]`), which is what they always meant and which reproduces Earth's
    geometry exactly on Earth.

    **ATHAD_COND already had this repair** — made when the fork was cut, unconditional, with a
    comment explaining that at 60 bar the broken triggers do not merely fail, they *fire in the
    wrong place*. ATHAD was the straggler, carrying the defect with a note saying so. Ported
    verbatim rather than re-derived, so the two files now agree line for line.

    **MEASURED IN ATHAD — 40 iterations, and it is a null:**

    | | absolute hPa | fractions of p_surf |
    |---|---|---|
    | dominant cloud base / LFS | level 37 / 38, **99 of 173 columns** | **identical** |
    | second most common | 37–38 / 38–39, 39 columns | **identical** |
    | cells with `\|M_u\|` > `coeff_recurr` | 7 of 7421 | **7** |
    | `max c_u` | 0.000000 | **0.000000** |
    | OLR @40 | 323.49 W/m² | 323.75 (+0.08 %, inside the item-61 envelope) |
    | `Psi_max` | 105588.89 | 105588.86 |

    **WHY IT IS A NULL, AND WHY THAT IS THE RESULT: INVARIANT 2 PINS THE CLOUD BASE, NOT THE
    GATE.** The cloud-base test is `q_v_u >= scale·q_sat`, and this column's own startup profile
    prints `condensation possible from i = 37 (236.4 km) upward` — water is supercritical below
    ~177 km and sub-saturated to 236 km, so **no threshold setting can put a cloud base lower**.
    The pressure gate was a second lock on a door the thermodynamics had already bolted.

    **Item 61's closing claim is corrected.** It ended "the blocker is the trigger levels, not
    the thermodynamics" and made fixing them the top microphysics job. Fixed, they move nothing:
    ATHAD's convective layer is one cell deep because its **condensing** layer is four cells
    deep and the LNB sits one level above the base. That is a property of a 250 bar supercritical
    column, not a defect.

    **THE SAME REPAIR IS LOAD-BEARING NEXT DOOR, WHICH IS HOW THE NULL BECAME INFORMATIVE.**
    Same code, same fractions, opposite outcome, because the regime differs:

    | | ATHAD (250 bar, supercritical) | ATHAD_COND (60 bar, 513 K sea) |
    |---|---|---|
    | cloud base | 236.4 km | **1786–2200 m**, at the sea |
    | LFS | 256.0 km — one level up | **20 517–28 130 m** |
    | convecting cells, zonal slice | 7 of 7421 | **5218 of 11041** |
    | `max \|M_u\|` | 0.135 kg/m²s | **3.0000 — pinned at `M_max`** |

    **AND OVER THAT 20 km COLUMN, `g·z` FINALLY DOES SOMETHING.** Item 52 found the geopotential
    missing from the static energy, item 53 predicted that adding it would let the updraft cool
    and condense, and item 61 added it and got nothing — in ATHAD, where there is no ascent.
    Run in ATHAD_COND (`ATM_MC_GEOPOTENTIAL=1`, 40 iterations, one binary, env-only A/B):

    | | off | on |
    |---|---|---|
    | **`max c_u`** | **0.000000 g/kg/s** | **7.19e-04 @ 24 854 m, 37°S** |
    | `max s_u`, and where | 2.5460 at **661 m** | 2.9111 at **24 854 m** |
    | `max MC_t` | 1.030e-03 K/s | 1.175e-03 (+14 %) |
    | `max q_c_u` | 7.0537 g/kg at 1017 m | 7.1392 at 1785 m |
    | OLR, photosphere, `Psi_max`, precipitable water | 273.70, 65.9 km, 40515.96, 160626.497 | **all identical** |

    **The updraft condenses for the first time in either fork**, and it does so 25 km up, which
    is where a parcel that has cooled at g/cp_l = 4.81 K/km for 23 km should first saturate. The
    `max s_u` location moving from the cloud base to the top of the column is the same statement
    read off the other field: with the geopotential in, static energy is largest where `g·z` is
    largest, as a static energy must be. **Nothing integrated moves** — OLR, photosphere, Ψ and
    the water column are identical to every printed digit — which is the albedo wall again, and
    expected: the deck is already saturated at `albedo_cloud`.

    **What this closes and what it opens.** Closed: the three-item chain 52 → 53 → 61 about the
    geopotential now has its demonstration, in the fork whose geometry allows one. Open: `M_u`
    is **pinned at `M_max` = 3.0 kg/(m²s)** in ATHAD_COND, so half that model's convective mass
    flux is a cap rather than a result — the next Earth constant to size, and the third time this
    file has found a cap standing where a measurement should be.

    **VTK bookkeeping in the same pass, both trees:** the six `Ubud_*` dumps per slice are
    commented out (a third of the file for a diagnostic read as a min/max in `Results_Atm`,
    which still prints it), and `Q_Latent`/`Q_Sensible` are written in all three writers —
    ATHAD_COND's **zonal** pair was still commented while its radial and longal were live, so
    the slice most looked at was the one missing them. Both trees now emit 84 scalar fields per
    zonal slice, verified against freshly written files.

64. **Item 16 instrumented: the saturation adjustment works, the supersaturation is re-created
    between calls, and freeing the profile makes it worse — my invariant-3 explanation is
    retracted in the same item that produced it.**

    Two print-only/default-off knobs. `ATM_SAT_TRACE=1` prints the Newton loop for one cell at
    levels 37/38; `ATM_SAT_NO_ALPHA=1` sets `alpha_entry` to 1.

    **(a) The adjustment is not stalling.** At level 38 (256 km), entry `q_v` = 0.6724 against
    `q_sat` = 0.0883 — **7.62× supersaturated** — the loop condenses 0.045 kg/kg on pass 1, the
    latent heat throws `T` from 266 to 342 K, past boiling at 0.02 bar so `q_sat` becomes
    **1.000000**, the target inverts and passes 2–9 evaporate it back. From pass 10 it settles
    into a **limit cycle**: `q_v` pinned at 0.6628 while `q_sat` oscillates 0.64 ↔ 0.67. It exits
    having condensed **0.47 kg/kg**, leaving the cell at `q_v/q_sat` ≈ 1.005. **One call reaches
    ~100.5 % RH.** The earlier guess that heavy damping stalls it far from saturation is wrong.

    **(b) But the exit test measures the STEP, not the residual.** `|q_v_b/q_v_hyp − 1| ≤ 1e-6`
    falls to 4.4e-06 by pass 20 while the residual still swings ±0.003. Here it happens to exit
    near the right answer; it would exit just as confidently anywhere on that cycle.

    **(c) THE PRESCRIBED PROFILE IS HOLDING THE SUPERSATURATION DOWN, NOT CAUSING IT.** I wrote
    that `densities()` re-imposing the adiabat discards the adjustment's 22 K of latent heating
    and re-creates the supersaturation each iteration. Measured over five iterations, entry
    `q_v/q_sat` at level 38:

    | call | prescribed (default) | `ATM_PROGNOSTIC_T=1` |
    |---|---|---|
    | 1 | 7.62 | 7.62 |
    | 2 | 16.50 | **975.97** |
    | 3 | 10.82 | **1933.97** |

    **Retracted.** Freed, the cell cools to **220–228 K** and `q_sat` collapses 0.088 → 0.00035,
    while `q_v` is flat at 0.683 to six digits. **The supersaturation is `q_sat` falling, not
    `q_v` rising** — which is why it reads as a moisture problem and is not one. Item 16 should
    be re-pointed at whatever cools the free-running top, the same unexplained behaviour as
    items 45–47.

    **(d) `alpha_entry` is a genuine defect and NOT the lever.** The −37 °C ice threshold is
    applied twice: inside the loop as the phase split (`CND`/`DEP`, correct — below −37 °C
    condensation becomes deposition), and again as `alpha_entry`, a **master gain on all five
    write-backs** (`S_c_c`, `c`, `cloud`, `ice`, `t`). A cell below ~236 K keeps a few per cent
    of what the loop computed, deposition included; below 213.2 K the entry test skips it
    outright. Invisible on Earth (~0.1 g/kg at −37 °C), live here at **683 g/kg with a gain of
    0.036**. But with `ATM_SAT_NO_ALPHA=1` the prognostic ratios go 7.62 → **2094 → 2977**,
    *worse* than 976 → 1934: the gain was masking the cooling, not causing the excess. Repair it
    because a threshold applied twice is a defect, not because it will move item 16.

    **Method note.** Both wrong explanations in this item were mine, both were stated with
    confidence, and both died to the same three printed columns. The trace cost one afternoon
    and is print-only; it should have existed before the first explanation, not after the second.

65. **What cools the free-running top is, in large part, the under-converged radiation solver —
    which means every prognostic-column number in this file rests on it.**

    Item 64 left the supersaturation at levels 37/38 attributed to a falling `q_sat` and asked
    what cools the cell. Same five-iteration trace, same cell, **only the solver changed**:

    | call | default (`ATM_N_LAMBDA` = 4) | `ATM_RAD_DIRECT=1` |
    |---|---|---|
    | 1 | 266.22 K | 266.22 K |
    | 2 | **219.72** | **228.92** |
    | 3 | **228.25** | **243.25** |
    | `q_v/q_sat` at call 3 | 1934 | **97.9** |

    The exact closed-form solve **halves the cooling and cuts the supersaturation ratio 20×**,
    and the trajectory turns around — 229 → 243 recovering, against 220 → 228 still wandering.
    Reproduce with
    `ATM_PROGNOSTIC_T=1 ATM_RAD_DIRECT=1 ATM_SAT_TRACE=1 cli/had c5.xml` (nm = 5,
    `moist_phys_start_iter` = 0).

    **This joins two threads that were being worked separately.** Item 30 established that four
    Lambda sweeps under-converge a 250 bar column by 4× and left `ATM_RAD_DIRECT` default-off
    pending one long measurement. Items 45–47 then measured the PROGNOSTIC column — 2927 W/m²,
    "still falling at 400 iterations" — with that same under-converged solver underneath. Item
    64's cooling top is the same defect seen from the humidity side. **Every prognostic-column
    number in this file was produced by a solver known since item 30 to be 4× under-converged**,
    which is a reason to redo them rather than to extend them.

    **Not the whole story**: with the converged solver the top still falls 266 → 229 K in one
    call, so radiation is a first-order part and not the only part. What remains is the open
    question, and it is now a smaller one.

    **The case for flipping `ATM_RAD_DIRECT` on by default has a second, independent argument
    now** — it is exact rather than under-iterated, 10× cheaper, changes the standard prescribed
    configuration by 0.15 %, and halves a cooling that has been corrupting the prognostic branch.
    What it still lacks is the long run item 30 asked for.

66. **Item 30's long measurement, finally run: the prognostic column CONVERGES, at 224 W/m²
    against item 45's 2927, and "10.8× absorbed and still falling" was a solver artefact.**

    400 iterations, 24 threads, `ATM_PROGNOSTIC_T=1 ATM_RAD_DIRECT=1`, otherwise the shipped
    configuration (`im` = 41, `moist_phys_start_iter` = 300) so it is directly comparable with
    items 45–47. OLR every 20 iterations:

    ```
    dry   (20-300):  474.45  363.26  303.35  271.01  253.17  242.98  236.89  233.05
                     230.47  228.63  227.25  226.15  225.25  224.48  223.83
    moist (320-400): 1021.93  766.57  751.00  1036.27  782.34
    ```

    **THE DRY PROGNOSTIC COLUMN CONVERGES.** Successive differences are 111, 60, 32, 18, 10, 6,
    4, 2.6, 1.8, 1.4, 1.1, 0.9, 0.8, 0.65 — a decelerating approach to **~223.8 W/m²**, not a
    trend. Against absorbed SW + geothermal of **271.57 W/m²** that is an imbalance of
    **−47.7 W/m²**, i.e. the column emits 82 % of what it takes in. It is the first prognostic
    number in this file that is the same order as the energy input.

    **ITEM 45 IS RETRACTED ON ITS HEADLINE.** With `n_lambda = 4` the same branch gave
    15 764 → 8 542 → 6 211 → 2 927 W/m², "10.8× the 271 absorbed and still falling at 400".
    That was the under-converged Lambda iteration, exactly as item 65 inferred from three trace
    lines. **What survives from item 45 is its mechanism**: `emission from isothermal skin` is
    **0.0 % at every diagnostic here too**, so the prognostic branch genuinely escapes the
    `t_skin` pin — and now it escapes it to a converged value rather than to a runaway.

    **`t_skin` still does not move** (262.99 → 263.01 K, target 263.07 from absorbed 271.57),
    because it is the fixed point of the absorbed flux and independent of the column. The
    difference is that the OLR is no longer σT_skin⁴ = 271.2: it is 223.8 and it got there by
    integrating. **For the first time the reported OLR is neither pinned to `t_skin` nor
    diverging.**

    **THE MOIST ONSET IS A SHOCK AND DOES NOT SETTLE IN 100 ITERATIONS.** At iteration 300 the
    moist physics switches on and the OLR jumps 223.83 → 1021.93, then oscillates 766.57,
    751.00, 1036.27, 782.34 — swings of ±35 % with no sign of a limit. **The dry column is
    converged; the moist one is not**, and nothing in the moist branch here should be quoted.
    That the shock is 4.6× the dry value is itself the finding: turning on condensation in a
    column that had settled dry is not a perturbation.

    Photosphere over the run: **237.4 → 255.3 km**, `T_ph` **388.8 → 342.3 K**, climbing and
    cooling steadily through both phases.

    **This is the third argument for making `ATM_RAD_DIRECT` the default**, after "it is exact
    rather than under-iterated" and "it is 10× cheaper": the branch it was gating turns out to
    converge, and every prognostic number measured without it — items 45, 46, 47 — was measured
    on a solver known since item 30 to be 4× wrong. Redo them, do not extend them.

67. **THE t_skin/OLR PROBLEM IS A MISSING FACTOR OF 2. `t_skin` holds the planet's EFFECTIVE
    EMISSION temperature where the top layer's SKIN temperature belongs, so the model assigns
    its transparent lid the whole planetary energy input as its emission — and "the energy
    balance closes" becomes an identity.**

    Both sites that set `t_skin` — the startup estimate in the composition printout and the
    fixed-point relaxation in `updateSkinTemperature()` — solve

        sigma * t_skin^4 = absorbed SW + geothermal = F

    That is the definition of `T_eff`, the temperature a black body needs to radiate the whole
    budget. It is not the temperature of the top of the column. An optically thin top sees no
    downward flux, so it absorbs only the upward stream and re-emits half up and half down:

        sigma * T_skin^4 = F / 2   ->   T_skin = T_eff / 2^(1/4)

    At this budget that is **263.07 K against 221.22 K — 41.9 K too warm, and a factor of two
    in the emission.**

    **THE MODEL ALREADY CONTAINS THE RIGHT RELATION, TWICE, IN THE FILE THAT THE PRESCRIPTION
    OVERWRITES.** `MultiLayerRadiation`'s flux sweep reduces at the top, where `dn -> 0`, to
    `sigma*T^4 = up/2`, and its comment says in as many words: *"the classical skin
    temperature, which this model has until now been PRESCRIBING as t_skin"*. The direct
    solver's validated `eps -> 0` limit is `T_i = T_s/2^(1/4)`, the same factor. So the
    radiation solver and the prescription that overwrites its answer every iteration disagree
    by 2^(1/4), and the prescription wins.

    **It is exact here, not an approximation being quibbled over.** All shortwave is deposited
    at the surface — `SW_abs` appears in exactly one place, `MultiLayerRadiation.h`'s surface
    energy balance, and nowhere in the column — so the atmosphere is a pure grey long-wave
    medium over a shortwave-heated surface. That is the classical skin problem, where
    `sigma*T_skin^4 = F/2` is the exact result rather than a fit.

    **WHY IT IS THE t_skin/OLR PROBLEM AND NOT A 42 K ERROR IN A CORNER.** `densities()` sets
    the whole upper column to `max(t_skin, T_ad)`, so the lid is exactly `t_skin`. Measured at
    the shipped configuration, 60 iterations, 24 threads:

    | | shipped | `ATM_SKIN_GREY=1` |
    |---|---|---|
    | lid emissivity | 0.0068 | 0.0066 |
    | `t_skin` | 262.95 K | **221.12 K** |
    | `sigma*T_lid^4` | **271.10 W/m²** | 135.55 W/m² |
    | absorbed SW + geothermal | **271.10 W/m²** | 271.10 W/m² |

    **A layer of emissivity 0.0068 is prescribed to emit the planet's entire energy input, to
    six figures, at every diagnostic.** That is the defect in one line. And since the OLR
    descends onto `sigma*t_skin^4`, choosing `t_skin` so that `sigma*t_skin^4 = F` makes
    "OLR -> absorbed" **arithmetic**. Item 25's −1.26 W/m² residual and item 43's "the
    imbalance is identically the OLR's distance from a constant" are the same event seen
    twice: the constant is `F`, and it is `F` because the lid was assigned the temperature at
    which the budget closes.

    **THE OLR MOVES WITH IT, AND THE IMBALANCE CHANGES SIGN.** Same binary, same 24 threads,
    the only difference the factor of 2:

    ```
    iteration        20        40        60
    shipped OLR   333.88    305.89    291.84     imbalance  -62.81  -34.80  -20.74
    grey    OLR   295.10    243.43    212.45     imbalance  -24.03  +27.66  +58.66
    ```

    −27.2 % at iteration 60, and the reported imbalance goes from −20.74 to **+58.66 W/m²**.
    **This is the first FORCING that moves the converged prescribed OLR.** It has failed to
    respond to kappa at 64x (0.10 %), the circulation at 500x (0.03 %), `im`, and the grid
    stretch at 6x (0.36 %) — items 29, 27, 39, 41. (Structural repairs have moved it before —
    item 22's initialisation order halved it — but nothing that was turned as a knob.) The
    lever was never opacity or resolution. It was the number the profile is clamped to.

    **AND IT IS NOT THE MECHANISM ITEMS 25 AND 43 BLAMED.** Those attributed the pin to the
    photosphere migrating into the isothermal skin. Here the photosphere does not move —
    237.4 → 237.3 km, `T_ph` 370.50 → 369.84 K — and `emission from isothermal skin` is
    **1.1 % of columns in BOTH arms**. So the OLR's sensitivity to `t_skin` is carried by the
    optically thin layers *above* the photosphere, which are prescribed too, not by the
    emission level sitting in the skin. Item 43's rising skin fraction is real; it is not what
    makes the OLR follow `t_skin`.

    **THE HISTORY IS IN THE REPOSITORY'S OWN COMMENTS.** `param.py` records that an earlier
    `t_skin` of 231.3 K came from `T_skin = (OLR/2sigma)^(1/4)` and was replaced because the
    OLR it read was a model output — circular, and correctly removed. **The circular part was
    reading a measured OLR. The factor of 2 was never the circular part, and it went out with
    it.** The repair is to keep the budget as the source of `F` and put the `/2` back.

    **WHAT THIS DOES NOT DO.** It does not free the OLR: the top is still prescribed, just at
    the right value, and the OLR still descends onto `sigma*t_skin^4`. It makes the closure
    *falsifiable* rather than automatic — and the honest consequence is that the corrected
    model is **much further out of balance**, which is invariant 3's warning arriving again:
    the better-looking of the two numbers was the more assumed one. Where the corrected OLR
    lands is not claimed here; at 60 iterations it is 212.45 and still falling, and this file
    has twice been caught extrapolating a monotone radiation trend (items 30 and 45/66).

    `ATM_SKIN_GREY` applies the factor, in `skinTargetFromFlux()`, covering both sites.
    **DEFAULT ON as of 2026-08-21**; `ATM_SKIN_GREY=0` restores the old branch, and that branch
    is a verified null — every printed radiative diagnostic identical to the pre-change
    binary's, though the rebuild reproduces item 61's ~3e-7 residuum divergence.

    **Why this one is flipped on a single pair, when this repo's convention is measure-then-
    wait.** The convention exists for knobs that are choices — `ATM_RAD_DIRECT`, `mc_M_max`,
    `n_cells_hemisphere` — where the measurement decides which value is better. This is not one
    of those. `sigma*T^4 = F` for a layer of emissivity 0.0068 is the wrong equation, and the
    right one is written in `MultiLayerRadiation`'s own comments, twice. There is no
    configuration in which the old branch is the physics. The pair was run to find out what the
    correction *costs*, not whether to make it.

    **WHAT IT COSTS, STATED PLAINLY.** The model is now far out of balance: 271 W/m² absorbed
    against 147 emitted, imbalance +123.76 W/m² and widening, where before it closed to −1.41.
    **Every OLR and imbalance figure recorded in this README and in CLAUDE.md before this date
    was measured with the knob off.** None of them is retracted — they are correct measurements
    of a branch still reachable with `ATM_SKIN_GREY=0` — but none of them describes the shipped
    model any more. The paired series above is the one that does.

    **And the `t_skin` = 254.0 default is still derived by the old formula.** The consistent
    start is 213.6 K. It is left alone deliberately: 254.0 is what the measured pair above began
    from, and changing it would mean the shipped default is described by no measurement in this
    file. It is a fixed-point *start* that relaxes to 221.12 within a handful of radiation
    calls, so it costs a short transient and nothing else — but it is an open inconsistency and
    the next person to touch `t_skin` should close it and re-measure together.

    **THE FOLLOW-UP THIS OPENS, AND IT IS NOW CHEAP.** An isothermal top is not the grey
    radiative-equilibrium solution either — that is `sigma*T^4 = (F/2)(1 + 3*tau/2)`, of which
    `F/2` is only the `tau -> 0` end. The prescription `max(t_skin, T_ad)` is the first term of
    a profile whose second term the model can now evaluate, because item 42 made `tau_above` a
    real array. Replacing the clamp with `max(T_rad(tau_above), T_ad)` costs one line, keeps
    the adiabat wherever it is warmer, and gives the prescribed branch a top that is *right*
    rather than merely *scaled correctly* — and unlike `ATM_PROGNOSTIC_T` it does not require
    the column to integrate its way there. That is the cheapest remaining move on invariant 3.

    ### The 200-iteration pair (same binary `70d6f1fe`, 24 threads, dry throughout)

    ```
    iteration      20      40      60      80     100     120     140     160     180     200
    shipped OLR  333.88  305.89  291.84  284.16  279.67  276.90  275.11  273.91  273.10  272.53
      imbalance  -62.81  -34.80  -20.74  -13.06   -8.56   -5.78   -3.99   -2.79   -1.97   -1.41
      skin %        1.1     1.1     1.1     9.9    29.8    42.0    50.8    53.0    53.0    53.0
    grey    OLR  295.10  243.43  212.45  192.21  178.20  168.10  160.62  154.98  150.68  147.37
      imbalance  -24.03  +27.66  +58.66  +78.90  +92.91 +103.02 +110.50 +116.14 +120.45 +123.76
      skin %        1.1     1.1     1.1     1.1     1.1     1.1     1.1     1.1     1.1     6.6
    ```

    **The shipped arm reproduces the record**: −1.41 W/m² at 200 against item 25's −1.26, so
    this is the same measurement, not a new configuration.

    **BOTH ARMS DESCEND ONTO `sigma*t_skin^4`, AND THAT IS MEASURED RATHER THAN EXTRAPOLATED.**
    The two constants are 271.11 and 135.55 W/m² — ratio 0.5000, the factor, exactly. Tracking
    the OLR's *distance* from its own arm's constant removes the need to guess a limit:

    ```
    excess over sigma*t_skin^4
      shipped   62.77  34.78  20.73  13.05   8.56   5.79   4.00   2.80   1.99   1.42
      grey     159.55 107.88  76.90  56.66  42.65  32.55  25.07  19.43  15.13  11.82
    decay ratio, successive
      shipped    0.554  0.596  0.630  0.656  0.676  0.691  0.700  0.711  0.714
      grey       0.676  0.713  0.737  0.753  0.763  0.770  0.775  0.779  0.781
    ```

    A clean geometric decay to zero in both, ratio settling at ~0.71 and ~0.78. The shipped arm
    is five e-foldings further along the identical curve and has demonstrably arrived. **So the
    OLR is `sigma*t_skin^4` in both arms; the shipped model closes its budget for the single
    reason that `sigma*t_skin^4` was set equal to the absorbed flux.** The grey arm is NOT at
    its limit at 200 — 11.82 W/m² of excess remains — and no limit is claimed for it.

    **THE SKIN FRACTION SETTLES THE MECHANISM, AND IT SETTLES IT AGAINST ITEMS 25 AND 43.** In
    the shipped arm it climbs 1.1 → 53.0 % and **saturates** — flat over the last three
    diagnostics, which retires item 43's "39.2 % is not a saturation". In the grey arm it stays
    at **1.1 % through 180 iterations** and reaches only 6.6 % at 200. **Yet the grey arm's OLR
    tracks its own `sigma*t_skin^4` just as tightly.** The emission level migrating into the
    isothermal skin is therefore not what pins the OLR: it is a *consequence* of a warm lid —
    a lid at 262.96 K sits within 1 K of the profile the photosphere climbs through, a lid at
    221.12 K does not — and the pin operates identically without it. Items 25 and 43 measured a
    real correlation and read the causality backwards.

    **WHAT THE CORRECTED ARM SAYS, PLAINLY.** The planet absorbs 271 W/m² and emits 147, and the
    gap is *widening* toward ~136. A prescribed adiabat with a correct skin is far too cold at
    the top to balance this budget, and it says so instead of hiding it. That is the finding:
    not that the corrected number is better-looking — it is much worse — but that it is now a
    measurement of the prescription rather than a restatement of it.
68. **The meridional streamfunction does not close at the ground, the violation carries TWICE
    the mass flux of the circulation it hides, and it is written by the INITIAL VELOCITY
    PROFILE — not by the dynamics, not by the solver, and not by `exp_rm`.**

    `Psi` is a mass streamfunction, so it must vanish at the ground: no mass flows through the
    surface. Both boundaries were checked rather than assumed — the radial velocity `u` at
    `i = 0` is **identically 0.0000 m/s at every latitude**, and `|Psi|` at the top three levels
    is ~0 — so the column-integrated meridional mass flux across every latitude circle is
    forced to be zero. It is not:

        Psi at the lid ......................... 0            (by construction)
        Psi at the ground ...................... 9.891e+13 kg/s   <- must be 0
        max |Psi| above 20 km (the real cells).. 4.733e+13 kg/s
        ratio .................................. 2.09

    So the field ParaView draws is dominated 2:1 by a continuity violation, peaked sharply at
    +-5 deg latitude. **That is why the cells never looked like they hand off to one another:
    there were no closed cells to see.** It decays at 0.04 %/iteration — 1.0636e14 to 9.891e13
    over 180 iterations — i.e. it is not a transient the run works off.

    **IT IS THERE BEFORE THE TIME LOOP RUNS.** A 4-iteration probe:

        iteration 1 ...... Psi(ground) 1.0703e+14   ratio 2.343
        iteration 40 ..... Psi(ground) 1.0567e+14   ratio 2.313

    **98.7 % of the final non-closure is present at iteration 1.** The dynamics add essentially
    nothing.

    **THREE NULLS, EACH KILLING A PLAUSIBLE STORY.** Every arm below is byte-identical apart
    from `output_path`; the knobs are environment variables.

    - **The time-loop pressure solver is not the cause.** `ATM_PRESS_SWEEPS` 1/4/16/64 at 40
      iterations moves the RMS of `Psi(ground)` from 2.8409e13 to 2.8502e13 — **+0.3 %, in the
      wrong direction.** The knob is live (output is not bit-identical and runtime rises), so
      this is a real null. Item 30's under-converged-solver shape does not repeat here.
      Incidentally the pressure solve is only **0.104 % of the step cost per sweep**: sweeps
      are nearly free, which is worth knowing on its own.
    - **`exp_rm` is not the cause, and this is the first measured consequence it demonstrably
      does NOT have.** A `zeta` scan sweeping the metric error from 13.4x to **exactly 1.00x**
      (`zeta = ln(1.5)`, where the quadratic-stretch Jacobian is correct by construction):
      RMS 2.841 / 3.011 / 3.029 / 2.966 e13 at `zeta` 3.0 / 2.0 / 1.5 / 0.405. No trend, and
      the arm where the metric is exact is **the worst of the four**. `exp_rm` remains live and
      unrepaired elsewhere (item 39); it is simply not what breaks continuity.
    - **The cell compression is not the cause.** `cell_lat_scale` 0.33 vs 1.0 crossed with
      `cell_amp_mode` 0 vs 1 moves absolute `Psi` by 4.5x while the closure ratio stays at
      **2.08-2.40**. The non-closure scales with the circulation instead of being caused by the
      compression, so it is in the profile's SHAPE at any amplitude.

    **THE LEVER IS THE INITIAL PROJECTION, AND IT WAS NEVER SCANNED BECAUSE IT IS A DIFFERENT
    KNOB.** `project_initial_velocity` exists precisely to remove "the unphysical dilatational
    artefact of the analytical profile", runs 200 passes x **1** sweep, and was deliberately
    decoupled from `ATM_PRESS_SWEEPS` (item 5c9a516) so that knob varies the time loop alone.
    Measured at iteration 1:

        ATM_PROJ_SWEEPS      RMS Psi(ground)     vs default
              1 (shipped)      2.910e+13            -
             10               1.383e+13          -52.5 %
            100               1.284e+13          -55.9 %

    **It plateaus rather than closing.** Ten sweeps buy 52 %, a further 10x buys 3.4 more, so
    ~44 % of the non-closure is **structural and outside the projection's reach entirely**.
    That residual is unexplained. The likeliest candidate — that the projection enforces
    `div(rho_bar u) = 0` on a BASE-STATE density while `Psi` integrates the actual `r_humid` —
    is a hypothesis and nothing more; five hypotheses died in the course of this item and it
    has not been tested.

    **`ATM_PROJ_SWEEPS` DEFAULTS TO 10 SINCE 2026-08-22, AND `Psi_max` NO LONGER MEANS WHAT IT
    MEANT.** 10 is the knee of the curve above, not a converged value; `ATM_PROJ_SWEEPS=1`
    restores the old branch exactly. Measured at 40 iterations, 24 threads:

        RMS Psi(ground) ........ 2.841e13 -> 1.325e13    -53.3 %
        max |Psi| above 20 km .. 4.569e13 -> 5.222e13    +14.3 %   the real circulation
        max |Psi| ANYWHERE ..... 1.057e14 -> 5.222e13    -50.6 %
        closure ratio .......... 0.622    -> 0.254       -59.2 %

    **The two maxima COINCIDE on the new default and did not on the old one.** The global
    maximum of `Psi` has moved off the ground and into the interior — so at one sweep the
    largest value in the streamfunction was the spurious surface flux itself, and **`Psi_max`
    was reporting the defect rather than the circulation.** That is the whole reason `Psi_max`
    falls 50 % while the circulation it is supposed to measure rises 14 %. Every Ψ figure in
    this file and in CLAUDE.md predates the change; do not compare across it. **The OLR is
    unmoved — 243.43 W/m² in both arms** — which is item 27's null again, not a surprise.

    **Grid-scale noise is ruled out** as a contributor: the 2-delta oscillation index of `u`,
    `v` and `Psi` is **0.006-0.2** across the column where a checkerboard would be ~4. Odd-even
    decoupling on the collocated grid would have explained the residual neatly. It is not there.

    **TWO METHOD NOTES, BOTH EARNED THE HARD WAY.**

    - **`Psi` is the trustworthy witness and it was already in the output.** It integrates with
      the TRUE layer thickness (`get_layer_height(i+1) - get_layer_height(i)`) and the real
      `r_humid`, touches `exp_rm` nowhere, is antisymmetric to **0.01 %** (invariant 1), and is
      **thread-count independent** — 8 vs 24 threads reproduce it to every printed digit, unlike
      the OLR, which item 61 found moving 0.51 % across a thread change. Every number in this
      item came from re-reading VTK files that had been sitting in `output_Hadean/` unexamined.
    - **A max is the wrong norm.** `Psi(ground)` was first measured as a max over latitude, and
      the max sits at +-5 deg inside the Hadley cell. A change localised to other cells then
      reads as *bit-identical to nine significant figures* while the field underneath it moves
      by 55 %. The RMS over latitude is the norm this quantity needs. Both nulls above were
      re-verified with it and both survive.

69. **Four of the five prescribed circulation cells turn the SAME WAY, because ATOM_Precipitation
    only has three cell templates and Earth's polar cell already carries the Ferrel's sense.**

    `VelocityInitializer::centreAmp` gives cell 0 the Hadley template, cell `n-1` the polar one,
    and every cell between them a **Ferrel copy**. At the shipped `n_cells_hemisphere = 5` the
    prescribed meridional amplitudes are

        cell     0        1        2        3        4
        v_trop  -3.0    +4.0    +4.0    +4.0    +0.5
                Hadley  Ferrel  Ferrel  Ferrel  polar

    Only cell 0 opposes the rest. The alternation cannot simply be read off the inherited table
    either, and `6afb2a1` said so at the time: **Earth's polar template is `+0.5`, the Ferrel's
    sign, not the Hadley's** — so the three-cell layout this was ported from is not a
    direct/indirect/direct sequence to begin with. The port is incomplete for the plainest of
    reasons: *"the formula only supplies edges Earth has no value for."*

    Measured in the field rather than inferred from the table: at iteration 1 the northern
    hemisphere is single-signed over almost its whole depth, with **3 sign bands at the surface,
    not 5**, the rest being 1-5 deg slivers at cell edges.

    **THE SAME FILE ALREADY ASSUMES THE ALTERNATION IT DOES NOT IMPOSE.** `edgeRadialCoeff`
    returns `(k % 2 == 0) ? mag : -mag` — the ascent and descent branches at the cell EDGES
    alternate by parity, which is only correct if the CORES alternate too. Shipped, `u`
    alternates and `v` does not, so the two components disagree about how many cells there are.

    **`ATM_CELL_ALTERNATE=1` imposes the parity**, matching `edgeRadialCoeff`: even `k` direct,
    odd `k` indirect. Each template keeps its own magnitude; only the sense is imposed, by
    flipping `v_trop` and `v_surf` together so the cell's vertical structure survives and only
    its rotation reverses. **`w` is deliberately NOT flipped** — the zonal wind is the jet
    structure, set by thermal-wind balance rather than by the overturning sense, and flipping it
    would confound the measurement. Default OFF; the off-branch is **bit-identical** at all four
    diagnostics.

    **IT SURVIVES.** Counting only bands at least 8 deg wide, so edge slivers cannot inflate the
    count, northern hemisphere:

                     9km  13km  23km  36km  55km  79km
        OFF  it10      3     3     3     1     1     2
        OFF  it40      3     3     2     1     1     1
        ON   it10      4     4     4     4     4     5
        ON   it40      4     4     4     4     4     3

    The shipped arm collapses to a **single hemispheric overturning above ~23 km** and keeps
    eroding. The alternating arm holds **four counter-rotating bands from 9 to 55 km, stable
    across all four diagnostics**. Band widths at 36.5 km: `+84` (one band, the whole
    hemisphere) against `-14 +18 -20 +29`. `max |Psi|` is **1.057e+14 in both arms** — this is
    a change in organisation, not in strength.

    **IT IS NOT THE CONTINUITY FIX.** Mean `|Psi(ground)|` improves 7 % and the RMS 0.66 %,
    because the dominant non-closure sits at 10 deg N in cell 0, whose sense was already
    correct. Item 68's residual is untouched.

    **What this does and does not establish.** It establishes that the co-rotating templates
    were why the cells never turned coherently, and that the corrected state is self-consistent
    rather than being torn apart the way the co-rotating one is. It does **not** establish that
    the model would maintain these cells: 40 iterations is minutes of physical time (item 47),
    and nothing here supplies the baroclinic eddy momentum flux that drives an indirect cell.
    The honest claim is that the initial state is now coherent, which it demonstrably was not.

70. **The zonal ParaView writer had three defects, and the worst of them has been distorting
    every latitude-height figure this project has produced.**

    All three are output-only — no physics reads these fields — but every zonal plot in this
    repository's history was made through them.

    - **`inv_u_0` where `u_0` belongs.** The arrays hold non-dimensional `u/u_0`, so m/s is
      `* u_0`, which is what `dump_zonal` writes for the u/v/w SCALARS three lines above. The
      VECTORS field divided instead, giving `u/u_0^2`. Measured in a shipped file: the scalar
      `v-Component` and the vector's y-component differ by **exactly 64 = u_0^2** at
      `u_0 = 8 m/s`. A uniform factor cannot rotate an arrow, so this never moved a direction —
      it made every magnitude 64x too small and inconsistent with the scalars plotted beside it.
      Present in `u-v-Cell`, `u-v-Updraft` and `u-v-Downdraft`.
    - **Velocities in m/s drawn on an index-space geometry.** The points were written at
      `x = i*dx` per LEVEL and `y = j*dy` per LATITUDE INDEX, while the vector was raw m/s. The
      two axes are compressed by wildly different factors, so glyph angles had no relation to
      the geometry they were drawn on. `uv_plot` is the same field converted into plot units per
      day, so glyphs and the stream tracer share one coordinate system.
    - **THE VERTICAL AXIS WAS LEVEL INDEX, NOT HEIGHT.** On this exponentially stretched grid
      that is not a height axis at all: layer 0 is 1.2 km and layer 39 is 22.8 km, so index
      space stretched the bottom of the atmosphere and squashed the top by **18.6x** — and the
      distortion varies with altitude, so no aspect-ratio setting in ParaView could undo it.
      Contours, glyph angles and streamline curvature all inherited it. The axis is now
      `get_layer_height(i)` mapped onto the same span: the plot-x to height ratio is constant to
      **1.000000** across the column, and the vertical exaggeration is one number (~30x) instead
      of a function of height.

    **VERIFIED AGAINST THE PHYSICS, NOT JUST AGAINST ITSELF.** The velocity field and
    `dPsi/dz` agree in sign at **100 %** of sampled points with magnitudes matching to **0.2 %**
    (0.63608 from the field against 0.63500 from the streamfunction). This family has had sign
    errors in its diagnostics while its dynamics were correct — the Coriolis case in the
    traps list — so the check was worth making explicitly.

    With all three fixed, one cell traced around its loop reads as a closed circulation for the
    first time: `v` reverses between surface and top on both flanks (-1.09 to +0.95 at 58 N,
    the mirror at 38 N), `u` carries opposite signs on the two flanks at every height
    (-0.11 against +0.10), and `Psi` crosses zero exactly where `v` does, at the cell core.


71. **`ATM_RAD_DIRECT` is the default. The Lambda iteration converges onto it, so the closed
    form is the ANSWER the sweeps were trying to reach — and the shipped 4 sweeps were 15 % high
    at no saving in time.**

    Item 30 wrote the closed-form solve and left it default-off "pending a longer measurement".
    Item 66 ran that measurement on the prognostic column; the flag never moved. This is the
    A/B on the *shipped prescribed* configuration, 40 iterations, 24 threads, arms identical
    apart from `output_path`:

        n_lambda = 4 (was shipped) .... 243.43 W/m2   +15.06 %   277.34 s
        n_lambda = 64 ................. 227.55 W/m2    +7.55 %   301.10 s
        n_lambda = 512 ................ 212.54 W/m2    +0.46 %   467.01 s
        ATM_RAD_DIRECT=1 .............. 211.57 W/m2     0.00 %   280.97 s

    **The convergence is the finding, not the difference.** Two solvers disagreeing says nothing
    about which is right; a Lambda iteration marching 243.43 -> 227.55 -> 212.54 onto 211.57 as
    its sweep count rises says the closed form is what the sweeps compute in the limit, and that
    4 of them do not get there. The exact answer costs **1.3 % more wall clock than the cheapest
    wrong one** and 40 % less than a still-unconverged 512 sweeps.

    This is the twentieth defect's bill arriving. `n_lambda = 4` is an Earth constant of the
    kind this file keeps finding, hidden where nobody looks for one — a loop bound does not read
    as a physical assumption. Items 45-47 and 65 were all produced with it, and item 66 already
    retracted item 45's 2927 W/m2 as its artefact.

    **`ATM_RAD_DIRECT=0` restores the old branch** and is verified against it to the digit
    (243.43 both). Every OLR in this file and in CLAUDE.md predates the flip.

    **ONE DISCREPANCY, RECORDED RATHER THAN EXPLAINED.** Item 30 measured this as a **0.15 %**
    change to the prescribed configuration. It does not reproduce: **13.09 %** on the shipped
    branch and **8.69 %** with `ATM_SKIN_GREY=0`. So item 67 amplified the effect but did not
    create it, and the first hypothesis — that `t_skin` had been pinning the OLR so the solver
    could not move it — is **wrong**, measured. Roughly fifteen defaults have changed since
    item 30 (items 38, 50-53, 57, 59, 60, 63, 67, 68, 70) and no run has been spent attributing
    it. It is a loose end, not a resolved one.

72. **Item 68's 44 % residual is STRUCTURAL, and the proof is that 64x the sweeps changes the
    divergence by nothing at all. The projection has converged — to a fixed point that is not
    divergence-free. NOT SOLVED; five mechanisms measured and excluded.**

    Item 68 left `ATM_PROJ_SWEEPS = 10` as the knee of a curve that plateaus rather than closing,
    with ~44 % of the non-closure beyond the projection's reach. This is the attempt to find out
    why, and it does not find out why. What it establishes is the SHAPE of the thing, which is
    worth more than the five guesses it killed.

    **THE ONE MEASUREMENT THAT MATTERS.** The solver's own residual diagnostic, at 4 iterations,
    24 threads, everything else identical:

        ATM_PRESS_SWEEPS = 1   ->  div(u) rms 3.984e-02   div(rho u)/rho rms 1.611e-02
        ATM_PRESS_SWEEPS = 64  ->  div(u) rms 3.984e-02   div(rho u)/rho rms 1.611e-02

    **Identical to four significant figures under 64x the work**, on a knob that is demonstrably
    live (item 68 verified it: output is not bit-identical and runtime rises). A relaxation that
    does not improve with 64x more sweeps is not converging slowly — it has CONVERGED, and its
    fixed point has non-zero divergence. That is a statement about the operator pair, not about
    the iteration: the discrete divergence and the discrete gradient are not adjoint on this
    stencil, so the Poisson solve finds a pressure whose correction leaves divergence behind no
    matter how exactly it is found.

    **THE MODEL HAS BEEN PRINTING THIS NUMBER ALL ALONG.** `div(rho u)/rho rms = 1.611e-02`,
    identical at the first diagnostic and the last, appears in every run log this project has
    ever produced. The information was never missing; it had simply never been connected to
    `Psi(ground)`. Same lesson as item 68's — the evidence was in the output, not in a run that
    still needed doing.

    **FIVE MECHANISMS, MEASURED AND EXCLUDED.**

    - **The base-state density is not it.** The projection enforces `div(rho_bar u) = 0` with
      `m_dlnrho_dr`, which is **1-D** — a function of radius only — while `Psi` integrates the
      full 3-D `r_humid`. Rebuilding `Psi` from the same slice with `rho_bar(r)` instead of
      `r_humid` moves the ground RMS by **1.0 %**. The density is horizontally uniform to
      0.3–3 % through the bulk of the mass; it varies 11 % only at 218 km, where rho = 0.25
      kg/m³ and carries no flux.
    - **`exp_rm` is not it, and this is the SECOND independent measurement exonerating it.**
      Item 68's zeta scan could not see the residual because the under-converged part dominated.
      Re-run at the new `ATM_PROJ_SWEEPS = 10` default, where the residual IS the signal:
      RMS 1.377e13 / 1.444e13 / **1.840e13** at `zeta` = 3.0 / 1.5 / 0.405, closure ratio
      0.273 / 0.358 / **0.486**. Making the metric EXACT makes the residual **worse**. (The arms
      differ in grid as well as metric — at `zeta = 0.405` the surface layer is 6.1 km against
      1.2 km — so the magnitude is confounded; the direction is not.)
    - **It is not a grid-scale mode in the velocity.** 2-delta oscillation index of `u`, `v` and
      `Psi` is **0.006–0.2** where a checkerboard is ~4.
    - **It is not a Nyquist mode in `Psi(ground)` or in the pressure either.** 2-delta index
      **0.076** for `Psi(ground)` across latitude and **0.002** for `p_dyn` at every level
      sampled, 0.000 in the vertical. Everything is smooth.
    - **It is not a boundary artefact.** Per-layer RMS contribution to `Psi(ground)` is
      distributed through the entire column — the bottom four layers (0–7 km) carry **11.1 %**,
      the 80–160 km band about **33 %**, and the top four layers (218–300 km) **0.5 %**. Every
      layer between contributes 1.5–4.5 %.

    **WHAT IT POINTS AT, AND WHY THAT IS NOT A CONCLUSION.** `PressureSolverAtm.h` already names
    the un-done work, in the note beside `ATM_POISSON_METRIC_FIX`: *"this repairs the metric
    POWER only; the collocated checkerboard (Rhie-Chow face reconstruction) is a separate,
    larger port not done here."* The metric-power half was repaired and bought 18 %
    (rms 2.625e-02 -> 2.153e-02). The Rhie-Chow half never was, and a non-idempotent collocated
    projection is exactly a projection whose fixed point is not divergence-free.

    **But the checkerboard that story predicts is not present** — see the two nulls above — so
    it is a partial match at best and is recorded as a lead, not a diagnosis. **The residual is
    not attributed.**

    **The consequence for the model is not small.** `div(rho u)/rho` at rms 1.6e-02 and max
    1.6e-01 is the state every ATHAD run has integrated in, and it is why `ATM_PROJ_SWEEPS = 10`
    is a knee rather than a cure. Curing it is a numerical-methods port — face-based velocity
    reconstruction — not a constant to correct, and it is the largest single piece of unbuilt
    machinery this file has named.

73. **THE SHIPPED MODEL'S OLR CONVERGES, AND WHAT IT CONVERGES TO IS `F/2` EXACTLY — so the
    imbalance is exactly half the energy input, and it has stopped widening. Item 67's
    follow-up, the grey radiative-equilibrium top, is MEASURED AND DEAD: `T_rad(tau)` exceeds
    the prescribed profile at all 41 levels, so there is no crossing to switch at.**

    Two things were unmeasured after items 67 and 71. Item 67's 200-iteration pair was run with
    `n_lambda = 4` and `ATM_PROJ_SWEEPS = 1`, both since replaced, so no long series described
    the shipped defaults; and item 67 named a one-line follow-up — replace the isothermal clamp
    `max(t_skin, T_ad)` with `max(T_rad(tau_above), T_ad)` — that nobody had run. This is both,
    as a pair: 200 iterations, 24 threads, dry throughout (`moist_phys_start_iter = 300`),
    19.5 min per arm, arms identical apart from `output_path` and the knob.

    ### Arm A: the shipped defaults, for the first time

    `ATM_SKIN_GREY` on, `ATM_RAD_DIRECT` on, `ATM_PROJ_SWEEPS = 10`:

    ```
    iter        20      40      60      80     100     120     140     160     180     200
    OLR     270.62  211.57  177.20  155.13  140.96  135.77  135.85  135.94  136.05  136.19
    sigma*t_skin^4
            135.90  135.58  135.61  135.66  135.70  135.78  135.85  135.95  136.05  136.20
    excess  134.72   75.99   41.59   19.47    5.26   -0.01   -0.00   -0.01    0.00   -0.01
    imbal     0.48   59.59   94.03  116.19  130.48  135.80  135.86  135.96  136.08  136.23
    F/2     135.55  135.58  135.62  135.66  135.72  135.78  135.85  135.95  136.06  136.21
    ```

    **THE DESCENT COMPLETES.** The OLR reaches `sigma*t_skin^4` at iteration 120 and stays
    within 0.01 W/m2 of it for the last five diagnostics. Item 67 measured 11.82 W/m2 of excess
    remaining at 200 and explicitly claimed no limit; that residual was the under-converged
    Lambda solver, not the physics. Item 71's flip is what closed it, and this is the first
    measurement that shows what the flip was worth over a long run: not a 15 % offset, but the
    difference between an asymptote approached and an asymptote reached.

    **AND IT LANDS ON AN IDENTITY.** `sigma*t_skin^4` IS `F/2` by construction under
    `ATM_SKIN_GREY` (item 67), so the converged state is

        OLR = F/2   and therefore   imbalance = F - F/2 = F/2

    The last two rows above are the same number to 0.02 W/m2 at every diagnostic. **The
    corrected model's imbalance is not "+123.76 and widening" — it is +136.23 and CONSTANT, and
    it is constant because it is definitionally half the absorbed flux.** Item 67 said the
    correction makes the closure falsifiable rather than automatic. It does; and the falsified
    version is just as arithmetic as the old one, in the opposite direction. A prescribed
    isothermal lid emits `sigma*t_skin^4` and the column delivers exactly that; whichever number
    `t_skin` is set to, the OLR is that number and the imbalance is what is left over.

    **Verified null, obtained for free.** Arm A's iteration-40 OLR is **211.57 W/m2**,
    digit-for-digit item 71's `ATM_RAD_DIRECT=1` figure measured on the pre-change binary at the
    same 24 threads. The `ATM_SKIN_TAU=0` branch introduced below is a null.

    ### Arm B: `max(T_rad(tau), T_ad)` does not work, in either form

    `ATM_SKIN_TAU` (default 0) implements item 67's follow-up. Under `ATM_SKIN_GREY`,
    `sigma*t_skin^4` is `F/2`, so the profile needs no separate flux argument:

        T_rad(tau) = t_skin * (1 + 3*tau/2)^(1/4)

    Two modes, because the literal form is not the intended one here. Mode **2** is the literal
    `max(T_rad(tau), T_ad)` over the whole column. Mode **1** is the minimal reading — the
    adiabat is kept wherever it already wins (`T_ad >= t_skin`, bit-identical to shipped there)
    and only the isothermal region follows `T_rad(tau)`.

    **Mode 2, 4 iterations, 71 s.** Item 67 wrote the replacement as "keeps the adiabat wherever
    it is warmer". `max tau_above` at the surface is **2.26e6** — printed in every run log this
    project has ever produced — so `T_rad(tau_s) = 221.4*(1 + 1.5*2.26e6)^(1/4)` = **9496 K**
    against a prescribed 1408 K. The literal form does not adjust a lid; it replaces the column.
    And it does not merely sit at 9496 K, it runs away: hotter column -> the hydrostatic
    integration carries more mass aloft -> `tau_above` rises -> hotter `T_rad`. Photosphere
    9668.8 K at iteration 2 and 15 598.8 K at iteration 4; OLR 4.8e8 then 3.0e9 W/m2.

    **Mode 1, 200 iterations.** Stable and sensible for ~150 iterations, then a runaway that
    ignites at the poles. The OLR is never below 1374 W/m2, is non-monotone
    (1809, 1374, 4054, 5554, 4678, 3915, 3205, 2595, 2075) and ends at **22 760**.

    ```
    236.4 km, 90N        it 20     it 100    it 200
      T                   229 K     241 K    2729 K
      tau_above          0.092      0.26        32
    218.2 km, 90N
      T                   327 K     274 K    2793 K
    ```

    Away from the poles the change is exactly what was intended: 236.4 km warms from a uniform
    221.4 K to 239-265 K, smooth, symmetric about the equator to 0.1 K (invariant 1 intact).
    **The ignition is geometric.** The pole holds FIVE levels at `t_skin` where the equator holds
    four, because its adiabat starts 50 K colder; as the column cools through the run the
    crossing descends onto level 36, where `tau` is **1850**, not the 0.15 of level 37. Mode 1
    sets that level to ~1600 K, which puts mass above it, which raises `tau` at levels 37 and 38,
    which sets those hot in turn. It propagates upward. At iteration 200, **18 of 181 latitudes
    at 236.4 km are hotter than their own surface** — which is how the defect was spotted, in
    ParaView, before any log had been read.

    ### Why no version of it works, measured rather than argued

    Arm A's own `tau_above` array, equator column, iteration 200, against
    `T_rad(tau) = t_skin*(1+3*tau/2)^(1/4)`:

    ```
    level    h        T        tau_above    T_rad     T_rad/T
       0      0.0 km  1407.5 K  2.26e+06    9495.8 K    6.75
      10     17.6 km  1330.7 K  1.14e+06    7997.4 K    6.01
      20     54.7 km  1161.3 K  2.28e+05    5352.5 K    4.61
      30    133.4 km   767.9 K  3.99e+03    1946.7 K    2.54
      36    218.2 km   279.7 K  1.63e+03    1555.9 K    5.56
      37    236.4 km   221.4 K  4.97e-01     254.5 K    1.15
      40    300.0 km   221.4 K         0     221.4 K    1.00
    ```

    **`T_rad(tau)` exceeds the profile at every one of the 41 levels** — 6.75x at the surface,
    2.5-5.6x through the middle, 1.15x at the switch, and equal only at the lid. The standard
    radiative-convective construction takes the adiabat below the crossing and the radiative
    profile above it. **There is no crossing.** Mode 1's spike and mode 2's runaway are one fact
    seen at two altitudes: a grey column of `tau_s` = 2.26e6 in radiative equilibrium at
    `F/2` = 136 W/m2 requires a **9496 K** surface, and the model prescribes 1408 K. Read the
    other way, the OLR a 1408 K surface can drive through `tau_s` = 2.26e6 in grey radiative
    equilibrium is `2*sigma*1408^4/(1+1.5*tau_s)` = **0.13 W/m2**.

    So item 67's "cheapest remaining move on invariant 3" is **withdrawn**, and what replaces it
    is a statement about the model rather than about a boundary condition: **the prescribed
    adiabat and the grey opacity are not compatible radiative descriptions of the same
    atmosphere, and they disagree by three orders of magnitude in flux.** One of the three —
    the 1500 K surface, `kappa`, or the greyness — has to give, and the choice cannot be made by
    a one-line change to the top boundary.

    ### One thing that belongs to item 39

    Between level 36 (218.2 km) and level 37 (236.4 km), `tau_above` goes **1626 -> 0.497**.
    That is `dtau` ~ 1600 across a single grid layer, where item 39 recorded 55 and called the
    photosphere "poor, not broken" at three layers. On the current defaults it is worse by 30x,
    the whole opaque-to-transparent transition happens inside one cell, and that is precisely
    why mode 1's switch lands on the wrong side of `tau` = 1. Item 39's conclusion — the lever
    is `zeta`, not `im` — is unaffected; its magnitude is.

    ### The knob

    `ATM_SKIN_TAU` in `ThermoAtm.h`, default **0**, off-branch verified null (see arm A's
    iteration-40 match above). `=1` is the skin-region-only form, `=2` the literal whole-column
    form. Both are kept: they are the two readings of item 67's proposal, and the record of why
    neither is available is worth more than the four lines they cost.

74. **`moist_phys_start_iter = 0` IS THE DEFAULT. The 300 was ATOM_Precipitation's, its
    justification is an Earth diagnosis, and its stated purpose is unattainable here — and
    running without it exposes a create-and-destroy cycle that had been invisible because the
    diagnostic samples on the wrong side of it.**

    ### Why the gate had no standing in this model

    Three findings, none of which needed a run:

    - **It was inherited, not chosen.** `moist_phys_start_iter` entered on the fork commit
      `9984458`; `ATOM_Precipitation/param.py` carries the identical parameter, the identical
      default 300 and a description string identical to the character. Nobody sized it for a
      250 bar water-vapour atmosphere.
    - **Its justification is an Earth surface classification.** The comment reads: *"too stiff
      during the initial velocity transient — a single tropical maritime column repeatedly
      drove q_c, q_i and S_s into runaway"*. Invariant 1 says there is no land, no sea and no
      maritime anything here.
    - **Its stated purpose is not attainable.** The gate waits for "the circulation to form on
      a dry field". Item 18 measured that it does not form: the pressure gradient that should
      balance the meridional wind is at 1.8 % of the Coriolis term after 400 iterations and
      growing linearly, putting geostrophic adjustment of order **1e4** iterations away. The
      gate waits 300 iterations for an event two orders of magnitude further off, then starts
      the stiff physics into a state no more settled than iteration 0 — which is plausibly why
      item 66 saw the onset swing the OLR 223.8 -> 1022 -> 767 -> 751 -> 1036 -> 782.

    And the gate has been open in ATHAD many times already: items 52, 53, 61, 63 and 64 all ran
    at `moist_phys_start_iter = 0` (nm = 5-40), which is how the moist repairs were measured at
    all. **ATHAD_COND has shipped 0 since it was forked.**

    ### The measurement: 200 iterations, 24 threads, no runaway

    23.6 min, exit 0, `q_c`/`q_i`/`S_s` finite throughout — the Earth-inherited runaway is not
    reproduced at five times the length it had ever been tested here.

    **WHAT IT COSTS, AND IT IS NOT SMALL. The moist column does not converge.** Against arm A
    of item 73, identical in every other respect:

    ```
                        iter 140    160     180     200
    dry   (300)  OLR      135.85  135.94  136.05  136.19   arrived at 120, excess 0.01 W/m2
    moist (0)    OLR      189.47  180.97  173.90  168.01   still falling, 32 W/m2 above sigma*t_skin^4
    ```

    The photosphere RISES (239.0 -> 242.8 km) where the dry arm's falls, and the albedo freezes
    at 0.4964 for all ten diagnostics where the dry arm's drifts 0.4986 -> 0.4927. **So every
    dry-column number in this file, item 73's identity included, now belongs to
    `moist_phys_start_iter = 300`**, which is one `<moist_phys_start_iter>300</...>` away.

    ### What the gate was hiding: condensate is created and destroyed every iteration

    The log reports `max cloud water = 0.000000 g/kg` from the first diagnostic onward, which
    reads as "the moist physics deleted the initial deck". **It is a sampling artefact — the
    print sits after the whole moist block.** Sampling inside it (`ATM_ICE_CENSUS=1`, five new
    probe sites) gives the cycle, per iteration:

    ```
    entering the moist block ...................  66 785 cells with cloud water
    post SaturationAdjustment .................. 131 404      (+65 000, it CONDENSES)
    post damp_wiggles .......................... 263 530      (+132 000, a SMOOTHER)
    post ice scheme ............................   6 486      (-257 000)
    post MoistConvection .......................   6 486
    post the cloud_cap clamp ...................   6 486
    ```

    Bisected inside the ice scheme, the loss is entirely in `computeColumns` (263 530 -> 6 486,
    **192 065 cells written to exactly zero**) with the peak value untouched at 37.500000 g/kg.
    Buffer aliasing was excluded first — all 14 array bases distinct — and the S-terms are zero
    at the affected cells, so it is not the microphysics.

    **The routine is `IceSchemeCommon::evaporateWhereImpossible`, called behind `canCondense`,
    and it is CORRECT.** There is no rate and no `dt` in it because there is no phase to relax:

        canCondense = (t_u < T_CRIT_H2O) && (qSatWater < 1.0)

    supercritical or superheated (`p_sat > p`), and in either state a droplet cannot exist at
    all, so the condensate returns to vapour in one step. That is invariant 2 being enforced,
    and the water budget confirms it — total H2O conserved to **-0.0003 %** with precipitation
    identically zero. The zeroed sample cell is at 218.2 km and ~280 K, far below 647 K, so it
    fails the *superheated* test: the air is too thin for water to saturate at any mixing ratio.

    **THE FINDING IS THE DISAGREEMENT UPSTREAM OF IT.** `SaturationAdjustment` condenses into
    ~65 000 cells that `canCondense` then rules impossible, and `damp_wiggles` — a *numerical*
    de-checkerboarding smoother, not a physical process — spreads condensate into another
    ~132 000 of them. Two routines in the same block disagree about where a condensed phase can
    exist, and the ice scheme is the one that is right. **The cycle is not a runaway and not a
    leak; it is work being done and undone every iteration**, and it is why the cloud field
    looks static while the albedo still sees condensate mid-iteration.

    **NOT FIXED, and the next step is to make the two agree** — either `SaturationAdjustment`
    gains the `canCondense` test at its entry, or `damp_wiggles` stops being applied to
    condensate in cells that cannot hold it. Which one is the smaller change is not yet
    measured. What is not in question is the direction: a smoother must not manufacture a phase
    the thermodynamics forbids.

    ### Method note, and it is the third instance

    `max cloud water = 0.000000` was the number on screen, and it was the wrong number: the
    diagnostic sits on the far side of the routine that does the work. Item 42 found this with
    `Psi`, item 68 with `Psi(ground)`'s max-versus-RMS, and here the axis is CALL ORDER rather
    than statistic or coordinate. **Before explaining a field that will not move, find out where
    in the iteration it is being read.**

75. **THE DISAGREEMENT IS INSIDE `SaturationAdjustment`, NOT BETWEEN IT AND THE ICE SCHEME:
    `clampAndFade` tests whether a cell can hold a condensed phase, then condenses, then adds
    the latent heat that makes the answer FALSE — and never re-tests. Fixed
    (`ATM_SAT_SUPERHEAT`, default ON). The OLR moves -49.4 %.**

    Item 74 left the attribution as "two routines disagree" and named `SaturationAdjustment` and
    `damp_wiggles` as the creators. Both halves of that were wrong, and measurement rather than
    source-reading corrected them.

    ### What the census actually says

    Condensate sitting where `IceSchemeCommon::canCondense` forbids it, per iteration:

    ```
    entering the moist block ..............   1 444 cells,    36.3 kg/kg
    post SaturationAdjustment ............. 130 678 cells, 6 520.0 kg/kg     <-- created here
    post damp_wiggles ..................... 263 530 cells, 6 520.0 kg/kg     <-- mass UNCHANGED
    post ice scheme .......................   6 486 cells,    94.8 kg/kg
    ```

    **`damp_wiggles` creates none of it.** The cell count doubles while the mass is identical to
    six figures: a smoother spreads the same condensate over twice as many cells. Item 74's "a
    numerical smoother manufactures a phase the thermodynamics forbids" is **withdrawn on mass**;
    what it does is redistribute, which matters for the cell count and not for the budget.

    **And the two predicates AGREE.** `canCondense` uses `saturationPressure`,
    `SaturationAdjustment` uses `saturationPressureAuto`; the obvious story is that they pick
    different curves either side of the triple point. Measured at the offending cells, they do
    not disagree at all:

    ```
    forbid sample [39][23][0]  T=276.08 K  p=1.225 hPa   E_liq=7.542     qsat_liq=1.000000  qsat_auto=1.000000
    forbid sample [37][0][0]   T=461.15 K  p=51.18 hPa   E_liq=1.201e+04 qsat_liq=1.000000  qsat_auto=1.000000
    ```

    Those cells are genuinely superheated — `p_sat` is 6x and 235x the local pressure — and every
    violation in the census is superheated, never supercritical, at every stage. The `T_CRIT`
    guards work; the `p_sat > p` half is the one that was missing.

    ### Where it is created, and why the guard that exists does not stop it

    `SaturationAdjustment::clampAndFade` **already carries the exact test** and calls
    `IceSchemeCommon::evaporateWhereImpossible` on cells that fail it. It is applied to the
    temperature the cell *arrives* with. Forty lines below, the "always-on supersaturation
    removal" does this:

    ```cpp
    if (c_row[k] > q_sat) {
        const double excess = c_row[k] - q_sat;
        c_row[k] = q_sat;
        cloud_row[k] += excess;
        T_dim        += latentHeat(T_dim) * excess / cp_of(...);   // <-- raises T
        t_row_nd[k] = T_dim * inv_t_0;
    }
    ```

    Raising T raises `E_sat`, and once `E_sat` exceeds the local pressure the cell is superheated:
    the cloud just written into it is a state the guard above would have refused. Nothing
    re-tests it, so the step stands. **The routine's own latent heating invalidates its own
    admissibility test.**

    This is the argument the file already makes at the critical-point cap — *"a condensation
    adjustment cannot legitimately heat one past the point where the phase it is condensing into
    ceases to exist"* — applied to only one of the two ways a cell stops being able to hold a
    condensed phase.

    **Staged measurement, which is what localised it** (`ATM_ICE_CENSUS=1`, counters inside
    `SaturationAdjustment::run`):

    ```
    call 2, before the fix:  entry 1 444 | post adjustSat 722 | post applyTopo 722 | post clampAndFade 130 678
    call 2, after  the fix:  entry 117 325 | post adjustSat 48 735 | post applyTopo 48 735 | post clampAndFade 0
    ```

    A first attempt put the guard in `adjustSaturation`'s Newton write-back. It fired **128 881
    times** and changed the census by **nothing** — the creator was the other loop. *A guard that
    fires is not a guard that matters; count the thing you are trying to remove, not the
    invocations.*

    ### The fix

    `ATM_SAT_SUPERHEAT` (**default ON**, `=0` restores the old behaviour) re-tests after the
    write-back and **rejects the whole step** — the cell returns to its entry state. Both loops
    carry it. A rejected cell keeps its supersaturation, which is the honest state of a parcel
    that cannot condense (item 64), rather than a manufactured phase the next routine deletes.
    Condensing only as far as `E_sat(T) = p` would need the joint solve `adjustSaturation`'s
    Newton loop does; that is a bigger change and is NOT attempted here.

    ### What it costs, and it is the largest effect in this file

    40 iterations, 24 threads, moist physics on from iteration 0 in both arms, identical
    otherwise. Item 61's reproducibility envelope is +-0.5 % on a 40-iteration OLR; this is 100x
    that, so it is a property of the model and not of the thread count.

    | | `ATM_SAT_SUPERHEAT=0` | ON (default) |
    |---|---|---|
    | OLR @20 / @40 | 297.71 / **267.84** | 135.85 / **135.52** |
    | imbalance @40 | 3.24 | 135.50 |
    | photosphere | 237.3 km / 363.4 K | **281.3 km / 225.7 K** |
    | emission from isothermal skin | 1.1 % | **79.0 %** |
    | max cloud water @40 | **0.000000 g/kg** | 40.09 g/kg |
    | max cloud ice @40 | 0.000000 g/kg | 8.15 g/kg |
    | max water vapour @40 | 690.34 g/kg | **778.11 g/kg** |
    | cloud cells surviving the ice scheme | 6 486 | **149 271** |

    **-49.4 % on the OLR**, against item 67's -27 % and the 0.02-0.11 % of every microphysics
    repair before it. The mechanism is opacity: condensate now survives the iteration, so
    `k_liq*LWP + k_ice*IWP` is real, the emission level rises 44 km and cools 138 K.

    **THREE HONEST CONSEQUENCES, none of them comfortable.**

    - **The model is MORE `t_skin`-dominated, not less.** 135.52 W/m2 is `sigma*t_skin^4` to two
      decimals, so item 73's `OLR = F/2` identity now arrives by iteration **40** instead of 120,
      and **79 % of columns radiate from the prescribed isothermal lid** against 1.1 %. Raising
      the photosphere pushed it into the skin — item 43's mechanism, driven this time by cloud
      opacity rather than by a warm lid.
    - **The retained supersaturation is real and it grows.** Max `q_v` reaches **778 g/kg**
      against the old branch's 690 and the design mass fraction of 672. Rejecting a step leaves
      the excess in the vapour, which is the correct state for a parcel that cannot condense, but
      it is item 64's supersaturation made larger and it has not been run past 40 iterations.
    - **`damp_wiggles` still spreads condensate into forbidden cells** — 167 504 cells and
      3 970 kg/kg after the fix — which the ice scheme still evaporates. A smaller residue of the
      same cycle, in a generic smoothing utility used for many fields. **Not addressed.**

76. **ZERO PRECIPITATION IS THE RIGHT ANSWER, AND FOR A BETTER REASON THAN "NOTHING FORMS":
    rain forms at its 260 mm/d CAP and is completely evaporated by the supercritical column
    before it can reach a 1490 K surface. The Earth phase bands were repaired anyway, and the
    repair is a measured no-op — the third time a suspected blocker has not been the blocker.**

    Every precipitation field in every ATHAD run reads `0.000000 mm/d`. Two questions hide behind
    that: should it be zero at the GROUND (yes, necessarily), and should it be zero EVERYWHERE
    (no — the conversion terms are live: `max S_r = 0.0160 g/kg/s` at 256 km and
    `max S_s = 0.0119` at 277 km).

    ### The phase bands are Earth's, and they did two jobs at once

    ```cpp
    m.P_rain    = (t_u >= t_0)                  ? (inherited + produced) : 0.0;   // >= 273.15 K
    m.P_snow    = (t_u < t_0 && t_u >= t_000)   ? (inherited + produced) : 0.0;   // >= 253.15 K
    m.P_graupel = (t_u < t_0 && t_u >= t_00)    ? (inherited + produced) : 0.0;   // >= 236.15 K
    ```

    Two defects in one expression. The bands **partition down to 236.15 K and no further**, and
    `t_00` = -37 C and `t_000` = -20 C are Earth's homogeneous-freezing and mixed-phase
    thresholds — the same defect class as item 63's convective triggers in absolute hPa. And the
    `: 0.0` conflates *production* with *transmission*: a flux arriving from the level above is
    DESTROYED rather than merely not added to, though falling ice does not cease to exist because
    the air it is passing through is cold.

    **Repaired** (`ATM_PRECIP_BANDS`, default ON, `=0` restores): the inherited flux always
    passes and only production is gated; snow production loses its `t_000` floor, since ice
    crystals form and fall at any temperature below `t_0`. Graupel keeps its `t_00` floor because
    that one is physical rather than Earth-specific — riming needs supercooled LIQUID, and below
    -37 C there is none.

    ### And it changes nothing, measured

    40 iterations, 24 threads: OLR **135.52 both**, max water vapour **778.109092 g/kg both**,
    max cloud water **40.087238 g/kg both** — bit-identical. Instrumented, the reason is flat:

    ```
    PRECIP PROBE: cells with S_s[i+1] > 0 .... 112 726
                  cells with prod_s > 0 ......       0
                  max P_snow ................. 0.000000e+00
                  max P_rain ................. 3.000000e-03   <-- P_max_flux, the 260 mm/d cap
    ```

    **`prod_s` is zero in every cell**, so every cell in which the flux loop evaluates production
    is at or above 273.15 K and the snow and graupel bands never applied there to begin with. The
    column's cold part — the 221-230 K condensing layer — is not where the flux recurrence is
    accumulating. **The bands were not the blocker.** Item 63 found the same shape (convective
    triggers repaired, nothing moved, because invariant 2 pinned the layer) and item 74 found it
    again (the `dt` that was not missing). *A defect can be real, worth repairing, and not be the
    thing you were looking for.*

    ### Why the answer is zero, which is the part that matters

    `max P_rain = 3.0e-3 kg/(m2 s)` INSIDE `ThreeCatIceScheme` — the flux is not merely nonzero,
    it is **pinned at its cap** — and `0.000000` by the time any diagnostic reads it. What
    happens in between is `IceSchemeCommon::evaporateWhereImpossible`, which converts an incoming
    flux to vapour at `P/(v*rho)` when it falls into a superheated or supercritical cell, mass
    accounted for. Everything below ~177 km is supercritical, so nothing survives the descent.

    That is confirmed by the budget rather than assumed: total H2O conserved to **-0.0004 %**
    with **0.000000 kg/kg** deleted at the `c` ceiling, in a 200-iteration run. The water that
    falls is not lost; it is returned to the vapour a few kilometres down.

    **So ATHAD's dry surface is a RESULT, not an absence.** ~~Rain forms in this atmosphere, at
    the maximum rate the scheme allows, and evaporates before it can arrive.~~ **CORRECTED BY
    ITEM 77: that describes the FIRST ice-scheme call only.** The cap binds on call 1, on
    `initCloudIce`'s deck, and never again — from call 2 onward no rain flux is produced at all,
    so there is nothing falling for the supercritical column to evaporate. The conclusion stands
    and its mechanism does not: the flux is never launched rather than launched and evaporated.
    A maximum quoted without asking WHEN it occurred. That is the physical
    distinction between this epoch and ATHAD_COND's, and the model computes it rather than
    asserting it.

    **Open, and NOT investigated here**: `P_rain` sitting exactly on `P_max_flux` is a cap
    binding, and this file's own rule is to look at what a cap is holding back (item 52). The
    normalisation beneath it is also suspect — `P_rain_0 = max(P_rain.x[0][j][k], 1e-6)` divides
    the flux by its SURFACE value, which in ATHAD is structurally zero, so the collection kernels
    see `P/1e-6`. `ATM_PRECIP_DIMENSIONAL` exists for exactly this and is default off. Neither is
    measured.

77. **THE P_rain CAP IS HOLDING BACK NOTHING — it binds on the FIRST ice-scheme call and never
    again. Which corrects item 76: rain is not "forming at its cap and evaporating", it is not
    being produced at all after the initial deck is consumed, because the only cells with a rain
    SOURCE are 50 K too cold for liquid.**

    Item 52's rule is that a binding cap must be interrogated, because `MC_t` sat pinned at
    `MCt_max` while updraft parcels ran at 46 000 K and every plot looked physical. Item 76 saw
    `max P_rain` sitting exactly on `P_max_flux` and named it as the next thing to check.
    `ATM_PRECIP_CAP` scales the cap (default 1.0 = the shipped 3.0e-3 kg/(m2 s) ~ 260 mm/d), and
    the census counts cells at it and the largest value the recurrence wanted.

    ### The measurement, 8 iterations, 24 threads

    ```
    ice-scheme call 1:  52 055 cells at the cap,  largest wanted 8.484496e-03   (2.83x the cap)
    ice-scheme call 2:       0 cells,             largest wanted 0.000000e+00
    ice-scheme call 3:       0 cells,             largest wanted 0.000000e+00
    ice-scheme call 4:       0 cells,             largest wanted 0.000000e+00
    ice-scheme call 5:       0 cells,             largest wanted 0.000000e+00
    ```

    **It binds once.** The first call autoconverts `initCloudIce`'s 37.5 g/kg deck and the
    recurrence asks for 2.83x the cap; from the second call onward the rain flux is not capped,
    it is **absent**. So this cap is not item 52's shape: nothing is hiding underneath it, and
    the field beneath it is trustworthy — which is a real answer, not a null result, because the
    alternative was a scheme running pinned against a backstop forever.

    ### The correction to item 76

    Item 76 reported `max P_rain = 3.0e-3` measured inside `ThreeCatIceScheme` and concluded that
    **"rain forms at the maximum rate the scheme allows and evaporates before it can reach a
    1490 K ground"**. That number came from the FIRST call. It is right about the first call and
    **wrong as a statement about the model's state**: after iteration 1 no rain flux is produced
    anywhere, so there is nothing for the supercritical column to evaporate.

    What is NOT retracted is item 76's conclusion that the dry surface is a computed result: the
    condensate is real, it is removed where it cannot exist, and the water budget closes to
    -0.0004 % over 200 iterations. What changes is the mechanism — **the flux does not fall and
    evaporate; it is never launched.**

    ### Why, and the new question it opens

    `prod_r` is gated on `t_u >= t_0` = 273.15 K, and the gate is correct: there is no liquid
    rain at 221 K. The rain SOURCE, however, is nonzero exactly there —

    ```
    max S_r = 0.019869 g/kg/s  at 0 N, 256 027 m    (T ~ 221-230 K)
    max S_s = 0.016272 g/kg/s  at 65 S, 277 192 m
    max S_c = 0.000032 g/kg/s
    ```

    **~~A rain source 50 K below freezing is the thing to explain.~~ WITHDRAWN BY ITEM 78:**
    decomposed at the maximum, `S_r` is entirely `S_c_au` at **280.51 K** — the cell is above
    freezing when the ice scheme reads it, and 221 K only after `densities()` overwrites `t` at
    the end of the iteration. There is no anomaly; there was a comparison between two different
    temperatures. The original claim follows.
    **A rain source 50 K below freezing is the thing to explain.** `S_c_au`, the autoconversion,
    carries its own `t_u >= m.t_0` guard and cannot be the contributor; `S_r` is a sum, so the
    nonzero part must come from accretion or from snow/graupel shedding — every one of which
    needs liquid water that cannot exist at 221 K. **Not diagnosed here.** The candidates are
    exactly the terms item 76 flagged as unexamined, and the normalisation beneath them
    (`P_rain_0 = max(P_rain.x[0][j][k], 1e-6)`, dividing by a surface flux that is structurally
    zero in this model) is still default-off behind `ATM_PRECIP_DIMENSIONAL`.

    **Method note.** The cap census was three lines and answered in 80 seconds a question that
    two runs and a README item had been circling. It also caught the item-76 error, which was
    mine: I quoted a maximum without asking WHEN it occurred. Item 68's lesson was that a max
    hides a field; this is the same lesson in time rather than space — **a maximum over a run is
    not a property of the run.**

78. **THERE IS NO RAIN SOURCE AT 221 K. There is a rain source at 280.5 K, in the same cell, in
    the same iteration — because the moist physics heats the top of the column by ~60 K and
    `densities()` discards that heating before anything else can read it.**

    Item 77 closed the cap question and opened this one: `max S_r = 0.0199 g/kg/s` at 256 km,
    where the diagnostics report ~221 K, with `S_c_au` provably gated off below `t_0`. Decomposed
    at the maximum, under `ATM_ICE_CENSUS`:

    ```
    max S_r = 1.649e-05 kg/(kg s) at [38][38][358]   t_u = 280.51 K
      S_c_au   1.649e-05      S_ac     0.000e+00     -S_ev     0.000e+00
      S_s_shed 0.000e+00      S_g_shed 0.000e+00     -S_r_cri  0.000e+00
      -S_r_frz 0.000e+00      S_s_melt 0.000e+00     S_g_melt  0.000e+00
    ```

    **The whole of `S_r` is `S_c_au`, at 280.51 K** — ordinary warm-rain autoconversion, in a cell
    above freezing, exactly where that term is licensed to run. Every other contributor is
    identically zero, as the gates say they must be. **Item 77's open question is withdrawn: it
    was not a finding, it was a comparison between two different temperatures.**

    The 221 K comes from `Results_Atm`, which reads `t` after `ThermoAtm::densities()` has
    re-imposed the prescribed adiabat and its isothermal skin. The ice scheme reads `t` in the
    middle of the same iteration. **Third instance this session of the same error** — after
    `max cloud water = 0.000000` (item 74, the print sits after the whole moist block) and item
    76's cap value (item 77, a maximum quoted without asking *when*). *Where a diagnostic samples
    decides what it means, and the axis can be call order, space, or time.*

    ### What is real, and it is a statement about invariant 3

    The same cell is **221 K in every diagnostic and 280.5 K inside the ice scheme**, and the call
    order in `run_3D_loop` says why:

    ```
    MultiLayerRadiation  (line 1446) ...... reads t as densities() left it: 221 K
    RK4 dynamics
    moist block          (line 1685) ...... SaturationAdjustment condenses, latent heat -> 280 K
                                            ThreeCatIceScheme autoconverts at 280 K
    ConvectiveAdjustment (line 2124)
    ThermoAtm::densities (line 2125) ...... overwrites t with max(t_skin, T_ad): back to 221 K
    diagnostics          (line 2129) ...... read 221 K
    ```

    **The latent heating exists in a window of one iteration and reaches nothing.** It is computed
    every iteration, it drives the microphysics, and it is discarded before the radiation, the
    diagnostics or the next iteration can see it.

    CLAUDE.md's invariant 3 says the prescription overwrites *what the radiation computed*. This
    is stronger and had not been stated: **it also overwrites what the LATENT HEAT RELEASE
    computed, so the moist physics cannot warm this column at all, by construction.** A 250 bar
    water-vapour atmosphere condensing 40 g/kg at its top is not permitted to raise its own
    temperature by one kelvin between iterations.

    And it closes the loop on items 74-75: `densities()` resets the cell to 221 K, which is far
    below saturation for the vapour there, so `SaturationAdjustment` condenses hard on the next
    pass, releases 60 K of latent heat, the ice scheme acts on the warmed cell, and `densities()`
    resets it again. **The same work done and undone every iteration, one level up from item 74's
    condensate cycle** — and this one is not a defect in a routine, it is the prescribed profile
    doing exactly what it is written to do.

    **`ATM_PROGNOSTIC_T=1` is the switch that stops it**, and item 66 measured that branch's dry
    column converging at 223.8 W/m2. What that branch does with the moist physics on, now that
    items 75-77 have made the condensate real, is not measured. **That is the run this line of
    work has been walking toward**, and it is the one to do next rather than another probe.

79. **THE FREE-RUNNING MOIST COLUMN: the OLR is decoupled from `t_skin` at last, the photosphere
    is a real emission level that barely moves, and the flux is HALF the dry column's — 103
    against item 66's 223.8 W/m2. It is NOT converged and no limit is claimed.**

    Every moist OLR this file carries is `sigma*t_skin^4` read back out (items 43, 67, 73), and
    item 75 made it worse rather than better: with the condensate real, **79 % of columns** radiate
    from the prescribed lid against 1.1 % before. `ATM_PROGNOSTIC_T=1` removes the lid, and item
    66 measured that branch's DRY column converging at 223.8 W/m2. What it does with the moist
    physics on — now that items 75-77 have made the condensate survive an iteration — had never
    been run. 400 iterations, `ATM_PROGNOSTIC_T=1 ATM_RAD_DIRECT=1`, `moist_phys_start_iter = 0`,
    24 threads, 44.5 min.

    ### The measurement

    ```
    iter   20    40    60    80   100   120   140   160   180   200
    OLR  92.61 91.85 91.39 91.12 90.58 90.63 91.32 91.96 92.70 92.95
    iter  220   240   260   280   300   320   340   360   380   400
    OLR  94.24 95.54 96.61 97.52 98.62 99.19 100.33 101.50 102.52 102.99

    photosphere  282.6 -> 283.1 km      T_ph  260.58 -> 260.13 K
    emission from isothermal skin       0.0 % at ALL TWENTY diagnostics
    absorbed SW + geothermal            270.96 W/m2      imbalance +167.98
    ```

    **THE PIN IS GONE, AND THIS TIME WITHOUT THE COLUMN FLYING OFF.** `sigma*t_skin^4` is
    135.5 W/m2 and this column sits at 103, never approaching it, with `skin%` at 0.0 throughout —
    the mechanism items 43 and 73 measured cannot operate where there is no isothermal lid. Item
    45 reached the same escape and found 2927 W/m2 still falling; item 66 showed that was the
    `n_lambda = 4` solver. With the converged solver AND the moist physics, the escape lands
    somewhere physical instead.

    **THE PHOTOSPHERE IS THE STRONGEST PART.** 282.6 -> 283.1 km and 260.6 -> 260.1 K across 400
    iterations: **half a kilometre and two kelvin**. That is a genuine radiating level, held by the
    column's own opacity, and it is what a converged photosphere should look like. Contrast the
    prescribed arm, where the "emission level" is 79 % lid.

    **AND IT IS NOT CONVERGED.** The OLR falls to a minimum of 90.58 at iteration 100, then rises
    monotonically for 300 iterations, still gaining **+0.47 W/m2 over the last 20**. A rising,
    slightly decelerating trend is exactly what items 30 and 45/66 were each caught extrapolating.
    **No limit is claimed. 102.99 is where it is at 400, not where it is going.**

    ### What it says, stated carefully

    - **The moist prognostic column radiates less than half the dry one** — 103 against 223.8 at
      the same iteration count and solver. That is the condensate items 75-77 made real: the
      photosphere sits at 283 km where the air is 260 K, instead of in the dry column's warmer,
      lower emission region.
    - **The imbalance is +168 W/m2**: absorbing 271 and emitting 103, so the column radiates
      **38 %** of what it takes in. Worse-looking than the prescribed arm's +135.5 and more
      honest than it, which is invariant 3's warning arriving for the fourth time.
    - **Stable, not merely finite.** Surface floats 1498.7 -> 1491.6 K, max cloud water 30.9 g/kg,
      cloud ice 14.0 g/kg, water conserved to **+0.0007 %** with 0.000000 kg/kg deleted at the
      `c` ceiling. Item 64's free-running top at supersaturation ratio 1934 does not recur.
    - **Open, and rising: `max q_v` reaches 794.7 g/kg** and is still climbing, where the
      prescribed arm plateaued at 778.15 (item 75). Retained supersaturation is the cost of item
      75's whole-step rejection, and on this branch it has not plateaued within 400 iterations.

    ### Two defects on the output side, found by running with plots on

    - **`restart_stride = 0` AND `checkpoint_save_iter = -1` do not disable the restart dump.**
      Four `atm_restart_0Ma_*.bin` at 593 MB each were written with both knobs off. Item 45
      recorded the first half of this as an aside — *"`restart_stride = 0` did NOT disable the
      restart dump"* — and it is a live defect costing **2.3 GB per run**, not a footnote.
    - **The output file indices do not match the iteration count.** A 400-iteration run wrote
      `panorama_900.vts` and `atm_restart_0Ma_800.bin`, so the writers key off `total_iter_count`
      rather than the loop index. Match a file to a diagnostic with care.
    - Sizes, for planning: 15 GB total, of which the panorama is **8.3 GB in 9 files**
      (`paraview_panorama_vts_flag`, independent of `checkpoint`), the zonal/radial/longal slices
      4.1 GB in 138 files, and the restarts 2.3 GB. Turning the panorama off cuts a run to ~6 GB
      without touching the slices anything here is actually read from.

80. **THE RADIAL METRIC IS REPAIRED AND THE REPAIR IS DEFAULT-OFF, BECAUSE THE CORRECT JACOBIAN
    MAKES THE PROJECTION WORSE — measured unconfounded this time. `checkRadialMetric()` goes from
    a 11.8x spread to 1.00x, the OLR does not move at all, and `div(rho u)/rho` triples.**

    Item 39 identified `exp_rm = 1/(rm+1)` as the Jacobian of a QUADRATIC stretch applied to an
    EXPONENTIAL grid, documented as the real thing in two files that agree with each other, and
    left the repair as an open decision: replace it (~91 sites, moves every number) or cut `zeta`
    to `ln(1.5)` so the wrong formula becomes accidentally right. This is the replacement.

    ### What the metric actually is

    `init_layer_heights()` builds `z(r) = (exp(zeta*(r - r0)) - 1) * L_atm`, so

        J(r) = dz/d(rad.z) = zeta * L_atm * exp(zeta*(r - r0))        [m per rad.z unit]

    and the core wants it dimensionless against its own length unit, `metricShellLength()`:

        exp_rm_true = metricShellLength() / J(r)

    which runs **6.36 at the surface to 0.317 at the top** against the shipped 0.5 -> 0.333.
    That is **12.7x at the bottom and 0.95x at the top** — the 11.8x spread item 39 measured,
    seen as an absolute error rather than a ratio.

    ### AND THE SECOND DERIVATIVE WAS MISSING A TERM IN BOTH METRICS

    Found while repairing the first. With `e = U/J`,

        U^2 * d2f/dz2 = e^2 * ( d2f/dr2 - (J'/J) * df/dr )

    and the core computes `d2f/dr2 * exp_2_rm` and stops. `J'/J` is **zeta** for the exponential
    stretch and `1/(rm+1)` for the quadratic one, so the term is the same order as what is kept
    under the true metric (zeta = 3) and merely small under the legacy one. **This is a second,
    independent defect in the same expression**, and it is why the repair is not a substitution.
    `metricCurv()` returns 0 on the legacy branch, so the legacy operator is left exactly as it
    was and its own missing term is recorded rather than silently changed.

    ### What was changed

    - `metricExpRm(rm)`, `metricCurv(rm)`, `metricExact()` on `cAtmosphereModel`, one definition
      for the whole model, replacing **11 sites** that each recomputed `1.0/(rm+1.0)` locally
      (RungeKutta, PressureSolver x3, ThermoAtm x2, TurbulenceAtm x3, UtilsAtm, checkRadialMetric).
      The other 84 `exp_rm` references are uses and needed no change.
    - The curvature term at all **11 live second-derivative sites** in `RHS_Atm_Turb.cpp`:
      `d2Xdr2 * exp_2_rm` -> `(d2Xdr2 - curv * dXdr) * exp_2_rm`.
    - The **Poisson operator**: the radial Laplacian's first-derivative coefficient folds into the
      existing anelastic off-diagonal, `num_a = exp_2_rm * (dlnrho - curv) * inv_2dr`. Same
      7-point stencil, no new solver; the added ratio to `num1` is `curv*dr/2` = 0.0375.
    - `TurbulenceAtm`'s two `exp_2_rm` are **dead** — declared, never used, and among the
      unused-variable warnings the build has always printed. No curvature term needed there.

    ### The off-branch null, stated precisely

    40 iterations, 24 threads: OLR, max temperature, max water vapour and max cloud water
    **identical to every printed digit**. `residuum_atm` differs in the last digit
    (1.34243666 against 1.34243665, 1e-8). That is item 61's documented "rebuild that re-arranges
    arithmetic in a physics loop" category — `(a - 0*b)*c` and `e*(x - 0)*y` reassociate — so this
    is a null to the printed precision of the physics, **not a bit-identity claim.**

    ### The repair is confirmed by the model's own instrument

    ```
    AGCM: radial metric check - core length unit runs 300287 m at the surface to
          300287 m at the top, spread 1.00x   (1.00 = exp_rm is the Jacobian)
    ```

    Against 11.8x before. `checkRadialMetric()` is unit-free and was written in item 39, before
    this fix existed, so it cannot have been tuned to it.

    ### And the physics says: don't ship it yet

    40 iterations, 24 threads, everything else identical:

    | | legacy | exact metric |
    |---|---|---|
    | OLR | 135.52 | **135.52** |
    | max water vapour | 778.109 | 775.538 g/kg |
    | max cloud water | 40.087 | 39.998 g/kg |
    | max v | 3.463 | 3.344 m/s |
    | max w | 25.056 | 25.051 m/s |
    | `div(rho u)/rho` rms | 2.739e-02 | **7.722e-02** |
    | `Psi_max` | 52 153 @ lat 19, **z = 36 470 m** | 158 088 @ lat 5, **z = 0 m** |

    **The OLR does not move at all** — a 12.7x change in the surface radial scaling, and the
    radiative answer is identical to six figures. That is items 29/41/60 confirmed from a fifth
    direction: the metric is not what the OLR was waiting for either.

    **But the projection gets 2.8x worse, and `Psi_max` jumps to the GROUND** — item 68's
    signature of a streamfunction that does not close. **This is item 72's zeta-scan observation
    reproduced UNCONFOUNDED**: that scan changed the grid as well as the metric and could only
    claim a direction; this changes one formula with the grid fixed, and the direction holds.

    **Part of it is solver stiffness and part is not.** `exp_rm` = 5.90 at the surface against
    0.169, so the radial coefficient `exp_2_rm*inv_dr2` grows ~1200x and the elliptic problem
    becomes far more anisotropic. Measured:

        exact metric, default sweeps ....... 7.722e-02
        exact metric, ATM_PRESS_SWEEPS=64 .. 5.258e-02      (-32 %)

    So unlike the legacy branch — where item 72 showed 64x the sweeps changes the divergence by
    **nothing** — this operator IS under-converged and responds to more work. But 5.258e-02 is
    still **~2x worse than the legacy 2.739e-02**, so stiffness is not the whole account.

    **DEFAULT OFF, and for a measured reason rather than caution.** The formula is right and the
    instrument agrees it is right; shipping it would triple the divergence the model integrates
    in. The open question is what the residual 2x is: a genuinely harder elliptic problem that
    needs a better solver, or a second inconsistency where something else in the discretisation
    still assumes the quadratic form. **Item 72's Rhie-Chow lead and this residual are now the
    same question**, and they should be worked together rather than separately.

81. **`ATM_CELL_ALTERNATE` IS MEASURED IN ATHAD AND ATHAD_COND NOW, AND IT IS A STRUCTURAL KNOB
    WITH NO SCALAR SIGNATURE: the sign-band count of Psi changes across the whole depth of the
    circulation while the OLR, the photosphere, the albedo, the water vapour and `Psi_max` do not
    move at all. In ATHAD_COND the off arm is a SINGLE overturning cell and the on arm is four.**

    Items 69 and `073e4ac`/`ff8fd94` left the knob default-ON in all three trees with an explicit
    "UNMEASURED here" banner on two of them: the effect had been measured only in ATHAD_PERID.
    This is the A/B in the other two. 4 runs x 40 iterations, 24 threads, one freshly relinked
    binary per tree, arms separated by the environment variable alone
    (`config_alt_off.xml` / `config_alt_on.xml`, panorama and checkpoints off, `nm = 40`).
    28 min wall: ATHAD ~5 min/arm at `im = 41`, ATHAD_COND ~7-10 min/arm at `im = 61`.

    The knob acted as designed in both trees, from the startup banner:

        off: direct(-3.0)/indir.(4.0)/indir.(4.0)/indir.(4.0)/indir.(0.5)   1 of 4 pairs counter-rotate
        on:  direct(-3.0)/indir.(4.0)/direct(-4.0)/indir.(4.0)/direct(-0.5) 4 of 4 pairs counter-rotate

    ### What moved: the band structure, and it moved through the depth

    Sign bands of Psi across the northern hemisphere at iteration 40 (`python/survival.py`),
    identical at iteration 20:

    | ATHAD, `im` = 41 | 0 km | 5 | 13 | 23 | 36 | 55 | 79 | 113 | 158 |
    |---|---|---|---|---|---|---|---|---|---|
    | off | 6 | 6 | 6 | 6 | 6 | 6 | **2** | **2** | **2** |
    | on  | 5 | 5 | 5 | 5 | 5 | 5 | **5** | **4** | **4** |

    | ATHAD_COND, `im` = 61 | 0 km | 2 | 5 | 9 | 15 | 22 | 32 | 45 |
    |---|---|---|---|---|---|---|---|---|
    | off | **1** | **1** | 2 | 2 | 2 | 2 | 2 | 2 |
    | on  | **4** | **4** | 5 | 5 | 5 | 5 | 4 | 4 |

    **The band count is a thresholded quantity, so it was checked against its threshold before
    being believed.** `bands()` ignores latitudes carrying less than `thr * max|Psi|`; over
    `thr` = 1e-3, 1e-2, 5e-2 and 1e-1 the separation is stable in both trees — ATHAD off is 2
    bands in the upper column at every threshold and ATHAD_COND off is 1 band at the ground at
    every threshold, while both on arms hold 4-6 through the circulation. What the threshold does
    move is the top of the column, where |Psi| is small and the count drops to 0; that is the
    instrument running out of signal, not the cells ending.

    **ATHAD_COND is where the knob pays.** Off, its zonal-mean circulation IS a single overturning
    at the ground — the co-rotating stack's mass fluxes adding, exactly as item 69 predicted they
    would — and the prescribed five cells are simply not present in Psi. On, four survive to
    iteration 40. ATHAD's own answer is weaker and mixed: it keeps a multi-cell structure below
    ~55 km either way, and the knob buys the upper half of the column (2 bands -> 4-5 above
    79 km). This is `ATM_MC_GEOPOTENTIAL`'s pattern again (item 63) — the machinery is right and
    ATHAD is the tree with less room for it to act.

    **Unexplained, and stated rather than smoothed over**: ATHAD's *off* arm shows **6** low-level
    bands against the on arm's 5 — more sign changes from the co-rotating stack, not fewer, stable
    across all four thresholds. ATHAD_COND does not do this. No mechanism offered.

    ### What did not move: everything scalar

    Iteration 40, 24 threads:

    | | ATHAD off | ATHAD on | COND off | COND on |
    |---|---|---|---|---|
    | OLR | 135.52 | **135.52** | 196.42 | **196.42** |
    | imbalance | 135.50 | 135.50 | 74.36 | 74.36 |
    | photosphere | 281.3 km / 225.74 K | identical | 65.8 km / 264.71 K | identical |
    | mean albedo | 0.4990 | 0.4990 | 0.5000 | 0.5000 |
    | max water vapour | 778.109092 | **778.109092** g/kg | 352.430272 | **352.430272** g/kg |
    | max `w` | 25.056255 | 25.056123 m/s | 27.736137 | **27.736137** m/s |
    | max `v` | 3.463144 | 3.507947 m/s | 3.591077 | 3.591004 m/s |
    | global max Psi | 5.21537e13 | 5.39721e13 (+3.5 %) | 2.74297e13 | **2.74295e13** |
    | max Psi(ground) | 3.33491e13 | 3.40814e13 | 2.03639e13 | **2.03637e13** |
    | RMS Psi(ground) | 1.32398e13 | 1.34996e13 (+2.0 %) | 1.03115e13 | **9.62170e12 (-6.7 %)** |
    | `div(rho u)/rho` rms | 2.739e-02 | **2.884e-02 (+5.3 %)** | 2.825e-02 | **2.928e-02 (+3.6 %)** |

    **The ATHAD_COND column is the whole methodological point of this item.** A change that turns
    one overturning cell into four leaves `Psi_max` agreeing to **five significant figures** and
    max `Psi(ground)` to five as well; only the RMS over latitude moves, by 6.7 %. The reason is
    item 69's geometry: the Hadley cell is `k = 0`, which is even and therefore already direct, so
    the knob cannot touch the cell the global maximum lives in. **Judging this knob by `Psi_max`
    returns a false null**, and the same trap caught ATHAD_PERID first. Read `survival.py` for the
    band count and `psicheck.py` — or the RMS over latitude of `meridional_streamfunction_*.csv`
    — for the closure metric. This is item 68's method note ("`Psi(ground)` must be read as an RMS
    over latitude, not a max") arriving a second time from a different direction.

    ### The divergence residual is not indifferent to which way the cells turn

    `div(rho u)/rho` rms is **4-5 % worse with the knob on in both trees**, same sign and same
    magnitude as ATHAD_PERID's 3.5 %. Three trees, three consistent measurements. That is small
    beside item 80's 2.8x, but it is not zero, and it says the non-closure of items 68/72 has some
    dependence on the prescribed cell parity — a counter-rotating stack is a harder projection
    problem than a co-rotating one. It does not change the conclusion that the residual is
    structural; it adds a term to the list of things it responds to.

    ### The runs are otherwise clean

    Enthalpy drift 8e-16 (off) and 1e-15 (on) in ATHAD's convective adjustment, water drift by
    level identical to three figures between COND's arms, `p_dyn_cap` clamping in 0 of 2 506 179
    fluid cells in every arm. No arm diverged, and the band counts at iteration 20 and 40 agree,
    so what is reported is a settled initial structure rather than a transient sampled once.

    ### Status

    **Default ON stays, on measured grounds in all three trees now.** The banner "MEASURED IN
    ATHAD_PERID / UNMEASURED here" on `073e4ac` and `ff8fd94` is superseded by this item.
    `ATM_CELL_ALTERNATE=0` restores the co-rotating stack exactly. **What is NOT claimed**: that
    the model sustains the alternating cells. 40 iterations is a settled initial structure, and
    item 18 puts geostrophic adjustment ~1e4 iterations away — the question of whether five cells,
    or four, survive to a converged state is untouched by this measurement.

82. **`ATM_GRID_PRESSURE` IS DEFAULT-ON AND MEASURED NOW, AND THE BRANCH'S OWN PRIOR — "expect
    no payoff here" — IS HALF WRONG. It buys two things it never claimed: the `exp_rm` metric
    spread falls 11.77x -> 2.19x and `div(rho u)/rho` falls 24.7 %. It costs two it never
    warned about: the photosphere stops being resolved at all (79 % -> 100 % of columns
    radiating from the prescribed lid) and the near-surface band structure halves. The OLR does
    not move, for the sixth time.**

    `ae25585` flipped the default with an explicit "MEASUREMENT STATUS: NONE IN THIS TREE"
    banner, the same one `073e4ac` carried for `ATM_CELL_ALTERNATE`. This discharges it.
    2 arms x 40 iterations, 24 threads, one binary, arms separated by the environment variable
    alone (`config_pgrid_legacy.xml` / `config_pgrid_new.xml`). **5 min 16 s legacy against
    4 min 56 s on the pressure grid** — 6 % cheaper, not dearer.

    ### The off-branch null is verified against a binary that predates the flip

    `ATM_GRID_PRESSURE=0` reproduces item 81's on-arm — a run made yesterday, before this branch
    was default — to every printed digit of the OLR, `div(rho u)/rho`, max `v`, max `w` and max
    `u`, with a single 8th-digit difference in min `u` (-0.318327 against -0.318326). That is
    item 61's documented reassociation noise. **The legacy path is untouched.**

    ### The grid, measured rather than quoted

    | | legacy | pressure grid |
    |---|---|---|
    | bottom layer `dz_0` | 1224.3 m | **3274.0 m** (2.67x coarser) |
    | top layer | 22.8 km | **10.6 km** (2.15x finer) |
    | lid | 300.0 km | 293.4 km |
    | `checkRadialMetric` spread | **11.77x** | **2.19x** |

    The first three reproduce the comment block's predictions exactly — 2.67x, and a lid that
    barely moves, so in ATHAD this is a redistribution and not a shell cut. **The fourth was not
    predicted by anything.**

    ### The unpredicted win: this grid makes the legacy metric nearly correct

    `exp_rm = 1/(rm+1)` is the Jacobian of a QUADRATIC stretch (item 39), and the model runs an
    exponential-in-height grid, which is why `checkRadialMetric()` — unit-free, and written in
    item 39 before any fix existed — reports an 11.8x spread. **The ln-p ladder placed on this
    column happens to sit far closer to quadratic-in-index than the exponential-in-height grid
    does**, so the same wrong formula is 5.4x less wrong on it. Core length unit 25 125 -> 295 650 m
    legacy, **66 849 -> 146 147 m** on the pressure grid.

    That is a **third route** to item 39's open decision, which until now had two: replace
    `exp_rm` (item 80: done, correct, and default-off because it triples the divergence) or cut
    `zeta` so the wrong formula becomes accidentally right. Changing the GRID is the third, and
    it is the only one of the three that improves the divergence rather than worsening it.

    ### And the divergence residual moves, which almost nothing does

    `div(rho u)/rho` rms **2.884e-02 -> 2.172e-02, a 24.7 % fall**, max 2.517e-01 -> 1.472e-01.
    Item 72 established that this residual is a converged fixed point of the projection — 64x the
    sweeps changes it by nothing — so the list of things that move it downward is short:
    `ATM_PROJ_SWEEPS` (item 68, 52.5 %, plateauing) and now the grid.

    **THIS SITS IN TENSION WITH ITEM 80 AND THE TENSION IS NOT RESOLVED HERE.** Item 80 made the
    metric EXACTLY right and the divergence got **2.8x worse**; this makes the metric 5.4x LESS
    WRONG and the divergence gets 25 % **better**. Both cannot be a simple story about metric
    accuracy. The difference is that item 80 changed the formula on a fixed grid, making the
    radial coefficient `exp_2_rm*inv_dr2` span ~1200x and the elliptic problem far stiffer, while
    this changes the grid under a fixed formula and makes the spacing MORE uniform in the core's
    own coordinate. **The lead that survives both is item 72's: the discrete div and grad are not
    adjoint on this collocated stencil, and grid uniformity — not Jacobian correctness — is what
    that non-adjointness responds to.** Rhie-Chow remains the named un-done work.

    ### What it costs: the photosphere stops existing

    | | legacy | pressure grid |
    |---|---|---|
    | OLR | 135.52 | 135.56 W/m2 (+0.03 %) |
    | imbalance | 135.50 | 135.56 |
    | mean albedo | 0.4990 | 0.4986 |
    | photosphere | 281.3 km, 225.74 K | 272.1 km, **221.12 K** |
    | emission from isothermal skin | 79.0 % | **100.0 % of columns** |

    221.12 K is `t_skin` to two decimals. **Every column now radiates from the prescribed lid**,
    so the model's one instrument for where the atmosphere actually emits reads nothing but the
    boundary condition it was given. Against item 75, which pushed this from 1.1 % to 79 %, this
    finishes the job. The OLR is unmoved to 0.03 %, which is items 29/41/60/80/81 confirmed from
    a sixth direction and is not interesting; the skin fraction is.

    ### What it costs: the near-surface bands

    Psi sign bands, iteration 40 (`survival.py`; note the sampled HEIGHTS differ between arms
    because the grid moved, so only level 0 is a like-for-like comparison):

        legacy   5 5 5 5 5 5 5 4 4      (0 / 5 / 13 / 23 / 36 / 55 / 79 / 113 / 158 km)
        new      2 5 5 5 5 4 4 4        (0 / 14 / 33 / 56 / 85 / 118 / 155 / 193 km)

    **At the ground the count halves, 5 -> 2**, which is what a 2.67x thicker bottom layer should
    do to structure that lives in the bottom few kilometres. Aloft the two are comparable.
    Psi interior max 5.397e13 -> 4.204e13 (-22 %), RMS Psi(ground) 1.350e13 -> 1.222e13 (-9.5 %),
    and the closure ratio therefore gets **worse**, 0.2501 -> 0.2906: the spurious surface flux
    fell less than the circulation it is measured against.

    ### PRECIPITATION APPEARS FOR THE FIRST TIME IN THIS MODEL — AND IT IS SITTING ON THE CAP

    `max precipitation total` goes **0.000000 -> 712.8 mm/d**, at 85 N and 240 km. Decomposed:

        rain     194.4 mm/d
        snow     259.200000 mm/d     <- P_max_flux = 3.0e-3 kg/(m2 s) = 259.2 mm/d EXACTLY
        graupel  259.200000 mm/d     <- the same cap, to six decimals

    **Two of the three categories are pinned at `P_max_flux`, so 712.8 mm/d is the cap's number
    and not the model's.** Item 76 said zero precipitation is a computed result rather than an
    absence, and a capped number does not refute it. What IS new and does not depend on the cap:
    production is nonzero at iteration 40 at all, where item 77 measured the cap binding on the
    FIRST ice-scheme call and never again. Cloud ice 8.13 -> 19.86 g/kg (2.4x) and cloud water
    40.07 -> 47.54 g/kg go with it. **The follow-up is `ATM_PRECIP_CAP` scaled up to find out
    what the scheme actually wants**, and until that is run no precipitation rate from this
    branch should be quoted. This is item 52's lesson arriving on schedule: look at what the caps
    are holding back.

    ### A CORRECTION TO ITEM 79, FOUND HERE: 794.7 g/kg IS A CEILING, NOT A TREND

    `max water vapour` reads **794.700000 g/kg** on this branch, six decimal zeros. It is
    `c_ceiling = 1 - q_CO2` at `UtilsAtm.h:266`, and this model's CO2 is well mixed at
    **0.205300 kg/kg**, so the ceiling is **0.794700** exactly. Item 79 recorded "max q_v =
    794.7 g/kg, **still climbing**" as an open risk of the free-running moist column. **It was not
    climbing; it had arrived at a clamp**, and the same number appearing here on a different
    branch is what gave it away.

    Worse, the instrument denies it: the census line prints `deleted by the c ceiling so far
    0.000000 kg/kg` in both arms, because `UtilsAtm.h:267` clamps silently while
    `cAtmosphereModel.cpp:2031` is the site that counts what it removes. **A field pinned to a
    ceiling by one clamp, with a second clamp's counter reporting zero deletions.** Water is
    conserved to -0.0008 % either way, so nothing is being lost — the vapour is being held, not
    deleted — but "the ceiling deleted nothing" must not be read as "no ceiling is binding".

    ### Status

    **Default ON, on instruction, and now measured rather than asserted.** `ATM_GRID_PRESSURE=0`
    restores the legacy grid exactly. **Every figure in README.md and CLAUDE.md recorded before
    2026-08-25 belongs to the legacy grid.** The open knob is `ATM_GRID_BETA`: ~4.33 restores
    ATHAD's legacy near-surface spacing, which is where both costs above live, while presumably
    giving back some of the metric and divergence gains, which are gains of UNIFORMITY. That
    trade has not been measured and is the obvious next arm.

83. **`ATM_GRID_BETA` DEFAULTS TO 4.33, NOT `zeta`, AND IT IS A TUNED CONSTANT. It recovers
    EVERYTHING item 82's coarse bottom cost — ground Psi bands 2 -> 5, Psi interior back to
    within 0.8 % of legacy, and the best closure ratio of the three grids — while keeping 60 %
    of the divergence gain. AND IT SPLITS ITEM 82's TWO COSTS APART: the near-surface one is a
    beta artefact, the unresolved photosphere is NOT, and no beta will fix it.**

    Item 82 named this arm as the obvious next one and did not predict its outcome. One run,
    40 iterations, 24 threads, `ATM_GRID_BETA=4.33` against the two grids already measured.
    ~5.5 min, water conserved to **-0.0000 %**.

    ### Three grids, one table

    | | legacy | beta = zeta = 3 | **beta = 4.33** |
    |---|---|---|---|
    | `dz_0` | 1224.3 m | 3274.0 m | **1227.5 m** (legacy to 0.3 %) |
    | lid | 300.0 km | 293.4 km | 293.4 km |
    | top layer | 22.8 km | 10.6 km | 14.4 km |
    | `checkRadialMetric` spread | 11.77x | **2.19x** | 7.20x |
    | `div(rho u)/rho` rms | 2.884e-02 | **2.172e-02** (-24.7 %) | 2.433e-02 (**-15.6 %**) |
    | OLR | 135.52 | 135.56 | 135.51 W/m2 |
    | photosphere | 281.3 km, 225.74 K | 272.1 km, 221.12 K | 277.0 km, 221.10 K |
    | **emission from skin** | 79.0 % | **100.0 %** | **100.0 %** |
    | Psi interior max | 5.39721e13 | 4.20390e13 (-22 %) | **5.35487e13** (-0.8 %) |
    | RMS Psi(ground) | 1.34996e13 | 1.22171e13 | 1.31337e13 |
    | closure ratio | 0.2501 | 0.2906 | **0.2453**, best of three |
    | ground Psi bands | 5 | **2** | **5** |

    ### The near-surface cost was beta, and it is fully recovered

    `beta = 4.33` was read off the comment block's table as the value making ATHAD's bottom layer
    match the legacy grid, and it does: **1227.5 m against 1224.3 m, 0.3 %.** Everything item 82
    recorded as a cost of the pressure grid at the ground comes back with it — the band count at
    level 0 returns to 5, and Psi's interior maximum returns to within 0.8 % of the legacy value
    after being 22 % down. **The closure ratio is then better than the legacy grid's**, 0.2453
    against 0.2501, which neither of the other two arms managed.

    ### And the trade is not one-for-one, which is the useful part

    Halfway back in metric spread (7.20x against 2.19x and 11.77x), but **60 % of the divergence
    gain survives**: -15.6 % against legacy, where `beta = 3` gave -24.7 %. So restoring the fine
    bottom does NOT hand back the projection improvement in proportion. Whatever the collocated
    stencil's non-adjointness (item 72) responds to, it is not simply the bottom spacing.

    ### WHAT NO BETA RECOVERS, AND THIS IS THE FINDING

    **Emission from the isothermal skin is 100.0 % of columns at BOTH beta values**, against
    79.0 % on the legacy grid, and the photosphere temperature is `t_skin` to two decimals in
    both (221.12, 221.10). Beta moved the bottom layer by 2.67x and the skin fraction by
    **nothing**.

    So item 82's two costs have two different causes, and only one of them was a knob:

    - the near-surface degradation was the coarse bottom layer — a beta artefact, now gone;
    - **the unresolved photosphere is the ln-p PLACEMENT itself** — the lid at 293.4 km and the
      finer top layers putting more levels inside the isothermal skin — and it is not reachable
      from beta.

    That is a sharper statement than item 82 could make with two arms, and it re-points the work:
    the knob to try is **`ATM_GRID_PTOP`**, which is what moves the lid, not `ATM_GRID_BETA`.
    Whether a lid placed differently can pull the tau = 1 crossing back out of the skin is
    unmeasured. **Note what is at stake — with 100 % of columns radiating from the prescribed
    lid, the model has no instrument left for where this atmosphere actually emits**, which is
    invariant 3's problem in its most acute form yet.

    ### Unchanged from item 82

    Precipitation is nonzero and still on its cap: snow and graupel at **259.200000 mm/d** each,
    `P_max_flux` to six decimals, with a differently-located total maximum (583.2 mm/d at 37 N,
    against 712.8 at 85 N) purely because the three category maxima sit at different latitudes.
    `max water vapour` reads **794.700000 g/kg** for the third arm running — the `1 - q_CO2`
    ceiling of item 82's correction to item 79. The OLR is unmoved across all three grids
    (135.51-135.56, a spread of 0.04 %), a seventh confirmation.

    ### Status, stated as bluntly as it deserves

    **`beta` defaults to 4.33 and that is a TUNED CONSTANT, not a law.** `beta = zeta` was the
    principled choice — the same exponential stretch on a different coordinate — and it lost on
    measurement. The value is read off this tree's row of the comment block's table and **differs
    per tree**: ATHAD_COND ~3.78, ATHAD_PERID ~2.83. A sibling porting this branch must re-read
    its own row. Given how much of this README is about Earth's constants surviving as literals
    inside physics kernels, a fitted per-tree constant is written down here as exactly that,
    with the table it came from three lines above it in the source.

84. **THE PHOTOSPHERE IS A CLOUD TOP, THE CLOUD DECK IS PLACED BY THE TEMPERATURE PROFILE AND
    BY NOTHING ELSE, AND FIVE SEPARATE ATTACKS ON IT MOVED NOTHING. The one thing that DOES
    move the emission level off the prescribed lid is turning item 75's correct repair OFF —
    so the repair that made the condensate real is also what buried the photosphere in an
    assumption. And the column is supersaturated by 10^2-10^5 %, on BOTH branches, cause not
    yet found.**

    Item 83 ended by naming `ATM_GRID_PTOP` as the next knob. That was wrong and is withdrawn
    inside this item. The chain of measurements below started from item 83's `dtau` probe and
    ended somewhere else entirely.

    ### First, the finding that reframes everything: the opacity at the top is ~100 % CLOUD

    Decomposing `tau_layer` into its gas and cloud terms on the prognostic column (equator):

    | h [km] | `tau_layer` | `tau_gas` | cloud share |
    |---|---|---|---|
    | 243.9 | 5313 | 0.317 | ~100 % |
    | 254.4 | 2844 | 0.059 | ~100 % |
    | 266.0 | 1099 | 0.008 | ~100 % |
    | 293.4 | 88.2 | 1.5e-5 | ~100 % |

    At 60 N the crossing is starker still: `tau_layer` = 553 at 254.4 km where cloud exists and
    **0.157** at 278.9 km where it does not. The gas contributes between 1e-5 and 4 where the
    total is thousands. **This model's emission level is the top of its cloud deck.**

    Two consequences, one of which killed a planned piece of work before it was written:

    - `buildReferenceColumn` is dry and cloud-free, so a `ln(tau)` ladder built from it would
      place levels where the GAS goes transparent — six orders of magnitude from where the
      atmosphere actually does.
    - And for the gas term that ladder is not even a new grid. `dtau = kappa*(dp/g)*(p/p_ref)`
      gives `tau ~ p^2`, so `ln tau = 2 ln p + const`: uniform in `ln(tau_gas)` IS
      `ATM_GRID_PRESSURE` with `Lambda` doubled. Analytic, no run needed.

    ### `ATM_GRID_TAU`: the two-pass rebuild, written, measured, and with a ceiling

    Since the only usable `tau` is the one the model computes, the ladder is built in two
    passes: run the radiation once on the initialised state, take the cos-lat-weighted mean
    `tau_above` it produces (clouds included), rebuild levels 0..im-2 uniform in `ln tau`,
    keep level im-1 at the old lid so the domain and `p_top` do not move, re-interpolate the
    29 `restart_arrays()` fields, rebuild the metric table, redo `init_tropopause_layers` and
    call `densities()`. Static thereafter. `ATM_GRID_TAU_MIX` blends against the grid
    `init_layer_heights` built; both sequences are strictly increasing so any convex
    combination is valid.

    | mix | `dz_0` | metric spread |
    |---|---|---|
    | 0 (base) | 1227.5 m | 7.20x |
    | 0.25 | 3040 m | **4.51x** |
    | 0.5 | 4853 m | 5.35x |
    | 1.0 | 8478 m | 31.98x |

    Pure `ln(tau)` gives a per-layer tau ratio of **0.728, constant by construction** — `dtau`
    ~ 0.3 at the crossing, which was the entire point — but starves the deep column, where
    99.99 % of the mass and all of the dynamics live and the radiation has nothing to resolve.

    **Measured at mix = 0.5, 40 iterations**, layer `dtau` at the equator:

        no regrid   15.32   63.63   390.7   1502
        regrid      23.46   43.79   81.06   202.5

    The upper-middle column improves — per-layer ratio 4-6x down to 2.2-2.9x — the top layer
    gets WORSE, and **`tau = 1` is still inside it**. `div(rho u)/rho` rms goes 2.433e-02 ->
    **4.076e-02 (+68 %)**; OLR 135.55, `skin%` 100.0, both unmoved.

    **THE CEILING IS THE DESIGN'S OWN.** The ladder is built from the `tau` profile at
    iteration 0, and the photosphere is a CLOUD TOP. Clouds move. By iteration 40 the grid no
    longer matches the profile it was built for — the same property that made the
    reference-column ladder impossible, reappearing one level up. A static rebuild cannot
    track a moving emission level. **Default off.**

    **A bug worth recording, because the guard is what caught it.** The height inversion used a
    two-sided bracket test; at `i = im-2` the target IS `ln(tau[im-2])` and rounding put it a
    hair below, so nothing matched, `lo` stayed 0 and the interpolation extrapolated 272 layers
    past the lid. The monotonicity check refused the grid and printed why, rather than shipping
    it. It now takes the last level at or above the target and clamps the weight.

    ### The noise floor, measured, because the next results are small

    The off-branch null was NOT bit-identical, so before reading anything the null arm was run
    **twice on the same binary at the same 24 threads**:

        same binary, twice ....... min q_v 554.746 vs 555.429 g/kg   (0.12 %)
        old binary vs new ........ min q_v 554.746130 vs 554.746131 (9th digit)

    **Run-to-run scatter exceeds the binary difference**, so `ATM_GRID_TAU` off is inert. This
    extends item 61: that item says fixed-thread run-to-run divergence shows "nothing in the
    printed scalars", and at 40 iterations with moist physics it demonstrably does — 0.12 % in
    min water vapour. Item 61's consequence stands and is now better evidenced. **Every null
    below is judged against this 0.12 % floor.**

    ### Attack 1 — the precipitation cap. NOT the blocker.

    | `ATM_PRECIP_CAP` | x1 | x10 | x100 |
    |---|---|---|---|
    | max precipitation | 583.2 | 5372.2 | **42 898.7 mm/d** |
    | max cloud water | 43.660 | 43.713 | **43.796 g/kg** |
    | max cloud ice | 21.671 | 21.627 | 21.625 g/kg |
    | max `q_v` | 794.700000 | 794.700000 | 794.700000 |
    | photosphere | 277.0 | 277.1 | 277.1 km |
    | OLR | 135.51 | 135.51 | 135.51 W/m2 |

    The cap WAS binding — item 82 measured snow and graupel pinned at `P_max_flux` to six
    decimals — and releasing it raises the through-flux **74x** while the deck stays within
    0.5 % at every level. **What moves is the opposite of draining**: layer-mean `q_v` just
    BELOW the deck rises, 692.9 -> 717.7 g/kg at 224.8 km and 710.3 -> 763.7 at 234.5 km. The
    precipitation falls out of the deck and evaporates straight back into vapour underneath it,
    in the supercritical column — item 76's `evaporateWhereImpossible`. **The deck sits on top
    of its own recycling loop**, so opening the drain wider only runs the loop faster. The deck
    is not drainage-limited.

    ### Attack 2 — `damp_wiggles` on the moisture fields. NOT the blocker either.

    Item 74 recorded this smoother taking cloud-bearing cells from 131 404 to **263 530**,
    spreading condensate into cells `canCondense` forbids, and left it as "not fixed".
    `ATM_DAMP_MOIST=0` skips it for `c`, `cloud` and `ice` at both call sites (`t` is
    deliberately excluded — its own comment records that an undamped 2-delta-t mode once ran it
    to 53 K).

    | | damping on | off |
    |---|---|---|
    | OLR | 135.51 | 135.50 W/m2 |
    | photosphere | 277.0 | 276.0 km |
    | skin % | 100.0 | 100.0 |
    | `div(rho u)/rho` | 2.433e-02 | 2.433e-02 |
    | max cloud water | 43.660 | **49.999580 g/kg** |
    | min water vapour | 554.746 | **0.004498 g/kg** |

    Layer-mean condensate shuffles WITHIN the deck (10.8 -> 8.6 at 234 km, 32.8 -> 26.7 at
    254 km) and the deck does not move, thin, or change the emission level.

    **But the smoother is LOAD-BEARING, not merely sloppy, and item 74's framing should not be
    acted on as though it were a defect.** Without it, max cloud water pins at
    **49.999580 g/kg = `cloud_cap` exactly** and the vapour develops a grid-scale hole at
    **0.0045 g/kg** with neighbours near 700. Removing it drives the field into two caps, not
    into a better state.

    ### What actually places the deck: the temperature profile, and only that

    Relative humidity, layer means, shipped arm:

    | i | h [km] | T [C] | mean RH | max RH | cloud |
    |---|---|---|---|---|---|
    | 33 | 214.2 | 149.4 | 18.5 % | 52.6 % | 0 |
    | 34 | 224.8 | 85.8 | **134 %** | 636 % | 0 |
    | 35 | 234.5 | 27.6 | **6 004 %** | 62 180 % | 10.77 |
    | 36 | 243.9 | -31.2 | **49 580 %** | 101 800 % | 35.04 |

    **The deck's base is exactly where RH crosses 100 %**, between 214 and 225 km. The
    condensate is placed by where the temperature profile first permits condensation — which is
    why the cap, the smoother, the grid, and before them `kappa` at 64x (item 29) and the
    circulation at 500x (item 27) all failed to move it.

    The instrument was checked before the number was believed: `HumidityRel`
    (`ThermoAtm.h:402`) is `e/E*100` from `SaturationH2O`, the same module
    `SaturationAdjustment` uses, uncapped, with the `NO_SATURATION` = -1 sentinel above the
    critical point. The comment there records that the old code capped it at 100 %, which is
    why a 10^5 % supersaturation had never been seen.

    ### Attack 3 — `ATM_SAT_SUPERHEAT`. Exonerated, and it is the one thing that moves the lid.

    The hypothesis was that item 75's re-test retains the supersaturation by rejecting the
    condensation step in cells that condensing would superheat. **Refuted — RH is essentially
    identical on both branches:**

        h [km]      224.8    234.5    243.9    254.4
        superheat ON   134.5    6 004   49 580   30 090 %
        superheat OFF  138.8    6 082   58 690   31 220 %

    What the knob DOES control is whether the condensate survives to the next radiation call in
    the top three layers, and the consequence is the uncomfortable result of this item:

    | | ON (shipped) | OFF |
    |---|---|---|
    | `tau_layer` at 266.0 km | 254.3 | **0.018** |
    | at 278.9 km | 31.95 | 0.589 |
    | at 293.4 km | 7.70 | 0.147 |
    | photosphere | 277.0 km, **T = `t_skin`** | **243.0 km, T = 270.46 K** |
    | skin % | **100.0** | **3.3** |
    | OLR | 135.51 | 181.56 W/m2 |
    | max cloud water | 43.660 | 0.000000 (item 74's annihilated state) |

    **Item 75's repair — which is physically correct — is precisely what buries the emission
    level in the prescribed isothermal lid.** Turning it off is the only thing besides
    `ATM_PROGNOSTIC_T` that has put the photosphere back on the real adiabat at a real
    temperature. It does so by keeping condensate in cells that cannot hold it, which is exactly
    what item 75 diagnosed, so **this is not a recommendation to flip the default** — it is a
    measurement of where the inconsistency lives. `albedo` barely moves either way (0.4990 vs
    0.4982), consistent with item 74's point that the radiation sees condensate mid-iteration.

    ### Attack 4 — the Newton loop's iteration count. NOT that either.

    `satIters()` is documented in its own comment as an Earth constant, with a measurement from
    ATHAD_COND showing cells left at 15-99.8 % of their own target. The loop is a DAMPED Newton,
    residual ~ (1-omega)^n, so 20 passes could plausibly close nothing. `ATM_SAT_ITERS=200`,
    6m19s against the control's 5m08s (+23 %, so the loop is not a large share of the runtime):

        h [km]      214.2   234.5   243.9   254.4   266.0
        n = 20      18.47    6 004  49 580  30 090   7 717 %
        n = 200     18.47    5 739  50 150  30 290   7 771 %

    **Unchanged.** The adjustment has CONVERGED — to a state with RH = 30 000 %. Max cloud water
    moves 43.66 -> 37.38 g/kg and the photosphere, `skin%` and `tau_layer` do not move at all.
    So the loop reaches its target and **the target is not saturation.**

    ### AND THE ANSWER IS ITEM 78's, QUANTIFIED: the prescribed profile throws away the latent heat

    Same levels, shipped prescribed arm against `ATM_PROGNOSTIC_T=1` (400 iterations, same grid):

    | h [km] | T presc. | T prog. | RH presc. | RH prog. |
    |---|---|---|---|---|
    | 234.5 | 27.6 C | 63.6 | 6 004 % | **91.6 %** |
    | 243.9 | -31.2 | **36.5** | 49 580 % | **112.3 %** |
    | 254.4 | -52.0 | 21.6 | 30 090 % | **107.1 %** |
    | 266.0 | -52.0 | 6.2 | 7 717 % | **104.1 %** |

    **The free-running column is SATURATED, 92-112 %, exactly as a cloudy layer must be — and it
    is 68 K warmer at 243.9 km.** That is item 78's finding measured in the humidity field:
    condensation releases ~60 K, `densities()` overwrites `t` with the adiabat at the end of the
    iteration, the cell is cold again on the next pass, and it is supersaturated again. **The
    prescribed profile REGENERATES the supersaturation every iteration.**

    So the 10^2-10^5 % is not a microphysics failure at all. Four separate microphysics
    suspects were excluded before the profile was — the cap, the smoother, the superheat
    re-test, the iteration count — and every one of those nulls is really the same null seen
    from a different side. **The condensate amount, the vapour ceiling and the deck's opacity
    are all downstream of invariant 3.**

    Note what survives: the DECK ITSELF is similar on both branches (layer-mean cloud water
    10.77/35.04/32.82 prescribed against 12.36/31.55/30.44 prognostic). The deck is real. Only
    the humidity around it is an artefact.

    ### What this item changes

    - **Item 83's "the knob to try is `ATM_GRID_PTOP`" is WITHDRAWN.** The lid is not the lever;
      the emission level is a cloud top and the cloud top is placed by the temperature profile.
    - **Item 74's framing of `damp_wiggles` as a defect should not be acted on.** It is
      suppressing a real grid-scale instability in the moisture fields.
    - **Item 75's repair is exonerated as the cause of the supersaturation** and simultaneously
      identified as what buries the photosphere in the prescribed lid. Both, at once.
    - **The RH field is a first-class instrument now and was not being read.** Its own comment
      records that the old code capped it at 100 %, which is why a 500x supersaturation sat
      unnoticed through items 74-83.
    - **`ATM_GRID_TAU` ships default-off** with its ceiling documented: a static rebuild cannot
      track an emission level that is a moving cloud top.

    **The open question is unchanged in kind and much sharper in form.** With
    `ATM_SAT_SUPERHEAT` on and the profile prescribed, 100 % of columns radiate from a lid whose
    temperature is an assumption. Turning the repair off moves the photosphere onto the real
    adiabat by keeping impossible condensate. Freeing the profile fixes the humidity and the
    photosphere together — and item 79's prognostic arm was still rising at 400 iterations. That
    is where the work goes; not the grid, and not the microphysics.

## Remaining work

- **The prescribed adiabat and the grey opacity are incompatible, and that is now the radiative
  question** (item 73). Item 67's one-line follow-up is written, measured and withdrawn:
  `T_rad(tau)` exceeds the prescribed profile at all 41 levels, so there is no
  radiative-convective crossing to switch at, and `tau_s` = 2.26e6 at `F/2` = 136 W/m2 wants a
  9496 K surface against the model's 1408 K. Meanwhile the shipped model's OLR converges to
  `sigma*t_skin^4` = `F/2` to 0.01 W/m2, so its imbalance is definitionally half the energy
  input. **Nothing radiative moves until one of the three gives — the prescribed surface
  temperature, `kappa`, or the greyness** — and a top-boundary change is not one of the three.
  `ATM_SKIN_TAU` (default 0) is kept as the record of what was tried.
  **And do not re-run the kappa scan on the prescribed branch to find out** — item 29 did it
  (0.10 % over 64x), item 60 confirmed it by another species (0.083 %), and the null is now
  structural rather than empirical: `t_skin` is a fixed point of a budget with no kappa in it,
  and item 73 measured the OLR sitting on `sigma*t_skin^4` to 0.01 W/m2, which is *tighter*
  than the branch those scans were run on. The kappa arm still worth its wall clock is under
  `ATM_PROGNOSTIC_T=1 ATM_RAD_DIRECT=1` (item 66's configuration, `skin%` = 0.0 throughout),
  where no kappa arm has ever been run and the baseline has to be re-run beside it because
  item 66 predates items 67-73.

- **Rhie-Chow face reconstruction is the largest piece of unbuilt machinery here** (item 72).
  The pressure projection has CONVERGED to a fixed point that is not divergence-free —
  `div(rho u)/rho` rms 1.611e-02 is bit-identical under 64x the sweeps — so the discrete
  divergence and gradient are not adjoint on this collocated stencil. Five mechanisms are
  measured and excluded (base-state density 1.0 %, `exp_rm` worse when exact, no checkerboard in
  velocity or pressure, not a boundary artefact). `ATM_PROJ_SWEEPS = 10` is a knee, not a cure,
  and item 68's 44 % residual is this. **It is a numerical-methods port, not a constant to fix**,
  and the residual is NOT attributed — the checkerboard the collocated story predicts is absent.

- **The updraft is one grid level deep, and that is now the top microphysics job** (item 61).
  Cloud base 236.4 km, LFS 256.0 km — adjacent levels in 99 of 173 convecting columns — so the
  recurrence `for(i = i_base+1; i <= i_LFS-1)` is an empty loop, no parcel ascends, and `c_u` is
  identically zero for a reason that has nothing to do with the thermodynamics. The
  deep-convection triggers (1000/970/900/800 hPa, absolute Earth surface pressures) are the
  suspect; they must become fractions of the local surface pressure. **`ATM_MC_GEOPOTENTIAL` is
  written and waiting for this**, and should be flipped after it, not before.
- **The 40-iteration OLR is reproducible to 0.51 % across a thread-count change** (item 61),
  which is 5–25× every OLR effect this file has reported. Nothing is retracted — those were
  controlled same-binary same-thread A/Bs — but **quote no OLR difference below ~0.5 % as a
  property of the model**, and state the thread count with any that is quoted. The amplifier is
  the convective trigger set, not the dynamics: Ψ, albedo and the photosphere are unmoved at
  1e-6 in every null pair. Curing it is the same ordered-reduction job item 18 names.


- **The background's opacity was N₂'s, applied to NH₃, CH₄ and SO₂ as well** (item 60). Split
  per species it is **1867× larger** and nearly 2× `kappa_CO2`, with NH₃ (49 %), SO₂ (37 %) and
  CH₄ (14 %) carrying all of it. The converged OLR moves **0.083 %** — items 29 and 41 confirmed
  a third time, by a different species and a different mechanism. The six κ are assumptions;
  `ATM_BG_LUMPED=1` restores the old value. **What is still open is the ordering**: NH₃ and CH₄
  being far from inert is certain, their individual κ are not, and a grey scheme cannot carry
  the window regions that would decide it.
- **~~The composition closes on the background, so CO₂ can never dilute~~ — FIXED in item 57**
  (`ATM_CO2_DILUTE=0` restores it). The CO₂ column gains 320× more spatial structure, the global
  CO₂ mass is conserved to 5.6e-07, and the OLR moves 0.012 %. What remains is second order: `c`
  is vapour only — **and that is fixed too, in item 59**, which also deleted `densities()`'s
  partial `water_factor` so the condensate is not counted twice. Original entry:
- **The composition closes on the background, so CO₂ can never dilute** (item 56). `split()`
  sets `q_b = 1 − q_v − q_c`, and water — 67 % of the mass — varies 536.9 to 739.0 g/kg across
  the domain, so the whole of that variation is absorbed by the background fraction (a **4.6×**
  spread) while CO₂ stays pinned at 0.205300 everywhere. Concentrating CO₂ and background
  together, which is what removing water from a parcel does, moves `R_mix` **−1.4 % at the
  wettest point and +2.6 % at the driest** — a systematic, water-shaped bias in `R_mix`, `cp_mix`,
  the density and the hydrostatic column. This is also the reason CO₂ shows no convective
  influence: the other two reasons (no `MC_co2` term exists; `rhs_co2 ≡ 0` for a uniform field
  with no source) are consequences of there being no CO₂ source, but this one is a choice.
- **`N²` measures the prescribed profile and the integrator, not the atmosphere** (item 55).
  It is formed at the end of `densities()` from the `t` and `p_stat` that same call just wrote,
  and it is identical to 0.03-0.5 % across runs whose `p_dyn` radial structure differs by
  **2500×**. Below the skin its value is a **second-order truncation error** of the adiabat
  integration — `N²/dz²` is constant to ±30 % over a 12× range of layer thickness while `N²/dz`
  varies 16× — so the monotonic rise up the column is the grid stretch, not stratification. The
  skin value is genuine (`g²/(cp·T)`, implying cp = 1543 against `cp_of`'s 1563, 1.3 %). It
  becomes a measurement of the atmosphere only under `ATM_PROGNOSTIC_T=1`.
- **`p_dyn` is 99.97 % a prescribed balance with no radial structure, and the two off-by-default
  metric terms are what supply one** (item 54). Turning `ATOM_CORIOLIS_NONTRAD` and
  `ATOM_METRIC_CURVATURE` on multiplies the radial range of `p_dyn` by **192× at j=45**, at
  iteration 0, and turns the radial momentum budget from one unopposed pressure gradient
  (rms 0.30, nothing against it) into a real balance (cor 18.9 against pgf 18.3, residual 2.2).
  Cost: Ψ_max −4.4 % at iteration 20 and **−9.1 % at 40** — a growing gap, so the endpoint is
  unknown. **Not flipped**, and two things must come first: `ubud_advh` is built from the raw
  advective pieces and does **not** capture the curvature term, so with those switches on the
  radial budget is missing a live term by construction; and a run long enough to say what Ψ does.
- **The updraft cannot condense any more, and that is the `g·z` bill coming due** (item 53).
  With the `s` arithmetic repaired, ATHAD_COND's `max c_u` goes 0.0678 → **0.000000**: the
  baseline's updraft condensation was computed from a parcel at 1128 K that the `/s_0` defect
  invented, and the correct 509 K parcel is sub-saturated. It will stay sub-saturated, because
  `s` is `cp·T` with no geopotential, so the parcel does not cool as it rises (item 52). ATHAD's
  `c_u` was already identically zero. **Adding `g·z` to `s` is now the top microphysics job** —
  it is no longer one defect among four, it is the only thing between this scheme and an updraft
  that condenses.
- **`MCt_max` and the 0.01 kg/m³ density floor are the next two Earth constants to size**
  (item 53). With the updraft denominator made consistent, `max MC_t` reaches its 0.01 K/s cap
  at iteration 10 — attributed to the denominator by `ATM_MC_UNBOUNDED_UPDRAFT=1`, which returns
  the pre-repair 8.3e-5 exactly. The cap is 36 K/hr, annotated "still 3× realistic" for Earth,
  and `inv_step_rh` floors the density at "≈ air at 50 hPa". At 256 km in a 250 bar column
  neither number was sized for what it is now bounding. Iterations 20 and 40 agree to 0.3 %
  across every build, so nothing downstream depends on it yet.
- **The moist-convection `s` fields are a normalised temperature wearing two other names, and
  three scaling errors sit on top of that** (item 52). `ATM_MC_S_CONSISTENT` repairs the three
  and is **default off** pending a longer run than the 40 iterations measured. The big one is
  that `MC_t` sat **pinned at its `MCt_max` cap** through the spin-up while updraft parcels ran
  at **46 000 K** — a safety cap absorbing an arithmetic error, which is the worst way for one
  to hide. **Open and not in the knob**: the missing `g·z` (dry static energy is `cp·T + g·z`
  and `g·z` is the *larger* term above ~157 km here, so the updraft never cools as it rises),
  `MC_t`'s second term carrying an extra `t_0`, `q_v_u` reaching 3.5 kg/kg, and the cubic lid
  extrapolation that gives `s` negative values — a negative absolute temperature — at 300 km.
- **~~Three microphysical corrections have now died at the same wall, and the wall is the
  albedo~~ — THE WALL WAS AN ANNIHILATION, REFUTED BY ITEM 75** (item 51). Making `cp` local in
  `MoistConvection` moves rain and snow production by
  **−18 %** and the OLR by **+0.02 %**, with Ψ, the photosphere and the albedo bit-identical or
  unchanged; `initCloudIce`'s H_crit repair was 0.18 % (item 42). ~~The reflectivity saturates on
  the **presence** of condensate, so nothing about condensate *amount or rate* can reach the
  radiation.~~ ~~**`albedo_cloud` is not merely the biggest lever after the opacities — it is the
  only path condensate has to the OLR**, and any further microphysics work should expect to be
  unmeasurable until that parameterisation responds to something.~~
  **REFUTED BY ITEM 75.** The ALBEDO half is right and the conclusion drawn from it was wrong:
  condensate reaches the OLR through the LONG-WAVE cloud opacity `k_liq*LWP + k_ice*IWP`, which
  responds to amount. That path looked dead because `SaturationAdjustment` was manufacturing
  condensate into superheated cells and the ice scheme was deleting all of it every iteration
  (item 74) — so the five repairs that "died at the albedo" were measuring an annihilation, not a
  wall. With the condensate alive, one repair moved the OLR **−49.4 %** and the photosphere 44 km.
  The lesson is not that the albedo is the only lever; it is that **five null results in a row
  should have prompted a check that the field being varied still existed downstream.**
- **`s_0` carries Earth's dry-air cp and is dead-code-protected** (item 51). `s_0` = 274515.75
  = 1005 × `t_0`, not ATHAD's 2040 × `t_0` = 557226, while `param.py` calls it "`cp_l * t_0`".
  It changes no result because `cp_l` cancels in the `s ↔ T` pair, and the one place the scale
  would matter (`s` = 1.0 at `BC_Atm.h:219,368`) is behind `if(!is_land) continue`, dead under
  invariant 1. Fixing `s_0` and the `s` definition together is a separate job with no measurable
  payoff — recorded so it is not rediscovered as a live defect.
- **Drag is eliminated as the cell-decay driver, and the timescale argument that eliminates it
  disqualifies more than drag** (item 49). A correctly-scaled surface drag explains 0.034 % of
  the Ψ decay. The reason is that 200 iterations is 39 s – 12.5 min of physical time (item 47),
  so **nothing with a timescale longer than minutes can explain a decay that completes inside
  200 iterations** — that removes drag, radiative relaxation and surface exchange together. The
  remaining search is the initialisation transient and geostrophic adjustment (items 26–28).
  **The attribution is the open part**: `vbud_*`/`wbud_*` are written to CSV and have not been
  analysed, and Ψ is built from `v` alone so it cannot see a radial failure (item 28).
- **`t_skin` blocks three separate lines of work, and that is now the top priority** (item 41).
  The OLR equals `σT_lid⁴` equals absorbed SW + geothermal, exactly, and it has now failed to
  respond to κ (64×, 0.10 %), the circulation (500×, 0.03 %), `im` (analytically) and the grid
  stretch (6×, 0.36 % and that albedo-driven). Item 29's opacity question, item 39's photosphere
  question and item 41's grid question are all unanswerable until the fixed point is broken.
  **Measuring anything radiative before then is spending runs to re-measure `t_skin`.**
- **The radial metric is wrong in shape and `zeta` is the lever, not `im`** (item 39).
  `exp_rm = 1/(rm+1)` is a quadratic-stretch Jacobian applied to an exponential grid, so every
  radial derivative is mis-scaled by a factor varying **11.8×** across the column — measured,
  unit-free, and printed at every startup by `checkRadialMetric()`. Of the three arguments for
  cutting `zeta` hard toward ~0.4, item 41 leaves **one standing and one untested**: the metric
  itself stands, joined by the new finding that only `zeta ≤ 1.5` reaches its own energy fixed
  point within 200 iterations; the photosphere argument is unconfirmed (the OLR cannot respond
  while `t_skin` pins it) and the diffusive-CFL argument is untested (blocked behind the `L_atm`
  normalisation below).
- **The `L_atm` normalisation blocks the CFL half of the `zeta` question** (item 41). The
  Held–Suarez relaxation, Rayleigh drag, `coeff_MC_*`, `coeff_S`, `coeff_L`, `nue_max` and
  `coeff_u_p` all divide by `L_atm`, which a shell-preserving `zeta` scan moves by 29×, so a
  raised-`dt` run cannot attribute a blow-up. `RHS_Atm_Turb.cpp:441` already documents these as
  wrong and names the fix — move them onto `metricShellLength()`, as `force_nd` already was.
  Same category error as item 39: a stretch amplitude used as a grid length.
  **The open decision is which repair**: replacing `exp_rm` with the true Jacobian is
  principled but touches ~91 sites and moves every number here, while cutting `zeta` makes the
  existing metric accidentally correct for ~30× less wall clock. That wants a measurement, and
  it should come before any further `im` or shell-depth work, because both are being chosen
  against a grid the dynamics and the radiation do not agree on. **Worse upstream**:
  `ATOM_Precipitation` at `zeta = 3.715` has a 23.2× spread.
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
- **The tropical cell decays even when the balance is right** (item 28), and item 31 split
  that decay in two: about a third was Earth's imposed cell latitude, and the rest tracks the
  unbalanced θ-force, of which mode 2 removes 61 %. At the regime's own latitude the decay
  is −7.3 % over 200 iterations and no longer improves with narrowing. Removing the residual
  force needs the rotational part balanced by buoyancy — a thermal-wind temperature
  perturbation — which `buoyancy_ramp` = 0 at iteration 0 and `densities()`'s overwrite of
  `t` both block. **Same task as the prescribed profile above.** `ATOM_METRIC_CURVATURE` was
  tested as the alternative and changes nothing (0.02–0.4 %).
- **`cell_lat_scale` is now a config parameter and its default is 0.33, not Earth's 1.0**
  (item 32), and it now scales the **Hadley edge only** (item 36). The Ro_T derivation picks
  the value; item 31's measurement corroborated it but cannot discriminate 0.33 from 0.50,
  and item 36 showed that scan ran on a layout where narrowing the cells mostly *deleted*
  the outer ones. Every result in this file from before item 32 was measured at 1.0; set
  `<cell_lat_scale>1.0</cell_lat_scale>` to reproduce them.
- **`n_cells_hemisphere` now defaults to 5, not Earth's 3** (item 38), on three grounds from
  items 36–37: the extratropical bands match the Rhines scale (20° against ~22°), the
  geometric artefact of the 40°-wide cells is gone, and at 400 iterations n = 5's Ψ maxima
  still sit exactly on its prescribed cores while n = 3's migrate equatorward off theirs.
  **None of the three is a converged measurement**, and the interior-boundary argument of
  item 36 is withdrawn outright — it reversed by 400. Everything measured before item 38 used
  3; set `<n_cells_hemisphere>3</n_cells_hemisphere>` to reproduce it.
- **The tropical cell decays linearly to 400 iterations with no turnover** (item 37),
  −16.7 % at n = 3 and −16.2 % at n = 5, converging to within 0.35 %. Cell count does not
  reach the Hadley cell at any run length tested.
- **The N–S asymmetry grows with run length** (item 37): ~5e-08 at 20 iterations, −3.7e-06
  at 400. Still 3× inside the 1e-5 tolerance, but it is a trend and invariant 1 is the one
  thing this model is not allowed to break. Worth a longer run purely as a check.
- **Every kinetic-energy figure measured before item 36 was taken with half the circulation
  missing.** Mean KE is 38–40 on the tiled layout against 19.50 on the untiled one. Item 34's
  `u_0` comparison is affected; item 33's mode comparison is internally consistent but its
  absolute numbers are not.
- **`cell_amp_mode` is written and default-off** (item 33). The tropical-cell decay tracks
  the prescribed **zonal jet**, not the overturning: modes 0 and 1 differ 3× in every
  amplitude and decay identically, while mode 2 changes only the jet and halves the decay.
  That is the first handle on item 28's residual force.
- **The buoyancy body force still carries an extra `*dt`, and now it also has a knob**
  (items 34, 50). The drag half was fixed in `3e2d78f`; this half is `ATM_BUOY_CONSISTENT`
  and is **default off**, because the consistent coefficient is **2e7×** the shipped one and
  **7.2× larger than the 336 that is recorded as having driven a polar vertical runaway** — on
  a body force, where the drag's "larger is stabilising" argument does not apply. The knob also
  repairs the Boussinesq reference temperature, which is the bigger surprise: the shipped term
  divides the anomaly by `t_0` = 273.15 K instead of `T_ref`, overstating the force **5.49× at
  the surface and 0.96× at the top** — a height-dependent distortion, not a rescaling. The open
  work is a long run with it on, and deciding what `buoyancy_ramp` should reach; a 20-iteration
  test at 6.7 % ramp and without the `T_ref` half is **not** evidence of stability.
- **`ATM_RAD_DIRECT` is written, exact and 10× cheaper than the wrong default, and is still
  off** (item 30). The case for flipping: it changes the standard configuration by 0.15 %,
  it is a closed-form solution rather than an under-iterated one, and it costs less. The case
  against: one 200-iteration measurement. `ATM_N_LAMBDA` exists to reproduce the old
  behaviour and to show what it was worth.
- **The prognostic column is relaxing, and where it lands is unknown** (item 30). With the
  converged solver its OLR falls 15 764 → 8 542 → 6 211 over 200 iterations and is still
  falling, against ~66 000 and motionless with 4 sweeps. **Do not quote 6 211 as a converged
  number** — a monotone trend is not a limit, which this file has now been caught on twice.
  A long prognostic run with the direct solver is the single most informative thing to do
  next, and unlike the κ scan it has a real chance of changing invariant 3.
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
- **`dt_visc` is tied to the grid, and raising `im` means lowering it.** The config ships
  1e-4, which is the value item 24 validated at the committed `im = 41`; at `im = 61` the same
  value is ~2.2× over the explicit diffusion limit, which goes as the physical surface spacing
  squared (1.22 km at 41 levels against 0.81 km at 61). *This bullet previously said the config
  still shipped 4e-5 and was stale — the 1e-4 restoration landed with the `im` default. Item 39
  is the live version of this concern: it argues the lever is `zeta`, not `im`, and cutting
  `zeta` loosens this limit rather than tightening it.*
- **The run is not converged, and 400 iterations is not close** (item 18). The meridional
  wind is in free acceleration under an unopposed Coriolis torque — the pressure gradient
  that should balance it is at 1.8 % of it after 400 iterations and growing linearly, which
  puts geostrophic adjustment of order 10⁴ iterations away. Ψ_max grows linearly throughout.
  **This is the largest open question about the dynamics**, it is common to both continuity
  formulations, and item 18 ruled out the elliptic solver as its cause.
- **The OLR is not grid-converged either**: 519 W/m² at a 260 km shell against 581 at
  300 km, with `im` fixed at 61. Both figures predate item 22, so the check has to be redone
  as well as extended — refine vertically and repeat.
- **~~`moist_phys_start_iter = 300`~~ — IT IS 0 SINCE 2026-08-23** (item 74). The 300 was
  ATOM_Precipitation's, inherited on the fork commit, justified by "a single tropical maritime
  column" — an Earth surface classification, in a model with no sea — and waiting for a
  circulation that item 18 puts 1e4 iterations away. Measured at 200 iterations: no runaway.
  **But the moist column does not converge where the dry one did** (OLR 168 and falling against
  136 arrived), so every dry-column number in this file belongs to the old default; set
  `<moist_phys_start_iter>300</moist_phys_start_iter>` to reproduce it. Original entry: it
  means a 400-iteration run is dry for three quarters of
  its length. Deliberate (it lets the circulation form before the stiff microphysics
  starts), but it must be stated whenever a run is quoted — and it is why the 20-iteration
  measurements above are all made with sedimentation and the ice schemes switched off.
- **~~Condensate is created and destroyed every iteration~~ — FIXED in item 75**
  (`ATM_SAT_SUPERHEAT`, default on): `clampAndFade`'s own latent heating invalidated its own
  admissibility test. OLR -49.4 %, cloud cells surviving the ice scheme 6 486 -> 149 271. What
  remains open is the residue — `damp_wiggles` spreading condensate into forbidden cells
  (167 504 cells, 3 970 kg/kg) — plus the retained supersaturation at 778 g/kg and the fact that
  79 % of columns now radiate from the prescribed lid. Original entry: (item 74). `SaturationAdjustment` condenses into ~65 000 cells that
  `IceSchemeCommon::canCondense` rules impossible, and `damp_wiggles` — a numerical smoother —
  spreads it into ~132 000 more; `evaporateWhereImpossible` then deletes 257 000 cells' worth,
  correctly, because the cells are superheated (`p_sat > p`) or supercritical. Water is
  conserved to -0.0003 %. **The fix is to make the two agree** — the `canCondense` test at
  `SaturationAdjustment`'s entry, or keeping the smoother off condensate in such cells — and
  which is smaller is not yet measured.
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
