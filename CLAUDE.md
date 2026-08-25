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
| `MultiLayerRadiation.h` | Grey optical depth from column mass with pressure broadening; surface energy balance; the **closed-form** two-stream equilibrium solve (`ATM_RAD_DIRECT`, item 30) beside the inherited Lambda iteration (`ATM_N_LAMBDA`, default 4 = under-converged by 4× at 250 bar) |
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

3. **Radiation runs in mode 2** (direct σT⁴), **and `ATM_RAD_DIRECT` is the default solver since 2026-08-22** (item 71): the Lambda iteration converges monotonically onto the closed form (243.43 → 227.55 → 212.54 → **211.57** W/m² at n_lambda 4/64/512/direct), so the shipped 4 sweeps were **15 % high** at no saving in time. `ATM_RAD_DIRECT=0` restores it; **every OLR in this file predates the flip**. Item 30's "0.15 %" does not reproduce (13.09 % shipped, 8.69 % with `ATM_SKIN_GREY=0`) and is unattributed. Modes 0/1/3/4/5 all lean on the Scotese
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

**THE ENERGY BALANCE NO LONGER CLOSES, AND THAT IS THE REPAIR, NOT A REGRESSION** (item 67,
default since 2026-08-21). This file used to say the balance closed to -1.26 W/m2 and that this
was the least trustworthy number in it, because closing was what the `t_skin` fixed point
guaranteed. It was worse than untrustworthy: it was an identity created by a missing factor of
2. `t_skin` solved `sigma*T^4 = F`, the planet's EFFECTIVE temperature, where the grey skin
value `F/2` belongs — so the lid was assigned the temperature whose emission equals the entire
energy input, and the OLR duly descended onto it. With the factor in, the model absorbs 272
W/m2 and emits **136** at 200 iterations, imbalance **+136.23**. That number measures the
prescribed adiabat; the old one restated the prescription.

**AND THE CONVERGED VALUE IS `F/2` EXACTLY, SO THE IMBALANCE IS EXACTLY HALF THE ENERGY INPUT**
(item 73, the first 200-iteration series on the current defaults — item 67's pair predates both
`ATM_RAD_DIRECT` and `ATM_PROJ_SWEEPS = 10`). The OLR reaches `sigma*t_skin^4` at iteration 120
and stays within **0.01 W/m2** of it; item 67's "+123.76 and still widening" was the
under-converged Lambda solver, and **is superseded, not merely extended**. `sigma*t_skin^4` IS
`F/2` by construction, so the converged state is `OLR = F/2`, `imbalance = F/2` — 136.23 W/m2,
constant. The corrected closure is as arithmetic as the old one, in the opposite direction:
a prescribed isothermal lid emits `sigma*t_skin^4`, the column delivers exactly that, and
whatever `t_skin` is set to is what the OLR becomes.

**EVERY OLR AND IMBALANCE FIGURE ELSEWHERE IN THIS FILE AND IN THE README WAS MEASURED WITH
`ATM_SKIN_GREY` OFF**, i.e. on the old default. They are not retracted — they are correct
measurements of the old branch, reachable now with `ATM_SKIN_GREY=0` — but none of them
describes the shipped model any more. The series in item 73 is the one that does.

The model reproduces its design targets exactly. None of this makes the outputs predictions.
These are inputs, in rough order of how much they move the answer:

| Parameter | Value | Status |
|---|---|---|
| `kappa_H2O` / `kappa_CO2` / `kappa_bg` | 0.01 / 0.001 / 1e-6 m²/kg | Factor-of-2 uncertain, and **not a lever on the converged OLR at all**: 64× moves it 0.10 % (item 29), confirmed with a converged solver — 0.14 % over 16× (item 30) |
| `geothermal_flux` | 150 W/m² | Open. The ≥195 W/m² argument is retracted, and item 25 makes it worse: it enters the `t_skin` fixed point, so it helps set the very flux it was being compared against |
| `t_surf_equator` / `t_surf_pole` | 1500 / 1450 K | **Prescribed, not solved** |
| `t_skin` | 254.0 K start, **relaxes to 221.2-221.4** (was 262.96) | Item 67, **default since 2026-08-21**. The fixed point used to solve σT⁴ = absorbed, which is T_eff, where the grey skin value T_eff/2^(1/4) belongs — the lid was assigned the planet's whole energy input as its emission, so the budget closed by construction. Now σT⁴ = F/2; `ATM_SKIN_GREY=0` restores the old branch. **The 254.0 start is still derived by the OLD formula** and the consistent value is 213.6 K; left alone deliberately so item 67's measured pair describes the shipped default exactly |
| insolation | 0.71 S₀ | Faint young Sun at 4.4 Ga |
| `omega` | 3.17e-4 (5.5 h day) | Earth–Moon angular-momentum inversion (item 40): 5.5 h **is** the Moon at 5.95 R_E, 4.35× modern. Robust for a good reason — the Moon from 3 to 10 R_E spans only 5.0–6.1 h, so the bracket needs no tidal chronology. **Two Hadean-specific torques are omitted and neither is bounded**: thermal atmospheric tides on ~250× Earth's air mass (on Venus they spin the planet the *other* way) and dissipation in a molten surface. And the high-angular-momentum impact scenarios (Ćuk & Stewart 2012, Canup 2012) break conservation outright, toward a *shorter* day. Nothing observational reaches 4.4 Ga |
| `cell_lat_scale` | 0.33, scaling the **Hadley edge only** | Config parameter since item 32; **default was Earth's 1.0 for everything measured before it**. Held–Hou puts the edge at 5.1° because Ro_T is **1/12.5 of** Earth's (0.0048 against 0.0598) — smaller, so the cells are narrower; read the other way the argument inverts. Item 36: scaling every anchor left cells 10/40/40° wide with the extratropics a bare ramp, so item 31's scan compared layouts that differed in more than width |
| `n_cells_hemisphere` | **5** (Earth is 3) | Default since item 38; **everything measured before it used 3**. Three unconverged arguments: 20° bands against the Rhines ~22°, no 40°-wide-cell artefact, and Ψ maxima still on their prescribed cores at 400 iterations where n = 3's migrate off. **The Rhines argument rests on `omega`** (item 40): at the nominal 5.5 h the estimate is 4.5 cells/hemisphere, selecting 4 or 5 equally, and only the short end of the 4–6 h range makes 5 clear. So n = 5 is inside the rotation uncertainty, not selected by it; the other two grounds are the stronger pair. Not a claim the model *sustains* five cells — nothing maintains an indirect cell here. **AND UNTIL ITEM 69 FOUR OF THE FIVE TURNED THE SAME WAY**: `centreAmp` fills every middle cell with a Ferrel copy and Earth's polar template already carries the Ferrel's sense, so the prescribed `v_trop` was −3.0/+4.0/+4.0/+4.0/+0.5 and the initial Ψ was single-signed over almost the whole hemisphere. `ATM_CELL_ALTERNATE` imposes the parity `edgeRadialCoeff` already assumed and is **default ON in all three trees since 2026-08-24** (`=0` restores the co-rotating stack exactly). So every cell-count argument above was made about a layout that was not a multi-cell circulation. **MEASURED HERE AND IN ATHAD_COND 2026-08-25** (item 81): it is a structural knob with NO scalar signature — Psi sign bands at 40 iterations go 6/6/6/6/6/6/2/2/2 -> 5/5/5/5/5/5/5/4/4 by height here, and **1/1/2/2 -> 4/4/5/5 in ATHAD_COND, whose off arm is a SINGLE overturning cell**, while the OLR (135.52), photosphere, albedo and max water vapour are identical to every printed digit. **`Psi_max` is a FALSE NULL for this knob** — bit-identical to five figures in COND, because the Hadley cell is `k = 0`, already direct; read `survival.py`'s band count and the RMS over latitude of `Psi(ground)` (-6.7 % in COND). `div(rho u)/rho` is 4-5 % WORSE with it on in both trees, matching PERID's 3.5 % |
| `cell_amp_mode` | 0 (off) | Item 33. Scaling latitudes without amplitudes multiplies the initial meridional shear by 1/s. The decay tracks the **zonal jet**, not the overturning |
| `mc_M_max` | 100 kg/(m²s) | Convective mass-flux ceiling, a parameter since 2026-08-20 (was a bare 3.0 in **two** places). Earth's 3.0 implies σw = 2.5 m/s; ×42.97 kg/m³ here gives 107. **Never binds in ATHAD** — largest \|M\| over a 40-iteration run is 0.654, 153× below it — unlike ATHAD_COND, where 3.0 was pinned through the whole convective column |
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
  absorbed is arithmetic, not a result. **ITEM 67 FINISHES THIS PARAGRAPH: the arithmetic is a
  missing factor of 2, the fixed point is now `sigma*T^4 = F/2` by default, and on that branch
  the imbalance is +136.23 W/m2 rather than -1.26. The migration into the skin is NOT the
  carrier — measured, the corrected arm's skin fraction stays at 1.1 % while its OLR tracks
  `sigma*t_skin^4` just as tightly.** The -1.26 belongs to `ATM_SKIN_GREY=0`. **ITEM 73 CLOSES
  IT: on the current defaults the OLR converges to `sigma*t_skin^4` = `F/2` to within 0.01
  W/m2 by iteration 120, so the imbalance is exactly half the absorbed flux and is constant.
  Both closures are arithmetic; only the sign of the leftover changed.**
  **THE KAPPA SCAN THIS PARAGRAPH USED TO DEMAND IS DONE, AND IT MUST NOT BE RE-RUN ON THE
  PRESCRIBED BRANCH** (item 29, four 200-iteration runs: 64x in `kappa_H2O` moves the converged
  OLR **0.10 %**; item 60 confirmed it a third time at 0.083 % by a different species). The
  answer cannot change with the current defaults, and the reason is structural rather than
  empirical: `t_skin` is a fixed point of absorbed SW + geothermal, in which **no kappa
  appears**, and item 73 measured the OLR sitting ON `sigma*t_skin^4` to 0.01 W/m2. Item 67
  changed what that fixed point equals, not what it depends on; item 71 shrank the only part
  kappa could ever have moved (the Lambda solver's 11.82 W/m2 of residual excess) to 0.01. So a
  kappa scan on the prescribed branch re-measures `t_skin`, which is exactly what the
  `t_skin`-blocks-everything bullet below forbids. **The one kappa measurement still worth its
  wall clock is on a branch with NO prescribed top** — `ATM_PROGNOSTIC_T=1 ATM_RAD_DIRECT=1`,
  item 66's configuration, where `skin%` is 0.0 at every diagnostic and the dry column converges
  to 223.8 W/m2. No kappa arm has ever been run there, and item 66's own run predates items
  67-73 so it is not reusable as the baseline. ~39 min per arm at 400 iterations, 24 threads.
- **Not grid-converged**: 519 W/m2 at 260 km against 581 at 300 km with `im` fixed at 61
  (both at 100 iterations, both pre-item-22 — the check has to be redone, not just extended).
  **Item 39 says why, and says the retry should not be a shell-depth or `im` scan.** At the
  committed `im = 41` the entire photosphere is ONE grid cell of `dtau = 55`: level 35 sits at
  `tau_above` = 10.8 and level 36 at 0.82, so the column goes opaque-to-transparent inside one
  layer, and the two-stream sweep is first order in `dtau`. `im = 61` gives 3 photosphere
  layers and `dtau` = 7.6 — poor, not broken. The lever is `zeta`, not `im`. **On the current
  defaults it is 30x worse than item 39 recorded** (item 73): `tau_above` goes **1626 -> 0.497**
  between level 36 (218.2 km) and level 37 (236.4 km), so the whole opaque-to-transparent
  transition is inside one cell. The conclusion is unaffected; the magnitude is.
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
- **MEASURED 2026-08-20 (item 61), and it is worse than "the last digit".** The last digit is
  where it *starts*; the moist-convection trigger set is where it ends up. Four null pairs at
  40 iterations with moist physics on, same physics in every arm:
  **same binary, same 24 threads, twice** — bit-identical through iteration 20, then run-to-run
  divergence from residuum #25 at 7e-7, reaching 1.5e-9 of scale in Ψ and nothing in the printed
  scalars (so item 18's "bit-identical at a fixed thread count" is right for ~20 iterations and
  not exactly true after that). **A rebuild that touches no physics code** — indistinguishable
  from that. **A rebuild that re-arranges arithmetic in a physics loop** (adding `+ phi_s(i)`
  with `phi_s ≡ 0` did) — 5e-8 at residuum #3, and by iteration 40 the OLR differs **0.11 %**
  and 795 grid points of water vapour by up to **64 g/kg**. **24 vs 23 threads** — 3.5e-7 at
  residuum #2, and by iteration 40 the OLR differs **0.51 %**, `max M_u` by 56 % of scale.
  **In every pair Ψ_max, albedo and the photosphere are unmoved** (Ψ to 1e-6 or better). So the
  dynamics are not the amplifier and the convective triggers are: a trigger flipping relocates a
  cell, which moves the water where the photosphere is. **Consequence: do not quote a
  40-iteration OLR difference below ~0.5 % as a property of the model**, and state the thread
  count with any OLR number at all. Items 51-53, 59 and 60 reported 0.02-0.11 % — controlled
  same-binary same-thread comparisons, so not retracted, but inside the envelope.

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
`ATSAT`, `ATURAN`, `ATNEPT` (giants), `ASTIM` (impacts), and `ATHAD_COND` (this model's
post-condensation epoch: 60 bar, a 513 K sea, water subcritical).

**ATHAD_COND carries this file's five trace-gas mole fractions since 2026-08-20** — CH₄,
NH₃, H₂, CO and SO₂ at 0.014 each, where they used to be zero there. Only the *fractions*
are shared: 1.4 % of ATHAD's 250 bar is 3.5 bar per gas and 1.4 % of COND's 60 bar is
0.84 bar, so the two forks do not carry the same trace-gas inventory and no conclusion
about one transfers to the other on amount. What does transfer is the background: with
the traces in, COND's `R_bg` is 318.10 against this model's 317.26, because in both the
traces are ~71 % of the background by mass.

C++ class, file and function names are kept **identical to `ATOM_Precipitation`** so fixes
cherry-pick in both directions; only the outer shell is renamed (`libathad.a`, `cli/had`,
`config_athad.xml`, `pyathad`). Preserve that.

**Twenty-eight defects found so far. The twenty-fourth is the first that is not inherited at all
— ATHAD wrote it** (item 67): `t_skin`'s fixed point solves `sigma*T^4 = F` where the grey skin
relation is `sigma*T^4 = F/2`, so the transparent lid is assigned the planet's entire energy
input as its emission and the energy budget closes by construction. It has the twenty-first's
shape rather than the Earth-constant shape — a wrong formula with the right one written in a
comment nearby, here in the very file the prescription overwrites — and it arrived through a
*correct* repair that removed a circular measured-OLR reference and took the `/2` out with it.
**The lesson: when you remove a bad term from an expression, check what else was in the
parentheses.**

**Twenty-three defects found in the inherited code so far, all latent on Earth and live here.**
The pattern is consistent and worth expecting: *Earth's numbers as bare literals inside
physics kernels, each with a comment justifying it by Earth's conditions.* Examples —
`dr = 0.025` silently tied to `im = 41`; a 333.15 K cap written back into the prognostic
temperature; `287.0` J/(kg·K) as the density gas constant; convective triggers as absolute
hPa; `p_stat` cubically extrapolated at the lid. The twentieth is not a physical constant but
an **iteration count**: `n_lambda = 4` radiation sweeps, adequate on a 1 bar column and 4×
wrong on a 250 bar one (item 30) — the same pattern in a place nobody thinks to look for it,
since a loop bound does not read like an Earth assumption. The twenty-second is subtler still and
is item 50's: the Boussinesq buoyancy divides its temperature anomaly by `t_0` = 273.15 K, the
**non-dimensionalisation constant**, where the physical reference temperature `T_ref` belongs. On a
288 K planet the two agree to 5 % and the defect is invisible; here it overstates the body force
5.49× at the surface and 0.96× at the top, so it is a height-dependent distortion rather than a
rescaling. The twenty-third (item 52) has a different hiding place rather than a different
shape: `MoistConvection`'s `s` recurrences divide a whole mass-flux bracket by `s_0` when only
one term in it needs the division, and the resulting garbage was **absorbed by a safety cap** —
`MC_t` sat pinned at `MCt_max` while updraft parcels ran at 46 000 K, so every plot of the field
looked bounded and physical. When something behaves oddly, look for a constant that was true at
1 bar and 288 K — including the ones that look like units, and **look at what the caps are
holding back**.

The three worst are worth stating because they show the pattern's worst form — an
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

**The twenty-first is the newest, and it is not a literal at all — it is a wrong formula
with a confident comment next to it** (item 39). `exp_rm = 1/(rm+1)` is documented in two
places as the Jacobian of the radial coordinate transformation (`TurbulenceAtm.h`;
`PressureSolverAtm.h` writes it as `dp/dr_physical = exp_rm * dp/d(rad.z)`). It is the
Jacobian of a **quadratic** stretch, `z ~ (rm+1)²/2`, while `init_layer_heights` builds an
**exponential** one. Every radial derivative in the core is therefore mis-scaled by a factor
varying **11.8×** across the column — the core's radial length unit is 25 km at the surface
and 296 km at the top. The test is unit-free (`ratio(i) = [inv_2dr·exp_rm]/[1/(z[i+1]−z[i−1])]`
must be constant in `i`), so it cannot be argued away as a units convention.
`checkRadialMetric()` prints the spread at every startup; `ATM_METRIC_CHECK=1` adds the
per-level table. **The lesson is narrower and worse than "look for Earth constants": the
comment asserting the invariant was right there, in two files, agreeing with itself, and
agreeing with the variable's name.**

**Fixes worth porting back upstream** (not yet applied to ATOM_Precipitation as of
2026-08-11): the `t.x[-1]` out-of-bounds in `MoistConvection::findCloudBaseLFS`; the
`m_node_weights` OpenMP race in `GetMean_2D/3D`; the UB in `get_temperatures_from_curve`;
and `-MMD -MP` header dependencies in the Makefile.

**The `exp_rm` metric defect is live in ATOM_Precipitation and ATHAD_COND, and is worse
upstream** (item 39) — checked, not assumed. ATOM_Precipitation has the identical
`exp_rm = 1/(rm+1)`, the identical exponential `init_layer_heights`, and `zeta = 3.715`
hard-coded at `cAtmosphereModel.h:273`, giving a **23.2×** spread against ATHAD's 11.8×;
ATHAD_COND is at `zeta = 3.0`, so ~12×. ATJUP/ATSAT/ATURAN/ATNEPT have no `exp_rm` at all.
`checkRadialMetric()` is the diagnostic to port first — it is print-only and self-silencing
once the spread drops below 1.05, so it costs one startup line and cannot change a result.

**The constant-density meridional streamfunction (`MinMax_Atm.cpp:175`,
`const double rho = r_air`) is live in ATOM_Precipitation and ATHAD_COND** — checked, not
assumed; ATJUP/ATSAT/ATURAN/ATNEPT/ASTIM have no such diagnostic at all. It labels a volume
flux as kg/s. **ATHAD_COND is the urgent one** (same 250 bar class, 15 km scale height over
a 120 km shell), where it can hide a cell outright, as it did here. On Earth's 16 km shell
ρ spans ~6× rather than 10⁴, so it distorts rather than inverts — but ATOM_Precipitation
uses that diagnostic to judge exactly this question (`24ff23a` "revive Hadley/Ferrel
cells", `1e59daa` jet spin-down), and an instrument used to measure cell strength should
not be ~2× overweight at the cell core. Fixed here in `dabbc94`.

Traps already solved elsewhere in the family — check before re-deriving:
Coriolis/centrifugal signs (ATURAN `8b284cb`, `4201957`; ATNEPT `024c37f`, `e412b1b` —
ATHAD's dynamics already agree, its *diagnostics* did not); mass- not mole-weighted mixture
properties (ATNEPT `c116d71`); in-place Gauss–Seidel as a threading defect (ATURAN
`ffd0e0e`); report failures and limits in the README (ATURAN `74b4ded`, ATNEPT `34286b8`).

## Open risks

- **THE RADIAL METRIC IS REPAIRED, `ATM_METRIC_EXACT=1`, AND IT IS DEFAULT-OFF BECAUSE THE
  CORRECT JACOBIAN MAKES THE PROJECTION WORSE** (item 80). `metricExpRm()` replaces the 11 sites
  that each recomputed `exp_rm = 1/(rm+1)`; the true value is `metricShellLength()/J` with
  `J = zeta*L_atm*exp(zeta*(r-r0))`, running **6.36 at the surface to 0.317 at the top** against
  0.5 -> 0.333. `checkRadialMetric()` — unit-free, written in item 39 before the fix existed —
  goes from a **11.8x spread to 1.00x**. **A second defect was found repairing the first**: the
  Laplacian omits the `-(J'/J)*f'` curvature term in BOTH metrics (`J'/J` = `zeta` = 3 exact,
  `1/(rm+1)` legacy), added at 11 RHS sites and folded into the Poisson operator's existing
  anelastic off-diagonal; `metricCurv()` returns 0 on the legacy branch so that branch is
  unchanged. **Measured at 40 iterations: the OLR does not move AT ALL (135.52 both) — items
  29/41/60 confirmed from a fifth direction — while `div(rho u)/rho` rms goes 2.739e-02 ->
  7.722e-02 and `Psi_max` jumps from 52 153 at 36 km to 158 088 AT THE GROUND**, item 68's
  non-closure signature. This is item 72's zeta-scan direction reproduced **unconfounded** (same
  grid, one formula). Partly solver stiffness — `ATM_PRESS_SWEEPS=64` brings it to 5.258e-02,
  where the legacy branch does not respond to 64x at all — but ~2x worse survives. **The
  residual and item 72's Rhie-Chow lead are now the same question.** Off-branch null: every
  printed physics figure identical at 40 iterations, `residuum_atm` last digit only.
- **THE FREE-RUNNING MOIST COLUMN ESCAPES THE `t_skin` PIN AND LANDS AT HALF THE DRY COLUMN'S
  FLUX** (item 79). 400 iterations, `ATM_PROGNOSTIC_T=1 ATM_RAD_DIRECT=1`, moist from iteration 0:
  OLR **102.99 W/m2** against `sigma*t_skin^4` = 135.5, with `skin%` **0.0 at all twenty
  diagnostics**, and a photosphere that moves **half a kilometre and two kelvin** in 400
  iterations (282.6 -> 283.1 km, 260.6 -> 260.1 K). That is the first moist emission level in
  this model held by the column's own opacity rather than by the lid — against the prescribed
  arm, where 79 % of columns radiate from the lid. Half of item 66's dry 223.8 W/m2, because
  items 75-77 made the condensate real. Imbalance **+168 W/m2**: absorbing 271, emitting 103.
  **NOT CONVERGED** — the OLR bottoms at 90.58 by iteration 100 and rises monotonically to 102.99,
  still +0.47 over the last 20, and this file has been caught extrapolating that shape twice.
  Stable otherwise: water conserved to +0.0007 %, surface floats to 1491.6 K. **Open and rising:
  `max q_v` = 794.7 g/kg, still climbing, where the prescribed arm plateaued at 778.15.**
- **`restart_stride = 0` AND `checkpoint_save_iter = -1` STILL WRITE RESTART DUMPS** (item 79) —
  four 593 MB `.bin` files with both knobs off, 2.3 GB per run. Item 45 recorded half of this as
  an aside; it is live. And the output file indices key off `total_iter_count`, so a
  400-iteration run writes `panorama_900.vts`. A full-output run is 15 GB, of which 8.3 GB is the
  panorama (`paraview_panorama_vts_flag`, independent of `checkpoint`).
- **THE PRESCRIBED PROFILE DISCARDS THE LATENT HEATING, NOT JUST THE RADIATION'S ANSWER**
  (item 78). Invariant 3 says `densities()` overwrites what the radiation computed. Measured, it
  also overwrites what the LATENT HEAT RELEASE computed: the same cell at 256 km is **221 K in
  every diagnostic and 280.51 K inside `ThreeCatIceScheme`**, because the order per iteration is
  radiation (line 1446) -> dynamics -> moist block (line 1685, condensation releases ~60 K) ->
  `densities()` (line 2125, overwrites `t` with `max(t_skin, T_ad)`) -> diagnostics. **So a
  250 bar water-vapour atmosphere condensing 40 g/kg at its top cannot raise its own temperature
  by one kelvin between iterations, by construction**, and the microphysics runs on a field
  nothing else ever reads. It also closes items 74-75's loop: the reset to 221 K is what makes
  `SaturationAdjustment` condense hard on the next pass. `ATM_PROGNOSTIC_T=1` stops it; **what
  that branch does with the moist physics on, now that items 75-77 made the condensate real, is
  the run this line of work has been walking toward and is NOT measured.**
- **`moist_phys_start_iter` IS 0 SINCE 2026-08-23, AND EVERY DRY-COLUMN NUMBER IN THIS FILE
  BELONGS TO THE OLD 300** (item 74). The gate was ATOM_Precipitation's, inherited on the fork
  commit, justified by "a single tropical maritime column" — an Earth surface classification in
  a model with no sea — and waiting for a circulation item 18 puts 1e4 iterations away. At 200
  iterations there is no runaway. **But the moist column does NOT converge where the dry one
  did**: OLR 168.01 and still falling, 32 W/m2 above `sigma*t_skin^4`, against the dry arm's
  arrival at 136.19 by iteration 120. Set `<moist_phys_start_iter>300</...>` to reproduce the
  dry series, item 73's identity included.
- **ZERO PRECIPITATION IS A RESULT, NOT AN ABSENCE** (item 76, mechanism CORRECTED by item 77).
  **The `P_max_flux` cap binds on the FIRST ice-scheme call and never again** — 52 055 cells at
  it in call 1, wanting 2.83x, then 0 cells in every later call — so rain is not "forming at its
  cap and evaporating", it is **not produced at all** after `initCloudIce`'s deck is consumed.
  The cap is therefore NOT item 52's shape and hides nothing. The dry surface is still computed
  rather than asserted, and the condensate that exists still reads `0.000000` at every
  diagnostic, because `evaporateWhereImpossible` converts the falling
  flux to vapour at `P/(v*rho)` in the supercritical column below ~177 km. Water conserved to
  -0.0004 % over 200 iterations with no ceiling deletions. **ATHAD's dry surface is computed,
  not asserted** — that is the physical line between this epoch and ATHAD_COND's. The Earth
  phase bands were repaired alongside (`ATM_PRECIP_BANDS`, default on: the inherited flux always
  passes, only production is gated, and snow loses its -20 C floor) and the repair is
  **bit-identical at 40 iterations** — `prod_s > 0` in 0 cells, because every cell where the flux
  loop evaluates production is already above 273.15 K. Third time a suspected blocker was not the
  blocker (items 63, 74, 76). **The "rain source 50 K below freezing" is WITHDRAWN** (item 78):
  decomposed, `S_r` is entirely `S_c_au` at **280.51 K**, and the 221 K is what `densities()`
  writes at the END of the iteration — the ice scheme reads the cell mid-iteration, 60 K warmer.
  Still open is `P_rain_0 =
  max(P_rain.x[0], 1e-6)`, which normalises the flux by a surface value structurally zero here
  (`ATM_PRECIP_DIMENSIONAL`, default off).
- **FIXED IN ITEM 75, AND IT IS THE LARGEST EFFECT IN THIS FILE: `SaturationAdjustment`'s
  `clampAndFade` tested whether a cell could hold a condensed phase, condensed, added the latent
  heat that makes the answer FALSE, and never re-tested.** The admissibility test is applied to
  the temperature the cell ARRIVES with; the supersaturation removal forty lines below raises
  that temperature, `E_sat` with it, and past `E_sat > p` the cell is superheated and the cloud
  just written into it cannot exist. `ATM_SAT_SUPERHEAT` (**default ON**, `=0` restores)
  re-tests after the write-back and rejects the step. Measured at 40 iterations, 24 threads,
  moist on: **OLR 267.84 -> 135.52 (-49.4 %)**, photosphere 237.3 km/363.4 K ->
  **281.3 km/225.7 K**, max cloud water **0.000000 -> 40.09 g/kg**, cloud cells surviving the
  ice scheme 6 486 -> **149 271**. Against item 67's -27 % and every microphysics repair's
  0.02-0.11 %. **Three consequences, none comfortable**: the model is MORE `t_skin`-dominated
  (item 73's `OLR = F/2` now arrives by iteration 40, and **79 % of columns radiate from the
  prescribed lid** against 1.1 %); the retained supersaturation reaches **778 g/kg** against a
  design 672 and is untested past 40 iterations; and `damp_wiggles` still spreads condensate
  into forbidden cells (167 504 cells, 3 970 kg/kg) for the ice scheme to evaporate — not
  addressed. **Every moist number recorded before 2026-08-23 predates this.**
- **CONDENSATE IS CREATED AND DESTROYED EVERY ITERATION, AND THE DIAGNOSTIC SAMPLES ON THE
  WRONG SIDE OF IT** (item 74, mechanism CORRECTED by item 75 — `damp_wiggles` creates none of
  it, the mass is identical either side of the smoother, and the two predicates agree). Per iteration: 66 785 cells with cloud entering the moist block,
  131 404 after `SaturationAdjustment` condenses, **263 530 after `damp_wiggles`** — a numerical
  smoother spreading condensate into cells that cannot hold it — and **6 486** after the ice
  scheme. The remover is `IceSchemeCommon::evaporateWhereImpossible` and it is CORRECT: the
  cells are superheated (`p_sat > p`) or supercritical, where no droplet can exist, so there is
  no rate and no `dt` — invariant 2 enforced, water conserved to -0.0003 %. **The defect is the
  disagreement upstream**: `SaturationAdjustment` and `damp_wiggles` both manufacture a phase
  `canCondense` forbids. Not fixed. And note what it did to the instrument — `max cloud water`
  prints **0.000000** because the diagnostic sits after the whole block, while the radiation
  sees real condensate mid-iteration and returns albedo 0.4964. **Before explaining a field
  that will not move, find out where in the iteration it is read.**
- **THE PRESCRIBED ADIABAT AND THE GREY OPACITY ARE NOT COMPATIBLE DESCRIPTIONS OF THE SAME
  ATMOSPHERE, AND THEY DISAGREE BY THREE ORDERS OF MAGNITUDE IN FLUX** (item 73). Item 67's
  named follow-up — replace the isothermal clamp with `max(T_rad(tau_above), T_ad)`, the grey
  radiative-equilibrium profile `sigma*T^4 = (F/2)(1 + 3*tau/2)` — is **written
  (`ATM_SKIN_TAU`, default 0), measured and WITHDRAWN**. Measured on the model's own
  `tau_above` array: `T_rad(tau)` exceeds the prescribed profile at **all 41 levels** (6.75x at
  the surface, 2.5-5.6x through the middle, 1.00x only at the lid), so there is **no
  radiative-convective crossing to switch at**. `tau_s` = 2.26e6 in radiative equilibrium at
  `F/2` = 136 W/m2 wants a **9496 K** surface against the model's 1408 K; read the other way,
  the OLR a 1408 K surface can drive through that opacity is **0.13 W/m2**. The two failure
  modes are one fact at two altitudes: whole-column (`=2`) runs away to 15 599 K in 4
  iterations, skin-region-only (`=1`) ignites a **polar** runaway at 236 km after ~150
  iterations — 2729 K, 18 of 181 latitudes hotter than their own surface. One of the three —
  the 1500 K surface, `kappa`, or the greyness — has to give, and not by a one-line change to
  the top boundary.
- **THE MERIDIONAL STREAMFUNCTION DOES NOT CLOSE AT THE GROUND, AND THE VIOLATION IS 2.1× THE
  CIRCULATION IT HIDES** (item 68). `Psi(ground)` must be 0 — `u` at `i = 0` is identically
  0.0000 m/s at every latitude and `|Psi|` at the lid is ~0, so the column-integrated meridional
  mass flux is forced to vanish. Measured: **9.891e13 kg/s against a 4.733e13 interior maximum**.
  **98.7 % of it is present at iteration 1**, so it is written by the initial velocity profile,
  not by the dynamics. Three nulls: the time-loop solver (`ATM_PRESS_SWEEPS` 64× → +0.3 %, wrong
  direction), **`exp_rm` (a `zeta` scan to the exactly-correct metric at `ln(1.5)` → no trend,
  and that arm is the worst of four)**, and the cell compression (ratio invariant at 2.08–2.40
  while absolute Ψ moves 4.5×). The lever is `ATM_PROJ_SWEEPS`, the **initial** projection,
  default **1**: 10 sweeps cut the RMS **52.5 %**, 100 sweeps 55.9 % — **it plateaus, so ~44 % is
  structural and unexplained**. Grid-scale noise is excluded (2Δ oscillation index 0.006–0.2, where a checkerboard is ~4). **ITEM 72 SHOWS THE RESIDUAL IS STRUCTURAL AND DOES NOT ATTRIBUTE IT**: `div(rho u)/rho` rms **1.611e-02 is bit-identical under 64× the sweeps**, so the projection has CONVERGED to a fixed point that is not divergence-free — the discrete div and grad are not adjoint on this collocated stencil. Five mechanisms excluded (base-state density 1.0 %; `exp_rm` **worse when made exact** — a second independent exoneration; no Nyquist mode in velocity, `Psi` or `p_dyn`; not a boundary artefact — it is distributed, 11 % in the bottom 7 km and 33 % at 80–160 km). The lead is **Rhie-Chow face reconstruction**, named as un-done in `PressureSolverAtm.h`, but the checkerboard that story predicts is absent, so it is a lead and not a diagnosis. **The model has been printing this residual in every run log all along.** **Method: `Psi(ground)` must be read as an RMS over latitude, not
  a max** — the max sits inside the Hadley cell, so a change localised elsewhere reads as
  bit-identical to nine figures while the field moves 55 %. **`ATM_PROJ_SWEEPS` DEFAULTS TO
  10 SINCE 2026-08-22** (`=1` restores the old branch): RMS Ψ(ground) −53.3 %, the real
  circulation above 20 km **+14.3 %**, and `Psi_max` −50.6 % — because at one sweep the global
  maximum of Ψ *was the spurious surface flux*, so **`Psi_max` was reporting the defect, not the
  circulation**. The two maxima coincide now and did not before. **Every Ψ figure in this file
  predates it.** The OLR is unmoved (243.43 W/m² both arms).
- **EVERY ZONAL (LATITUDE–HEIGHT) FIGURE MADE BEFORE 2026-08-22 WAS DRAWN ON A DISTORTED AXIS**
  (item 70). `paraview_vtk_zonal` wrote the vertical coordinate as **level index**, which on this
  exponentially stretched grid stretches the bottom and squashes the top by **18.6×**, varying
  with altitude so no ParaView aspect setting could undo it. Contours, glyph angles and
  streamline curvature all inherited it. Two smaller defects went with it: the glyph vector was
  scaled `1/u_0` where the scalars beside it used `u_0` (**64× = u_0², direction unaffected**),
  and it was raw m/s on an index-space geometry. Fixed: the axis is true height (plot-x/height
  constant to 1.000000), and `uv_plot` carries the field in plot units per day. Verified against
  physics — velocity and `dPsi/dz` agree in sign at **100 %** of sampled points, magnitudes to
  0.2 %. **No physics reads these fields**, so no computed result changes; but any conclusion
  drawn by eye from an older zonal plot should be re-examined.

- **Drag is eliminated as the cell-decay driver, and the argument that eliminates it removes
  more than drag** (item 49). The four-arm scan returned `Psi_max` agreeing to **eight
  significant figures** across 100× in `rayleigh_kf` — which was the **reproducibility noise
  floor, not a result**: `surf_drag` carried its own `*dt` and RK4 multiplies every RHS by `dt`
  again, so the drag entered as **dt²**, 4.5e-08 over a 200-iteration run. Item 34's pair, and
  the Held–Suarez block thirty lines above had made the identical repair in 2026-07. Fixed in
  `3e2d78f` (ATHAD_COND `a609d2b`); **every Ψ, KE and wind number in these files predates it.**
  Correctly scaled, drag explains **0.034 %** of the Ψ decay, 0.31 % at ten times the rate,
  linear in `kf` to 9.25 against 10. **The reason is a timescale**: `rayleigh_kf` is 1/86400 s⁻¹
  and 200 iterations is 39 s – 12.5 min of physical time (item 47), so **nothing slower than
  minutes can explain a decay that completes in 200 iterations** — drag, radiative relaxation
  and surface exchange all go together, leaving the initialisation transient and geostrophic
  adjustment. **Two method notes worth more than the result.** A null at the noise floor is not
  a measurement, and the scan had to be diagnosed before it could be read. And the term was
  checked to be *connected* — `kf` = 1.0 moved Ψ, `surf_drag` has no enclosing conditional —
  because the same day a claim in ATHAD_COND was retracted for the opposite error, liveness
  inferred from a grep instead of read from the control flow.
- **The photosphere is measured now, and it is 21 km higher than this file used to say**
  (item 42). `tau_above` (cumulative LW optical depth from the lid) is a real array; before it
  the photosphere was computed **nowhere in the model**, and every figure quoted for it came
  from an offline dry, cloud-free Python column. Measured with moist physics on it is
  **239.4 km / 325.2 K**, against the reconstruction's 218 km / 402.6 K — clouds add
  `k_liq·LWP + k_ice·IWP` and push the crossing up. **Re-read item 39's "one grid cell of
  dtau = 55" with that 21 km in mind.** `printPlanetaryBalance` now prints the height, the
  temperature there, and the fraction of columns radiating from within 1 K of `t_skin`. New
  companions: `tau_layer` (per-layer dtau, the resolution measure), `ubud_*` (the **radial**
  momentum budget, absent while θ and φ both had one — item 28's spurious 293-rms radial
  acceleration was invisible for exactly this reason), `brunt_N2` (measures "neutrally stratified by
  construction" instead of asserting it, and tests invariant 4), and `Psi` as a field. In
  ParaView and Results — except `ubud_*`, whose six VTK dumps per slice were dropped on
  2026-08-20 (a third of the file for a diagnostic that is read as a min/max, not as a field).
  The arrays are still computed and still printed by `Results_Atm`; only the dumps are
  commented out, in one block per writer, so restoring them is a four-line edit.
  **Field names changed in item 61**: the Brunt-Vaisala frequency squared is `brunt_N2` in the code
  and `BruntVaisala_N2` in the VTK files — it used to be `N2`, in the same file as nitrogen's
  `q_N2`. The eight species are written under bare names (`H2O CO2 N2 CH4 NH3 H2 CO SO2`), all
  mass fractions from `AtmMixture::split()`, so what is plotted is what is integrated; the raw
  transported CO2 array is `CO2_tracer` and answers a different question (what the transport
  advances, not what the physics reads — they differ by item 57's dilution, 0.1359-0.2805 against
  a uniform 0.2053). **`initCloudIce`'s H_crit is on a fraction of surface pressure now**,
  which is right but measured as a 0.18 % effect on OLR, not the lever it was billed as: the
  albedo cannot respond to cloud *amount* at all, because reflectivity saturates on presence.
  **Item 53 makes it five for five** — every repair in it (three `s` scalings, `MC_t`'s `t_0`, the lid
  and surface BCs, the bounded recurrence, `cc_factor`) moves the OLR by 0.03 % in total and leaves Ψ
  identical to eight digits.
  **Item 51 makes that three for three**: a local mixture `cp` in `MoistConvection` (the constant `cp_l` = 2040 is the
  mixture cp at ~1101 K, so +30.5 % wrong in the skin where all the condensation happens) moves rain and snow production
  **−18 %** and the OLR **+0.02 %**, with Ψ and the photosphere bit-identical. ~~**`albedo_cloud` is the only path condensate
  has to the radiation**, so expect microphysics work to be unmeasurable until that parameterisation responds to something.~~
  **REFUTED BY ITEM 75 (2026-08-23): condensate has a SECOND path and it is the larger one.** The
  long-wave cloud opacity `k_liq*LWP + k_ice*IWP` in `MultiLayerRadiation` responds to condensate
  AMOUNT, not merely to its presence. It looked inert only because the condensate was being
  annihilated every iteration (item 74); with `ATM_SAT_SUPERHEAT` keeping it alive, the OLR moves
  **−49.4 %** and the photosphere climbs 44 km. **The albedo half of the claim still stands** —
  reflectivity does saturate on presence, `albedo_mean` is 0.4964 against `albedo_cloud` = 0.50 —
  but "expect microphysics work to be unmeasurable" is wrong, and it is the reason five
  consecutive repairs were read as confirming a wall that was not there. **Do not use this
  paragraph to decide that a microphysics change cannot matter.**
  **Item 52 is the fourth**: repairing the `s` scaling in `MoistConvection` cuts `MC_t` by 2.9–3.5× and snow production
  by 5.9 %, and moves the OLR by −0.02 % with Ψ, albedo and photosphere identical.
- **`s`, `s_u`, `s_d` are a normalised temperature, not entropy and not dry static energy** (item 52,
  audited and measured). They hold `cp_l·T/s_0`; `param.py` calls `s_0` "entropy at 0 °C" and
  `cAtmosphereModel.h` called the arrays "dry static energy" — DSE is `cp·T + g·z` and **the `g·z` is
  absent**, which on a 300 km shell is the *larger* term above ~157 km (2.7e6 against 4.2e5 J/kg at the
  lid) and is worth 5 K per layer at the surface rising to 140 K at the top. So the updraft parcel never
  cools as it rises, while `findCloudBaseLFS`'s θ_e parcel — which sets the trigger and CAPE — does.
  Three scaling errors sit on top: the updraft recurrence divides by `s_0` when nothing needs it, the
  downdraft divides the whole bracket when only its latent term does, and `MC_t` converts back to kelvin
  with `t_0` where `s_0/cp_l` = 134.57 K is meant. **`ATM_MC_S_CONSISTENT=1` repairs all three, default
  off, off-branch bit-identical.** Measured at 40 iterations: shipped runs the spin-up with `MC_t` **at
  its cap** and parcels at **46 000 K**; repaired, `MC_t` is 120× below the cap at iteration 10 and 2.9–3.5×
  smaller at 20–40, snow −5.9 %, **OLR −0.02 % and Ψ/albedo/photosphere identical**. Still open: the
  missing `g·z`, `MC_t`'s second term carrying an extra `t_0`, `q_v_u` reaching 3.5 kg/kg, and `bcRadius`'s
  cubic lid extrapolation giving `s` **negative** values at 300 km.
  **ALL OF THOSE ARE REPAIRED AS OF ITEM 53 except the `g·z`**, and the `s` repair is now the default
  (`ATM_MC_S_LEGACY=1` restores the old behaviour). A fifth defect turned up during the repair and was
  the one that actually produced the impossible `q_v_u`: **`cc_factor`'s reference temperature was
  Earth's 288.15 K**, so the moisture-seed cap `q_v_u_add·q_sat(1500 K)/q_sat(288 K)` came out at
  **3.47 kg/kg** — a mass fraction of 3.5. Repairing the recurrence did not move it at all; referencing
  the ratio to the model's own surface temperature did (3564 → 725 g/kg). All of it is ported to
  ATHAD_COND, where `s_0`/`cp_l` makes the `MC_t` error 1.342× rather than 2.03×, and where the surface
  half of the BC defect was measured (`s_u` = −1020 K at 0 m).
  **THE `g·z` IS WRITTEN NOW (`ATM_MC_GEOPOTENTIAL=1`, default off) AND IT WAS NOT THE BLOCKER**
  (item 61). `s` becomes a real static energy — 1.567–11.108 → 10.229–14.268, nearly uniform, as a
  column on its own adiabat must be — but `c_u` stays **identically zero and the parcel temperature is
  302.7 K in both arms, to the digit**, because the updraft **never rises a single level**: cloud base
  is level 37 and the LFS level 38 in 99 of 173 convecting columns, and the recurrence
  `for(i = i_base+1; i <= i_LFS-1)` is an empty loop when they are adjacent. `|M_u| > coeff_recurr`
  holds in **one level and 7 of 181 columns**. So item 53's "the only thing between this scheme and an
  updraft that condenses" is **retracted: the blocker is geometry, not thermodynamics.** The one place
  the geopotential does act is the DOWNDRAFT, which descends from the LFS to the surface, and it moves
  `MC_t` by 5× — peak heating −81 %, peak cooling ×33.
  **ITEM 63 FIXED THE TRIGGER LEVELS AND ATHAD DID NOT MOVE — invariant 2 is what pins the base here,
  not the gate. The demonstration is in ATHAD_COND**, whose 20 km convective column turns `max c_u`
  from 0 to 7.19e-04 g/kg/s at 24.9 km with the knob on, nothing integrated moving. So the knob is
  right and ATHAD is simply not a place it can pay.
- **The radial momentum balance is the pressure gradient and nothing else** (item 42, measured): `ubud_pgf` 1.15 against `ubud_cor` **exactly 0** (non-traditional Coriolis is genuinely off, not just documented off), `ubud_advh` 0.0085, and `ubud_buoy` **1e-6** — item 34's extra `*dt` measured, the buoyancy body force is ~1e6 down and effectively absent. **`ATM_BUOY_CONSISTENT=1` corrects it and is default off** (item 50): the consistent coefficient is 2e7× the shipped one and 7.2× the 336 recorded as having driven a polar vertical runaway. It also repairs the Boussinesq reference temperature — the shipped term divides the anomaly by `t_0` = 273.15 K instead of `T_ref`, overstating the force 5.49× at the surface and 0.96× at the top. With (a) alone at 20 iterations, `ubud_buoy` goes 0 → −34.6 (~4× the pressure gradient), the radial wind rises 19.5 % and Ψ moves 0.0002 %; `buoyancy_ramp` was 0.067 there, so that is **not** a stability result. `brunt_N2` — the Brunt-Vaisala frequency squared, renamed from `N2` in item 61 so it cannot be read as nitrogen — confirms invariant 4: 0 through the column, 2.74e-4 s⁻² only in the isothermal skin — **but read item 55 before quoting it**: `brunt_N2` is formed at the end of `densities()` from the `t` and `p_stat` that same call just wrote, so it re-reads the prescribed profile and is identical to 0.03–0.5 % across runs whose `p_dyn` radial structure differs by **2500×**. Below the skin its non-zero part is a **second-order truncation error** of the adiabat integration (`brunt_N2/dz²` constant to ±30 % over a 12× range of `dz`), i.e. the grid stretch, not stratification. It measures the atmosphere only under `ATM_PROGNOSTIC_T=1`.
- **The imbalance is not a second number** (item 43). `t_skin` is **262.95 K at every one of ten
  diagnostics across 200 iterations** while the OLR falls 338.85 → 281.59, so σT_skin⁴ = 271.2 is
  constant and the reported imbalance is *identically* the OLR's distance from it (281.59 − 271.2
  = 10.4 against a printed −10.51). **"The energy balance closes" and "the OLR arrives at a number
  that never moved" are the same event.** Any small imbalance reports only how far the OLR has
  drifted toward σT_skin⁴. Item 25's mechanism is confirmed in direction — the skin-emission
  fraction rises 1.1 → 39.2 % as the imbalance closes and the photosphere cools 368 → 270 K — but
  **39.2 % is not a saturation and is not extrapolated**; both rates are decelerating.
- **The model cannot reach radiative equilibrium by integration, and that reframes invariant 3**
  (item 47). Column heat capacity `p/g*cp` = 5.20e9 J/m2/K, so the measured -1130 W/m2 net drains
  it at **6.85 K per YEAR** — 100 K takes 14.6 years. A 901-iteration run covers **under an hour**
  of physical time (177 s or 3375 s depending on which L the time unit uses — see below), i.e.
  **5-7 orders of magnitude short**. So item 46's "plateau" is a transient, not an endpoint, and
  its mechanism was wrong: **the surface is NOT prescribed during a run.**
  `MultiLayerRadiation.h:452` already floats it via a surface energy balance and it settles at
  1498.7 K, radiatively locked to a column that radiates 280 kW/m2 down onto it (suppression
  x236). The reservoir is the 250 bar of gas. **Invariant 3 treats the prescribed adiabat as a
  deficiency; it is in fact the only equilibrium obtainable here** — the alternatives are
  operator-splitting the thermal equation, cutting the column heat capacity during spin-up, or
  solving the RC column directly. A Neumann `p_stat` BC would NOT help: `p_s` is the weight of the
  air above, not a free value, and item 19 pinned it for exactly that reason.
- **No run in this file has a stated physical duration** (item 47). The time unit is "L/u_0" and
  which L is ambiguous — `L_atm` (15.7 km) against `metricShellLength` (300 km), a **19x**
  difference. Item 41's defect in the time coordinate. Every "N iterations" is a count, not a time.
- **There is no functioning surface boundary-condition layer** (item 48). All three sites --
  `BC_Atm.h:433`, `BC_Atm.h:1108`, `ThermoAtm.h:1540` -- are written as "copy from
  `i_topography`", which invariant 1 fixes at 0, so all are identities `x[0] = x[0]`. Surface
  values are whatever the last physics routine wrote (MLR for `t`, `densities()` for
  `p_stat`/`r_dry`/`r_humid`) — chosen by call order, not by design. **Any surface condition must
  be written, not enabled.** The dead-land-branch pattern at the scale of a whole layer.
- **Freeing the column does not free the model: the prescribed SURFACE is the next pin down**
  (item 46 — mechanism SUPERSEDED by item 47; the plateau is real, the explanation was not). Continued to 901 iterations, the prognostic OLR **stops falling** — a turning point
  near 660-680, then a plateau at **~1400 W/m2 (+-6 %, flat to slightly rising)**, with `T_ph`
  bottoming at 357.89 K and returning to 360.2. That is **5.2x the 271 absorbed**, a persistent
  -1130 W/m2 imbalance that *cannot* close: `t_surf_equator` = 1500 K is prescribed, so the
  surface is an infinite reservoir and the steady state is the flux it drives through the column,
  not a TOA balance. **A genuine radiative equilibrium needs the surface temperature SOLVED, not
  just the column freed.** Restart fidelity was attempted and is INCONCLUSIVE (+-5-6 % deltas
  against +-5 % intrinsic scatter, a one-iteration diagnostic phase shift, and a binary that
  differed between runs); `t_skin` matched exactly and the 29 serialized arrays are the full
  prognostic set, so evidence favours fidelity without verifying it. **Sequence runs after
  builds.**
- **The `t_skin` pinning belongs to the PRESCRIBED profile and is escapable** (item 45). With
  `ATM_PROGNOSTIC_T=1 ATM_RAD_DIRECT=1` there is no isothermal top for the radiating level to
  migrate into, and `skin%` is **0.0 at all twenty diagnostics** over 400 iterations — the
  mechanism item 43 measured cannot operate. **But what escapes it radiates 2927 W/m², 10.8x the
  271 absorbed, and is still falling at 400** (extends item 30: 6211 at 200 was not a floor, it
  fell a further 53 %; 2927 is not claimed as one either). `t_skin` is still 262.93 throughout,
  because it is a fixed point of the ABSORBED flux and independent of the column — so the
  imbalance is again exactly the OLR's distance from a constant, in a regime 10x different.
  **Open**: the effective emission temperature is 476 K against 368 K at `tau_above = 1`, so
  emission comes from deeper and hotter than the photosphere. Understand that before spending
  more integration. Also: `restart_stride = 0` did NOT disable the restart dump.
- **THE `t_skin` PIN HAS A CAUSE, IT IS A MISSING FACTOR OF 2, AND THE FIX IS THE DEFAULT AS OF
  2026-08-21** (item 67). `ATM_SKIN_GREY=0` restores the old branch; every OLR and imbalance
  number recorded before that date was measured on it. Both sites that
  set `t_skin` solve `sigma*t_skin^4 = absorbed SW + geothermal`, which is the planet's
  EFFECTIVE EMISSION temperature; the top layer of a grey atmosphere is at `T_eff/2^(1/4)`,
  because it sees no downward flux and re-emits half of what passes up through it. 263.07 K
  against 221.22 K. **All shortwave is deposited at the surface** (one site, MLR's surface
  balance), so this is the classical skin problem and `sigma*T_skin^4 = F/2` is exact, not a
  fit. **`MultiLayerRadiation` already contains the right relation twice** — the flux sweep's
  top boundary reduces to `sigma*T^4 = up/2` and its comment says so, and the direct solver's
  validated `eps -> 0` limit is `T_s/2^(1/4)` — so the solver and the prescription that
  overwrites it disagree by 2^(1/4) every iteration. **Consequence: a lid of emissivity 0.0068
  is prescribed to emit 271.10 W/m2 against an absorbed 271.10 W/m2, six figures, every
  diagnostic.** That is why the budget closes: item 25's -1.26 W/m2 and item 43's "the
  imbalance is the OLR's distance from a constant" are one identity, not two findings.
  `ATM_SKIN_GREY` applies the factor (**default on**; `=0` is a verified null against the
  pre-change binary). Measured at
  the 200-iteration pair: OLR 272.53 against **147.37**, imbalance -1.41 against **+123.76** —
  the **first FORCING that moves the converged prescribed OLR** (structural repairs have, e.g.
  item 22's initialisation order; no knob had). The shipped arm reproduces item 25 (-1.41
  against -1.26), so it is the same measurement.
  **BOTH ARMS DESCEND ONTO THEIR OWN `sigma*t_skin^4`** — 271.11 and 135.55, ratio 0.5000
  exactly — and this is measured, not extrapolated: the OLR's *excess* over that constant decays
  geometrically in both, ratio settling at 0.71 and 0.78, with the shipped arm five e-foldings
  further along the identical curve. **So the shipped model closes its budget for one reason
  only: `sigma*t_skin^4` was set equal to the absorbed flux.** No limit is claimed for the
  corrected arm — 11.82 W/m2 of excess remains at 200.
  **THE SKIN FRACTION KILLS ITEMS 25/43's MECHANISM.** Shipped, `skin%` climbs 1.1 -> 53.0 % and
  **saturates** (retiring item 43's "39.2 % is not a saturation"); corrected, it stays at
  **1.1 % through 180 iterations**, 6.6 % at 200 — **yet its OLR tracks `sigma*t_skin^4` just as
  tightly**. The emission level migrating into the skin is a CONSEQUENCE of a warm lid, not the
  cause of the pin: a lid at 262.96 K sits within 1 K of the profile the photosphere climbs
  through, one at 221.12 K does not. Items 25 and 43 read a real correlation backwards.
  **The corrected model is much further out of balance** — absorbing 271 and emitting 147 — which
  is invariant 3's warning again, that the better-looking number was the more assumed one.
- **`t_skin` blocks every radiative measurement, and that is the top priority** (item 41).
  `OLR = σT_lid⁴ = absorbed SW + geothermal`, exactly, and the OLR has now failed to respond to
  four separate 6×–500× forcings: κ (64×, 0.10 %), the circulation (500×, 0.03 %), `im`
  (analytically) and the grid stretch (6×, 0.36 % — and that driven by albedo, not resolution).
  Items 29, 39 and 41 are all unanswerable until this is broken. **Do not spend runs on opacity
  or grid questions first; they re-measure `t_skin`.** Corollary from item 41: the famous
  **−1.26 W/m² residual is a `zeta = 3.0` artefact** — at `zeta ≤ 1.5` the imbalance closes to
  0.00 exactly, so the shipped grid simply fails to reach its own fixed point in 200 iterations.
  And `initCloudIce`'s pressure-keyed deck (`p_crit = 1000` hPa, Earth's surface pressure) is
  the **only** live path from the grid to the OLR: max cloud water moves 37.3 → 49.9 g/kg
  between `zeta` 3.0 and 2.0, with moist physics off.
- **The dynamics and the radiation do not agree where the levels are** (README item 39).
  `exp_rm = 1/(rm+1)` is a quadratic-stretch Jacobian applied to an exponentially stretched
  grid, so the core's radial derivatives are mis-scaled by a factor varying **11.8×** from
  surface to top, while `get_layer_height()` — which the radiation, `ConvectiveAdjustment`
  and every diagnostic use — has the true heights. The core's effective bottom layer is
  ~15 km where the radiation's is 1.2 km. **It is not a function of `im`** (11.8× at 41,
  12.3× at 61), so do not reach for the level count; `zeta` is the lever and
  `zeta = ln(1.5) = 0.405` would make `exp_rm` correct by construction. Three independent
  arguments — photosphere resolution, the diffusive CFL, and this metric — all say cut `zeta`
  hard, and the CFL one says it is ~13× *cheaper*, not more expensive. **The open decision is
  which repair**: replacing `exp_rm` is principled but touches ~91 sites across 10 files
  including the Poisson operator and moves every number in the README, while cutting `zeta`
  makes the existing metric accidentally correct. Settle that before any further `im` or
  shell-depth work. `checkRadialMetric()` prints the spread every startup;
  `ATM_METRIC_CHECK=1` adds the table. What is measured is the *shape*; the absolute
  consequence (1.6× at the surface, 18.8× at the top against `L_atm`) is indicative only —
  the full non-dimensionalisation has not been traced.
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
- **The initial circulation is balanced against the model's own momentum equations —
  BOTH of them** (`initBalancedState`, README items 26-28, on by default since item 27,
  two-component since item 28). Without any balance the prescribed cells are buried within
  ~5 iterations by a linearly accelerating drift, and at 200 iterations that drift is twice
  the strength of the cell it replaced. Three things go with it. The balance must be written
  to **`p_dyn`**, because `p_stat` appears nowhere in the momentum equations — the entire
  meridional pressure-gradient force is `p_dyn`, and the 50 K equator-to-pole contrast exerts
  none of it directly. The radius in it is **`metricRadius(rm)`, not `rad.z[i]`**, worth a
  factor of ~21 and got wrong the first time. And **balancing one component is not a
  balance**: the θ-only version of item 27 created a radial acceleration of 293 rms where the
  model had none, because with the default switches `rhs_u` at u = v = 0 is `−dp_dyn/dr·exp_rm`
  and nothing else — no `nontrad` Coriolis, no curvature term, and `buoyancy_ramp` = 0 at
  iteration 0. It drove the vertical wind from 0.11 to 11 m/s over 200 iterations with the
  tropics sinking over the 1500 K equator, and Ψ_max never saw it. **Read the switches, not
  just the RHS**: `metric_curvature()` and `coriolis_nontraditional()` are both false by
  default, so terms that are written in `RHS_Atm_Turb.cpp` are not necessarily terms the model
  applies — the item-27 balance was built on one of them.
- **`p_dyn_ceiling` is 2000, not the inherited 10/3**, and is **untested beyond 200
  iterations**. The Earth values were keyed to steep orography this model does not have and
  forbade the balanced state outright (94.7 % of columns clipped). One accessor,
  `cAtmosphereModel::pDynCeiling()`, so the solver and the diagnostic cannot disagree.
- **The radiation solver was never converging, and now has a closed form** (item 30).
  `n_lambda = 4` is an Earth constant: the Lambda iteration is Jacobi on a 41-link chain and
  needs O(N²) sweeps, so 4 gives an OLR 4× too high on a prognostic column (65 982 against a
  converged 15 769). It does not need iterating — `a_i + b_i = 1` makes the net flux constant,
  which closes the system in two O(N) passes, exact to 8 figures against 10 000 sweeps and
  3200× faster. **`ATM_RAD_DIRECT` is default-off pending a longer measurement.** It changes
  the standard prescribed configuration by 0.15 %, so no quoted OLR is at risk — but with it
  on, the PROGNOSTIC column stops sitting on a hot branch and starts relaxing (OLR
  15 764 → 8 542 → 6 211 over 200 iterations, still falling). Where that lands is the open
  question, and **a monotone trend is not a limit** — this file has now been caught on that
  twice in one day, both times on radiation numbers.
- **The OLR does not respond to anything, now measured in both directions.** A **64×** change
  in `kappa_H2O` moves the converged OLR by **0.10 %** (item 29, four 200-iteration runs) and
  a 500× change in the circulation moves it 0.03 % (item 27). σT_lid⁴ is 271.11 W/m² in every
  one of those runs and `t_skin` is 262.96 in every one. Treat any OLR number as a statement
  about `t_skin` until that is broken. (The 0.9 %-over-8× figure this file used to carry was
  written before any scan existed — see item 29.)
- **THE PROGNOSTIC COLUMN CONVERGES, AND ITEM 45's 2927 W/m² IS RETRACTED** (item 66). 400
  iterations with `ATM_PROGNOSTIC_T=1 ATM_RAD_DIRECT=1`: the dry OLR falls 474 → **223.8 W/m²**
  with successive differences 111 → 0.65, a decelerating limit rather than a trend, against
  absorbed + geothermal of 271.57 — an imbalance of −47.7 W/m², the first prognostic number here
  of the same order as the energy input. Item 45's "10.8× absorbed and still falling" was the
  `n_lambda = 4` solver, as item 65 inferred. Its MECHANISM survives: skin emission is 0.0 %
  throughout, so the branch does escape the `t_skin` pin — now to a converged value. `t_skin`
  still never moves (262.99 → 263.01), but the OLR is no longer σT_skin⁴. **The moist onset at
  iteration 300 is a different story**: 223.8 → 1022, then 767, 751, 1036, 782 — ±35 % with no
  limit in 100 iterations, so the dry column is converged and the moist one is not.
- **A large part of the free-running top's cooling is the UNDER-CONVERGED radiation solver**
  (item 65). Same trace, same cell, only the solver changed: T at level 38 goes 266 → 219.7 →
  228.3 K with the default four Lambda sweeps and 266 → **228.9 → 243.3** with
  `ATM_RAD_DIRECT=1`, and the supersaturation ratio at call 3 falls **1934 → 97.9**. So item 30's
  `n_lambda = 4` and item 64's cooling top are one defect seen from two sides, and **every
  prognostic-column number in this file (items 45-47's 2927 W/m² and "still falling") was
  produced with that solver underneath** — redo them rather than extend them. Not the whole
  story: converged, the top still falls 266 → 229 K in one call.
- **The supersaturation at 236/256 km is `q_sat` FALLING, not `q_v` rising** (item 64,
  instrumented with `ATM_SAT_TRACE=1`). `SaturationAdjustment` reaches ~100.5 % RH within one
  call — it condenses 0.47 kg/kg at level 38 — so it is not stalling; and the prescribed profile
  HOLDS THE RATIO DOWN rather than causing it. Entry `q_v/q_sat` at level 38 over five
  iterations: **7.62 → 16.5 → 10.8** prescribed, against **7.62 → 976 → 1934** under
  `ATM_PROGNOSTIC_T=1`, where the cell free-runs to 220–228 K and `q_sat` collapses two orders of
  magnitude while `q_v` stays flat to six digits. Two things follow: the loop's exit test is on
  the STEP, not the residual, so it cannot tell convergence from a damped oscillation; and
  `alpha_entry` applies the −37 °C ice threshold a SECOND time as a master gain on all five
  write-backs (`ATM_SAT_NO_ALPHA=1`, default off) — a real defect, but ungating it makes the
  ratio worse (2094 → 2977), so it is not the lever. **Re-point item 16 at whatever cools the
  free-running top.**
- **The convective triggers are fractions of surface pressure now (item 63), and the convective
  layer is STILL one grid level thick — because invariant 2 pins it, not the triggers.** The
  thresholds were 1000/970/900/800/700/200 hPa, absolute Earth surface pressures, so every scan ran
  in the top few levels of a 250 bar column; they are `f_stat_* * p_stat.x[0][j][k]` now, ported from
  ATHAD_COND which had the repair from the day it was forked. **Measured: nothing moves** — cloud base
  level 37, LFS 38, the same 99 of 173 columns, the same 7 active cells, `c_u` still 0. The cloud-base
  test is `q_v_u >= scale*q_sat` and this column cannot saturate below 236.4 km, so no threshold can
  put a base lower. **Item 61's "the blocker is the trigger levels" is withdrawn**: ATHAD's convective
  layer is thin because its CONDENSING layer is thin. The same repair is load-bearing in ATHAD_COND
  (base at the sea, LFS at 20-28 km, 47 % of columns convecting), and it is there that
  `ATM_MC_GEOPOTENTIAL` finally produces updraft condensation — `max c_u` 0 -> 7.19e-04 g/kg/s at
  24.9 km. **The geopotential needed an ascent to act on, and ATHAD has none.**
- `time_start/end/step` remain because the time-slice loop is still structural, though only
  one slice ever runs.
