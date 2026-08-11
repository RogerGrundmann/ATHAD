# ATHAD — Atmosphere of the Earth in the Hadean Eon

An atmospheric general-circulation model of the Earth as it was in the Hadean
(~4.4 Ga): a ~250 bar, water-vapour-dominated atmosphere over a molten or
quenching surface with no known topography.

**Forked from `ATOM_Precipitation` @ `1e3f319` (2026-07-28), atmosphere half only.**
Started 2026-08-11. There is no hydrosphere, no paleogeography, and no time-slice
series — ATHAD is one epoch.

## Build and run

```bash
make had
cd cli && OMP_NUM_THREADS=1 ./had config_athad.xml
```

`make` regenerates the parameter bindings from `param.py` first. The generated
files (`atmosphere/*.inc`, `python/atmosphere_pxd.pxi`, `python/pyathad.pyx`,
`cli/config_athad.xml`, `python/config_athad.xml`) are **tracked on purpose**, so a
`param.py` change shows its full effect in the diff. Regenerate and commit them
together.

## Atmospheric composition

Mole fractions are the input; the model works in mass fractions. The residual
7 % is split evenly across the five trace gases.

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
- **Background** (everything except H₂O and CO₂): x_bg = 0.100, q_bg = 0.1223,
  M_bg = 26.207 g/mol, **R_bg = 317.3 J/(kg·K)**. This is what ATHAD's `R_Air`
  parameter now means — it is no longer "air".
- **p_surf = 250 bar**, **T_surf ≈ 1500 K**, **ρ_surf = 43.0 kg/m³** (Earth: 1.2)
- Scale height R_mix·T/g: **59 km** at 1500 K, 24 km at 600 K. Reaching the
  ~0.1 bar radiating level takes ln(2500) ≈ 7.8 e-foldings, hence a shell
  ~300 km deep rather than Earth's 16 km.
- cp ≈ 2040 J/(kg·K) and strongly T-dependent across 300–1500 K — roughly 2×
  Earth's 1005. A constant cp misplaces the lapse rate everywhere.

Only **H₂O (`c`) and CO₂ (`co2`) are prognostic**. The other six gases are a
fixed well-mixed background entering R_mix, cp_mix and the opacity.

## Three invariants — do not silently break these

1. **There is no topography.** `h ≡ 0`, `i_topography ≡ 0`, `Topography ≡ 0`
   everywhere, so `AtomUtils::is_land()` is false at every point. The Hadean
   surface is unknown; a featureless global surface is the deliberate choice, not
   a missing data file. Do not reintroduce a bathymetry read, and do not "fix"
   `is_land()` — the land branches are dead by construction and that is correct.

2. **Water is supercritical below the condensation level.** The critical point is
   647.096 K / 220.64 bar; p_H₂O at the surface is 0.8 × 250 = 200 bar at 1500 K.
   The inherited Magnus formula (`hp·exp_func(T, 17.2694, 35.86)`) is capped at
   ~101 °C and is invalid here — use the IAPWS curve in `SaturationH2O.h`. The
   dilute approximation `q_sat = ep·E/(p − (1−ep)E)` is also invalid, because H₂O
   *is* the bulk gas, not a trace: use the exact mass-fraction form. Every
   condensation path must be guarded on `T < T_crit`, and be a genuine no-op below
   the condensation level rather than a clamped Magnus value. `lv` is not a
   constant — it vanishes at the critical point.

3. **Radiation runs in mode 2 (direct σT⁴).** The inherited default was mode 5,
   Newtonian relaxation toward the Scotese paleo-temperature curve. **No Scotese
   curve exists for 4.4 Ga.** The empirical emissivity fits ATHAD inherited
   (Bignami 1995 `eps_dry = 0.684 + 0.0056·e_surf`, a Mediterranean sea-surface
   regression; CO₂ via Atwater & Ball) are calibrated on present-day Earth columns
   and are meaningless at 250 bar — they are replaced by a column-mass optical
   depth with pressure broadening.

## Relationship to the family

Sibling models live beside this directory: `ATOM_Precipitation` (modern Earth),
`ATJUP`, `ATSAT`, `ATURAN`, `ATNEPT` (giants), `ASTIM` (impacts).

C++ class, file and function names are kept **identical to `ATOM_Precipitation`**
(`cAtmosphereModel`, `ThermoAtm.h`, `RHS_Atm_Turb.cpp`, …) so fixes can be
cherry-picked in both directions; only the outer shell is renamed (`libathad.a`,
`cli/had`, `config_athad.xml`, `pyathad`). Preserve that property.

Known traps already solved elsewhere in the family — check these before
re-deriving:

- Coriolis / centrifugal sign and projection errors: ATURAN `8b284cb`, `4201957`;
  ATNEPT `024c37f`, `e412b1b`.
- Mixture properties must be weighted by **mass** fraction, not mole: ATNEPT `c116d71`.
- In-place Gauss–Seidel in the pressure solver is a threading defect: ATURAN `ffd0e0e`.
- Report failures and limits in the README, including the measurements that did
  not work out: ATURAN `74b4ded`, ATNEPT `34286b8`.

## Stated assumptions, open for revision

These were chosen to get a running model and are not settled results:

- **ω = 3.17e-4 rad/s** (a 5.5 h Hadean day, 4.35× modern)
- **Insolation 0.71 S₀ = 966 W/m²** (faint young Sun at 4.4 Ga)
- Albedo of a global cloud deck over a molten surface (~0.3–0.5)
- Geothermal / magma-ocean bottom heat flux — at 1500 K plausibly comparable to
  or larger than absorbed solar, so it cannot be omitted
- Tropopause height, currently inherited from Earth and certainly wrong
- The cp(T) polynomial fit range

## Known open risks

- **Boussinesq.** The solver rests on the Boussinesq buoyancy approximation, but
  density varies by ~2 orders of magnitude across a 250 bar column. This may force
  an anelastic or fully compressible formulation. The family's existing partial
  answer is the ATJUP hydrostatic split (ported in ATURAN `302a51e`) — and it did
  not cure the giants' problem.
- **Radiation is the project.** The domain, composition and surface changes are
  mechanical; the emissivity rewrite decides whether ATHAD produces a meaningful
  Hadean climate. The sharpest single check is OLR at the top of the atmosphere:
  it should sit near the **280–310 W/m² Nakajima / Komabayashi–Ingersoll runaway
  limit**, not the ~287 kW/m² that σT⁴(1500 K) through a transparent atmosphere
  would give.
- Inherited from the parent (`ATOM_Precipitation/docs/co2_sensitivity.md`):
  radiation mode 2 gives the correct-sign hydrological response but a
  structurally weak temperature response on Earth. Hadean forcing is ~250× larger,
  so this may not bite here — but it is not yet demonstrated.
