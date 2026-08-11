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
   Physics at this point is still the inherited Earth physics.

*(Further entries are added as each phase is measured.)*
