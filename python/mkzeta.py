#!/usr/bin/env python3
"""Generate one config per zeta for the zeta scan at im = 41.

The shell is held at 300 km (L_atm = 300000/(exp(zeta)-1)) and dt_visc is held FIXED, so
every run is at the same physical time per iteration and the comparison isolates zeta.
Holding dt fixed is also the safe direction: lower zeta widens the surface spacing, so the
diffusive CFL only gets looser. The 13x-larger-dt claim is a separate run.

usage: mkzeta.py <nm> <zeta> [<zeta> ...]
"""
import math, sys, re

SHELL = 300.0e3

def make(nm, zeta):
    L_atm = SHELL/(math.exp(zeta) - 1.0)
    tag   = f"z{zeta:g}".replace('.', 'p')
    src   = open("config_athad.xml").read()
    out   = src
    subs = [
        (r"<output_path>[^<]*</output_path>", f"<output_path>output_zeta_{tag}/</output_path>"),
        (r"<nm>[0-9]*</nm>",                  f"<nm>{nm}</nm>"),
        (r"<zeta>[0-9.eE+-]*</zeta>",         f"<zeta>{zeta}</zeta>"),
        (r"<L_atm>[0-9.eE+-]*</L_atm>",       f"<L_atm>{L_atm:.4f}</L_atm>"),
        # moist physics OFF: coeff_MC_*/coeff_L ride on L_atm, which varies 29x across this scan
        (r"<moist_phys_start_iter>[0-9]*</moist_phys_start_iter>",
                                              f"<moist_phys_start_iter>{nm+1}</moist_phys_start_iter>"),
        # kill the heavy I/O: 925 MB restart dumps and the panoramas
        (r"<restart_stride>[0-9-]*</restart_stride>",   "<restart_stride>0</restart_stride>"),
        (r"<panorama_print>[0-9-]*</panorama_print>",   f"<panorama_print>{nm+1}</panorama_print>"),
        (r"<checkpoint>[0-9]*</checkpoint>",            f"<checkpoint>{nm}</checkpoint>"),
        (r"<diagnostic_stride>[0-9]*</diagnostic_stride>", "<diagnostic_stride>20</diagnostic_stride>"),
    ]
    for pat, rep in subs:
        out, n = re.subn(pat, rep, out)
        if n != 1:
            raise SystemExit(f"pattern {pat!r} matched {n} times")
    fn = f"config_zeta_{tag}.xml"
    open(fn, "w").write(out)
    print(f"{fn}  zeta={zeta}  L_atm={L_atm:.1f} m  shell={(math.exp(zeta)-1)*L_atm/1e3:.1f} km")

if __name__ == "__main__":
    nm = int(sys.argv[1])
    for z in sys.argv[2:]:
        make(nm, float(z))
