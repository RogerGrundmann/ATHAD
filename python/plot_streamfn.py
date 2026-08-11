# Plot the meridional mass streamfunction Psi(lat, z) for two iterations side by side
# from the write_meridional_streamfunction() CSV dumps, on a shared color scale, so the
# vertical (depth) collapse of the overturning is directly visible.
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "output_0Ma"
iters = [20, 100]

def load(it):
    df = pd.read_csv(f"{OUT}/meridional_streamfunction_{it}.csv")
    lat = np.sort(df["lat_deg"].unique())          # ascending: 90S -> 90N
    z   = np.sort(df["height_m"].unique())         # ascending: surface -> lid
    P = df.pivot(index="height_m", columns="lat_deg", values="psi_kg_per_s")
    P = P.reindex(index=z, columns=lat)
    return lat, z, P.values / 1.0e9                # -> 1e9 kg/s

data = [load(it) for it in iters]
vmax = max(np.nanmax(np.abs(P)) for _, _, P in data)
levels = np.linspace(-vmax, vmax, 21)

fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
for ax, it, (lat, z, P) in zip(axes, iters, data):
    cf = ax.contourf(lat, z / 1000.0, P, levels=levels, cmap="RdBu_r", extend="both")
    ax.contour(lat, z / 1000.0, P, levels=levels, colors="k", linewidths=0.3, alpha=0.5)
    ax.set_title(f"iter {it}")
    ax.set_xlabel("latitude [deg N]")
    ax.set_xlim(-90, 90)
    ax.set_xticks(np.arange(-90, 91, 30))
    ax.axvline(0, color="0.4", lw=0.6, ls="--")
axes[0].set_ylabel("height [km]")
cb = fig.colorbar(cf, ax=axes, shrink=0.9, pad=0.02)
cb.set_label(r"$\Psi$  [$10^{9}$ kg s$^{-1}$]")
fig.suptitle("Meridional mass streamfunction  Psi(lat, z)  -  Ma=0  (v/w drag cut 10x)", y=0.98)
out = f"{OUT}/streamfn_lat_z_iter{iters[0]}_vs_{iters[1]}.png"
fig.savefig(out, dpi=130, bbox_inches="tight")
print("wrote", out)

# also report the column-max overturning depth (highest level where |Psi| > 25% of peak)
for it, (lat, z, P) in zip(iters, data):
    thr = 0.25 * np.nanmax(np.abs(P))
    strong = np.abs(P) > thr
    zmax = z[np.where(strong.any(axis=1))[0].max()] if strong.any() else 0.0
    print(f"iter {it:3d}: peak|Psi|={np.nanmax(np.abs(P)):6.1f}e9 kg/s   "
          f"overturning reaches up to {zmax/1000.0:4.2f} km (|Psi|>25% of peak)")
