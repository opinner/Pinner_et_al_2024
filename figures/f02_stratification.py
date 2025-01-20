import numpy as np
import pandas as pd
import scipy.stats as ss
import matplotlib
import matplotlib.pyplot as plt
import cmocean as cmocean
import warnings
# Suppress specific RuntimeWarning related to mean of empty slice
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*Mean of empty slice.*")

ONE_COLUMN_WIDTH = 8.3
TWO_COLUMN_WIDTH = 12
GOLDEN_RATIO = 1.61
cm = 1/2.54  # centimeters in inches

plt.style.use('./paper.mplstyle')

binned_neutral_density = pd.read_csv("../derived_data/binned_neutral_density.csv", index_col = 0)
binned_neutral_density.columns = binned_neutral_density.columns.astype("float") #convert column names from strings to floats

fig, ax = plt.subplots(1, figsize=(TWO_COLUMN_WIDTH*cm, 0.5*TWO_COLUMN_WIDTH*cm))
mpp = ax.pcolormesh(
    binned_neutral_density.columns,
    binned_neutral_density.index,
    binned_neutral_density.values,
    cmap=cmocean.cm.haline_r,
    vmin=27.8,
    rasterized=True # optimize the drawing for vector graphics
)

cb = plt.colorbar(mpp, ax=ax, extend="min", location="top")  # pad=0.02, aspect=12

water_mass_boundaries = [28.26, 28.40]  # + 28.00 boundary
# gravity_current_boundary = [28.40]  # from Garabato et al 2002
CS = ax.contour(
    binned_neutral_density.columns,
    binned_neutral_density.index,
    binned_neutral_density.values,
    levels=water_mass_boundaries,
    linestyles=["dashed", "solid"],
    colors="k",
    linewidths=3,
)

fmt = {}
strs = ['WSDW', 'WSBW']
for l, s in zip(CS.levels, strs):
    fmt[l] = s

# Label every other level using strings
ax.clabel(
    CS,
    CS.levels,
    inline=False,
    fmt=fmt,
    fontsize=10,
    colors = "white"
)


# ax.annotate('bottom\ncurrent', xy=(-51.69, 184), xytext=(-51.8, 230),
#              arrowprops=dict(facecolor='black', width=2, shrink=0.05), ha="center", color="white",
#              bbox=dict(facecolor='black', alpha=0.8, edgecolor='black', boxstyle='round'))
# ax.annotate('gravity current\nboundary', xy=(-48.8, 130), xytext=(-48, 270), #fontsize=9,
#             arrowprops = dict(facecolor='black', width = 2, shrink=0.05), ha = "center", va = "center", color = "white", bbox=dict(facecolor='black', alpha = 0.8, edgecolor='black', boxstyle='round, pad = 0.5'))

mooring_info = pd.read_csv("../scripts/IDEMIX_parameterization/method_results/eps_IGW_IDEMIX_results.csv")
moorings_mabs = mooring_info["rounded mab"]
moorings_lons = mooring_info["lon"]

# draw measurement positions
ax.plot(moorings_lons, moorings_mabs,
        "D",
        label="rotor current\nmeters",
        color="tab:red",
        markersize=10,
        markeredgecolor="k",
        zorder=6)

ax.set_ylim(-10, 500)
xlim = ax.get_xlim()
ax.set_xlim((xlim[0] - 0.2, xlim[1] + 0.8))
ax.set_xlabel("Longitude (°)")
ax.set_ylabel("Meters above Seafloor")
cb.set_label(r"Neutral Density $\gamma^\text{n}\,$(kg$\,$m$^{-3}$)")
cb.ax.plot([water_mass_boundaries[0], water_mass_boundaries[0]], [0, 1], 'k--', lw=2)
cb.ax.plot([water_mass_boundaries[1], water_mass_boundaries[1]], [0, 1], 'k', lw=2, ls = "solid")

# draw sea floor
x = list(ax.get_xlim())
y1 = [0, 0]
y2 = 2 * list(ax.get_ylim())[0]
ax.axhline(0, c="k")
ax.fill_between(x, y1, y2, facecolor="xkcd:charcoal grey", zorder=5)  # , hatch="///")
ax.legend(loc = "upper right")  #,facecolor='k', framealpha=0.8, edgecolor = "black", labelcolor = "white")
fig.tight_layout()
#fig.savefig("./stratification.pdf")
# fig.savefig("./stratification.png", dpi = 400, bbox_inches = "tight")
plt.show()
