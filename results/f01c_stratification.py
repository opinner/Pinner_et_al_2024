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

thorpe_gamma_n_df = pd.read_pickle("../scripts/thorpe_scales/method_results/Thorpe_neutral_density_df_with_mab.pkl")
lons = thorpe_gamma_n_df.columns.to_numpy()
mab = thorpe_gamma_n_df.index
max_lon = max(lons)
min_lon = min(lons)
# half a degree bins
BIN_EDGES = np.arange(min_lon - 1e-3 * min_lon, 0.5+max_lon + 1e-3 * max_lon, 0.5)

rows = []
for index, row in thorpe_gamma_n_df.iterrows():
    values = row.to_numpy()
    bin_means = ss.binned_statistic(x=lons, values=values, statistic=np.nanmean, bins=BIN_EDGES)[0]
    new_row = pd.DataFrame([bin_means], columns=BIN_EDGES[:-1])
    rows.append(new_row)
binned_thorpe_gamma_n_df = pd.concat(rows, sort=False).reset_index(drop=True)

fig, ax = plt.subplots(1, figsize=(TWO_COLUMN_WIDTH*cm, 0.5*TWO_COLUMN_WIDTH*cm))
mpp = ax.pcolormesh(
    binned_thorpe_gamma_n_df.columns,
    mab,
    binned_thorpe_gamma_n_df,
    cmap=cmocean.cm.haline_r,
    vmin = 27.8,
    rasterized = True # optimize the drawing for vector graphics
)

cb = plt.colorbar(mpp, ax=ax, extend="min", location="top")  # pad=0.02, aspect=12

water_mass_boundaries = [28.26, 28.40]  # + 28.00 boundary
# gravity_current_boundary = [28.40]  # from Garabato et al 2002
CS = ax.contour(
    binned_thorpe_gamma_n_df.columns,
    binned_thorpe_gamma_n_df.index,
    binned_thorpe_gamma_n_df,
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

mooring_info = pd.read_csv("../scripts/IDEMIX_parametrization/method_results/eps_IGW_IDEMIX_results.csv")
moorings_mabs = mooring_info["rounded_mab"]
moorings_lons = mooring_info ["lon"]

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
cb.set_label(r"Neutral Density $\gamma_n\,$(kg$\,$m$^{-3}$)")
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
fig.savefig("./stratification.svg")
# fig.savefig("./stratification.png", dpi = 400, bbox_inches = "tight")
plt.show()
