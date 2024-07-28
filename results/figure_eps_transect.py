import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import pandas as pd
import scipy.stats as ss
from matplotlib.markers import MarkerStyle
import cmocean
import warnings
# Suppress specific RuntimeWarning related to mean of empty slice
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*Mean of empty slice.*")

plt.rcParams.update({
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "font.size": 9
})

ONE_COLUMN_WIDTH = 8.3
TWO_COLUMN_WIDTH = 12
GOLDEN_RATIO = 1.61
cm = 1/2.54  # centimeters in inches

# read thorpe results data
thorpe_eps_df = pd.read_pickle("../scripts/thorpe_scales/method_data/Thorpe_eps_df_with_mab.pkl")
thorpe_mab = thorpe_eps_df.index
thorpe_gamma_n_df = pd.read_pickle("../scripts/thorpe_scales/method_data/Thorpe_neutral_density_df_with_mab.pkl")

# add background dissipation, but only where there is temperature data
BACKGROUND_DISSIPATION = 1e-10 #value taken from Hirano et al 2015
thorpe_eps_df.fillna(value=BACKGROUND_DISSIPATION, inplace=True)
thorpe_eps_df.where(cond=~thorpe_gamma_n_df.isna(), other=np.nan, inplace=True)

thorpe_lons = thorpe_eps_df.columns.to_numpy()
max_lon = max(thorpe_lons)
min_lon = min(thorpe_lons)


# half a degree bins
BIN_EDGES = np.arange(min_lon-1e-3*min_lon, 0.5+max_lon+1e-3*max_lon, 0.5)

# depth-level-wise (row-wise) arithmetic averaging
rows = []
for index, row in thorpe_eps_df.iterrows():
    values = row.to_numpy()
    bin_means = ss.binned_statistic(x=thorpe_lons, values=values, statistic=np.nanmean, bins=BIN_EDGES)[0]
    new_eps = bin_means
    new_row = pd.DataFrame([new_eps], columns=BIN_EDGES[:-1])
    rows.append(new_row)
binned_thorpe_eps_df = pd.concat(rows, sort=False).reset_index(drop=True)

rows = []
for index, row in thorpe_gamma_n_df.iterrows():
    values = row.to_numpy()
    bin_means = ss.binned_statistic(x=thorpe_lons, values=values, statistic=np.nanmean, bins=BIN_EDGES)[0]
    new_row = pd.DataFrame([bin_means], columns=BIN_EDGES[:-1])
    rows.append(new_row)

binned_thorpe_gamma_n_df = pd.concat(rows, sort=False).reset_index(drop=True)


#-------------------------------------------------------------------
# read eps_IGW results from IDEMIX method
eps_IGW_IDEMIX_df = pd.read_csv("../scripts/IDEMIX_parametrization/method_data/eps_IGW_IDEMIX_results.csv")

fig, ax = plt.subplots(1, figsize=(TWO_COLUMN_WIDTH*cm, 0.8*TWO_COLUMN_WIDTH*cm))

######################################################################## 
######################################################################## 
##################### Axis 0 ###########################################
######################################################################## 
"""
def fexp(f):
    return int(np.floor(np.log10(abs(f)))) if f != 0 else 0
def fman(f):
    return f/10**fexp(f)
def generate_pattern(start, end):
    num_steps = (end - start) * 2  # Each step/order of magnitude has a multiplier of 1 and 5 e.g. 5 \times 10^{-8}
    return [multiplier * 10**(start + i // 2) for i, multiplier in enumerate([1, 5] * (num_steps + 1)) if start + i // 2 <= end]

# Start and stop order of magnitude:
start_point = -10
end_point = -7
bounds = generate_pattern(start_point, end_point)[:-1]
ncolors = len(bounds) - 1
#cmap = plt.cm.get_cmap('viridis', ncolors)
norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors= 256)
"""
cmap = cmocean.cm.tempo
mpp = ax.pcolormesh(
    binned_thorpe_eps_df.columns, 
    thorpe_mab, 
    binned_thorpe_eps_df,
    norm=mcolors.LogNorm(vmin=1e-10, vmax=1e-7),
    cmap=cmap,
    rasterized=True
)
"""
cb = plt.colorbar(mpp, ax = ax, location = "right")#, aspect = 40, shrink = 0.8)
cb.set_label(r"Dissipation rate $\varepsilon\,$(W kg$^{-1}$)")
print(cb.ax.get_xticklabels(), len(cb.ax.get_xticklabels()))
assert cb.ax.get_xticklabels()
cb.ax.set_xticklabels([f'{fman(b):.0f}$\\times10^{{{fexp(b):.0f}}}$' for b in bounds])
"""

cb = plt.colorbar(mpp, ax=ax, location="top", extend = "max")  # , aspect=40, shrink=0.8)
cb.set_label(r"Dissipation rate $\varepsilon\,$(W$\,$kg$^{-1}$)")

water_mass_boundaries = [28.26, 28.40]  # + 28.00 boundary
gravity_current_boundary = [28.40]  # from Garabato et al 2002
CS = ax.contour(
    binned_thorpe_gamma_n_df.columns,
    thorpe_mab, 
    binned_thorpe_gamma_n_df,
    levels=water_mass_boundaries,
    linestyles=["dashed", "solid"],
    colors="k",
    linewidths=3,
    zorder=10
)
fmt = {}
strs = ['WSDW', 'WSBW']
for l, s in zip(CS.levels, strs):
    fmt[l] = s

# Label every other level using strings
clabels = ax.clabel(
    CS,
    CS.levels,
    inline=False,
    fmt=fmt,
    colors="black",
    fontsize=10,
    zorder=11
)
# adjust bboxes for better readability
[txt.set_bbox(dict(facecolor='lightgrey', alpha = 0.8, edgecolor='darkgrey', boxstyle="round", pad=0)) for txt in clabels]

ax.scatter(
    eps_IGW_IDEMIX_df["lon"],
    eps_IGW_IDEMIX_df["rounded_mab"],
    c=eps_IGW_IDEMIX_df["eps_IGW"],
    cmap=cmap,
    norm=mcolors.LogNorm(vmin=1e-10, vmax=1e-7),
    edgecolor="black",
    marker=MarkerStyle("o"),
    s=300,
    zorder=10
)

ax.set_facecolor('lightgrey')
ax.set_ylabel("Meters above bottom")
ax.set_xlabel("Longitude (°)")

# for the legend
# eps_IGW icon
ax.scatter(
    eps_IGW_IDEMIX_df["lon"],
    eps_IGW_IDEMIX_df["rounded_mab"],
    #c=energy_levels["eps"],
    color="tab:gray",
    edgecolor="black",
    marker=MarkerStyle("o"),
    s=200,
    zorder=-10,
    label='$\\varepsilon_{\\mathrm{IGW, IDEMIX}}$',
)

# artificial eps_total icon
ax.scatter(
    eps_IGW_IDEMIX_df["lon"],
    eps_IGW_IDEMIX_df["rounded_mab"],
    color="tab:gray",
    edgecolor="black",
    marker=MarkerStyle("s"),
    s=80,
    zorder=-10,
    label="$\\varepsilon_{\\mathrm{total, Thorpe}}$"
)


ax.legend(loc = "upper left", ncol=3, columnspacing=1).set_zorder(50)
# ax.annotate('gravity current\nboundary', xy=(-48.8, 130), xytext=(-48.5, 270), #fontsize=9,
#             arrowprops = dict(facecolor='black', width = 2, shrink=0.05), ha = "center", va = "center", color = "white", bbox=dict(facecolor='black', alpha = 0.8, edgecolor='black', boxstyle='round, pad = 0.5'))

ax.set_ylim(-10,500)
fig.tight_layout()
fig.savefig("./eps_transect.pdf")
#fig.savefig("./eps_transect.png", dpi = 400, bbox_inches = "tight")
plt.show()