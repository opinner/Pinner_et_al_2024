import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import pandas as pd
from matplotlib.markers import MarkerStyle
import cmocean
import warnings

# Suppress specific RuntimeWarning related to mean of empty slice
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*Mean of empty slice.*")

#plt.style.use('./paper.mplstyle')

ONE_COLUMN_WIDTH = 8.3
TWO_COLUMN_WIDTH = 12
GOLDEN_RATIO = 1.61
cm = 1 / 2.54  # centimeters in inches

# read thorpe results data
thorpe_eps_df = pd.read_pickle("../scripts/thorpe_scales/method_results/Thorpe_eps_df_with_mab.pkl")
thorpe_mab = thorpe_eps_df.index
thorpe_gamma_n_df = pd.read_pickle("../scripts/thorpe_scales/method_results/Thorpe_neutral_density_df_with_mab.pkl")

# add background dissipation, but only where there is temperature data
BACKGROUND_DISSIPATION = 1e-10  #value taken from Hirano et al 2015
thorpe_eps_df.fillna(value=BACKGROUND_DISSIPATION, inplace=True)
# use gamma_n as a mask
thorpe_eps_df.where(cond=~thorpe_gamma_n_df.isna(), other=np.nan, inplace=True)

#but then use the already binned version
binned_thorpe_gamma_n_df = pd.read_csv(
    "../scripts/preprocessing/method_results/binned_gamma_n.csv", index_col=0)
#binned_thorpe_gamma_n_df.set_index(keys='Unnamed: 0', inplace=True)
binned_thorpe_gamma_n_df = binned_thorpe_gamma_n_df.drop(
    binned_thorpe_gamma_n_df[binned_thorpe_gamma_n_df.index > 600].index)
binned_thorpe_gamma_n_df.columns = binned_thorpe_gamma_n_df.columns.astype(pd.Float64Dtype())

# # bin dissipation rates
# thorpe_lons = thorpe_eps_df.columns.to_numpy()
# max_lon = max(thorpe_lons)
# min_lon = min(thorpe_lons)
#
# # half a degree bins
# BIN_EDGES = np.arange(min_lon - 1e-3 * min_lon, 0.5 + max_lon + 1e-3 * max_lon, 0.5)
# BIN_CENTER = BIN_EDGES[:-1]-0.25
# # depth-level-wise (row-wise) arithmetic averaging
# rows = []
# for index, row in thorpe_eps_df.iterrows():
#     values = row.to_numpy()
#     bin_means = ss.binned_statistic(x=thorpe_lons, values=values, statistic=np.nanmean, bins=BIN_EDGES)[0]
#     new_eps = bin_means
#     new_row = pd.DataFrame([new_eps], columns=BIN_CENTER)
#     rows.append(new_row)
# binned_thorpe_eps_df = pd.concat(rows, sort=False).reset_index(drop=True)
# binned_thorpe_eps_df.to_csv("../derived_data/binned_thorpe_dissipation.csv")

binned_thorpe_eps_df = pd.read_csv(
    "../scripts/thorpe_scales/method_results/binned_thorpe_dissipation.csv", index_col=0)
#-------------------------------------------------------------------
# read eps_IGW results from IDEMIX method
eps_IGW_IDEMIX_df = pd.read_csv("../scripts/IDEMIX_parametrization/method_results/eps_IGW_IDEMIX_results.csv")

fig, ax = plt.subplots(1, figsize=(TWO_COLUMN_WIDTH * cm, 0.8 * TWO_COLUMN_WIDTH * cm), layout="constrained")

######################################################################## 
######################################################################## 
##################### Axis 0 ###########################################
######################################################################## 

cmap = cmocean.cm.tempo
mpp = ax.pcolormesh(
    binned_thorpe_eps_df.columns,
    binned_thorpe_eps_df.index,
    binned_thorpe_eps_df.values,
    norm=mcolors.LogNorm(vmin=1e-10, vmax=1e-7),
    cmap=cmap,
    rasterized=True
)

cb = plt.colorbar(mpp, ax=ax, location="top", extend="max", aspect=26)  # , aspect=40, shrink=0.8)
cb.set_label(r"Dissipation rate $\varepsilon\,$(W$\,$kg$^{-1}$)")

water_mass_boundaries = [28.26, 28.40]  # + 28.00 boundary
gravity_current_boundary = [28.40]  # from Garabato et al 2002

CS = ax.contour(
    binned_thorpe_gamma_n_df.columns,
    binned_thorpe_gamma_n_df.index,
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
[txt.set_bbox(dict(facecolor='lightgrey', alpha=0.8, edgecolor='darkgrey', boxstyle="round", pad=0)) for txt in clabels]

ax.scatter(
    eps_IGW_IDEMIX_df["lon"],
    eps_IGW_IDEMIX_df["rounded mab"],
    c=eps_IGW_IDEMIX_df["eps_IGW"],
    cmap=cmap,
    norm=mcolors.LogNorm(vmin=1e-10, vmax=1e-7),
    edgecolor="black",
    marker=MarkerStyle("o"),
    s=300,
    zorder=10
)

# for the legend
# eps_IGW icon
ax.scatter(
    eps_IGW_IDEMIX_df["lon"],
    eps_IGW_IDEMIX_df["rounded mab"],
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
    eps_IGW_IDEMIX_df["rounded mab"],
    color="tab:gray",
    edgecolor="black",
    marker=MarkerStyle("s"),
    s=80,
    zorder=-10,
    label="$\\varepsilon_{\\mathrm{total, Thorpe}}$"
)

with np.load("./flowfield.npz") as data:
    xi = data['xi']
    yi = data['yi']
    avrg_v = data['avrg_v']

levels = [0.1, 0.2]
CS = ax.contour(
    xi, yi, avrg_v,
    levels=levels,
    colors='yellow',
    linestyles=["dashed", "solid"],
    linewidths=3,
    alpha=0.8
)

# fmt = {}
# for l, s in zip(CS.levels, levels):
#     fmt[l] = f"{s:.1f}"
#
# clabels = ax.clabel(
#     CS,
#     CS.levels,
#     inline=True,
#     fmt=fmt,
#     fontsize=10,
#     colors='yellow',
#     zorder=11
# )

ax.legend(loc="upper left", ncol=3, columnspacing=1).set_zorder(50)
ax.set_ylim(-10, 500)
#ax.xlim()
ax.set_facecolor('lightgrey')
ax.set_ylabel("Meters above bottom")
ax.set_xlabel("Longitude (°)")

fig.savefig("./eps_transect.pdf")
#fig.savefig("./eps_transect.png", dpi = 400, bbox_inches = "tight")
plt.show()
