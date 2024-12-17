import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import pandas as pd
import scipy.stats as ss
from matplotlib.markers import MarkerStyle
import cmocean
cmap = cmocean.cm.amp
import warnings
# Suppress specific RuntimeWarning related to mean of empty slice
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*Mean of empty slice.*")

ONE_COLUMN_WIDTH = 8.3
TWO_COLUMN_WIDTH = 12
GOLDEN_RATIO = 1.61
cm = 1/2.54  # centimeters in inches

plt.rcParams.update({
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "font.size": 9
})

# -------------------------------------------------------------------
# read neutral density data
thorpe_gamma_n_df = pd.read_pickle("../scripts/thorpe_scales/method_results/Thorpe_neutral_density_df_with_mab.pkl")
thorpe_mab = thorpe_gamma_n_df.index

thorpe_lons = thorpe_gamma_n_df.columns.to_numpy()
max_lon = max(thorpe_lons)
min_lon = min(thorpe_lons)

# half a degree bins
BIN_EDGES = np.arange(min_lon-1e-3*min_lon, 0.5+max_lon+1e-3*max_lon, 0.5)
BIN_CENTER = BIN_EDGES[:-1]-0.25

# depth-level-wise (row-wise) arithmetic averaging
rows = []
for index, row in thorpe_gamma_n_df.iterrows():
    values = row.to_numpy()
    bin_means = ss.binned_statistic(x=thorpe_lons, values=values, statistic=np.nanmean, bins=BIN_EDGES)[0]
    new_row = pd.DataFrame([bin_means], columns=BIN_CENTER)
    rows.append(new_row)

binned_thorpe_gamma_n_df = pd.concat(rows, sort=False).reset_index(drop=True)


#-------------------------------------------------------------------
# read eps_IGW results from strain-based finestructure analysis
eps_IGW_strain_df = pd.read_csv("../scripts/shear_strain_parametrization/method_results/binned_strain_eps2.csv")
eps_IGW_strain_df.set_index('Unnamed: 0', inplace=True)

#-------------------------------------------------------------------
# read eps_IGW results from IDEMIX method
eps_IGW_IDEMIX_df = pd.read_csv("../scripts/IDEMIX_parametrization/method_results/eps_IGW_IDEMIX_results.csv")

fig,ax = plt.subplots(1, figsize=(TWO_COLUMN_WIDTH*cm, 0.8*TWO_COLUMN_WIDTH*cm), layout="constrained")
######################################################################## 
######################################################################## 
##################### Axis 0 ###########################################
########################################################################
norm = mcolors.LogNorm(vmin=1e-10, vmax=1e-8)
mpp = ax.pcolormesh(
    eps_IGW_strain_df.columns.astype(float),
    eps_IGW_strain_df.index,
    eps_IGW_strain_df,
    cmap=cmap,
    norm=norm,
    shading="nearest"
)
      
cb = plt.colorbar(mpp, ax=ax, location="top", extend="max", aspect=26)
cb.set_label(r"Wave-induced dissipation rate $\varepsilon_{\mathrm{IGW}}\,$(W$\,$kg$^{-1}$)")

# continental_slope = BIN_EDGES[3]
# deep_sea = BIN_EDGES[-1]
# ax.fill_between(
#     [continental_slope, deep_sea],
#     [0,0],
#     [250/2, 250/2],
#     hatch="xx",
#     facecolor='None',
#     edgecolor='darkgrey',
#     alpha=0.5
# )

binned_regions = pd.read_csv("../scripts/preprocessing/method_results/binned_regions.csv", index_col=0)
binned_regions.columns = binned_regions.columns.astype("float") #convert column names from strings to floats
binned_regions = binned_regions.iloc[0:600]
levels = [2.5,3.5]  # Border between IL and BL
plt.rcParams['hatch.color'] = 'xkcd:charcoal'
ax.contourf(
    binned_regions.columns,
    binned_regions.index,
    binned_regions.values,
    levels=levels,
    hatches=["xx"],
    colors="None",
    zorder=1
)
ax.contour(
    binned_regions.columns,
    binned_regions.index,
    binned_regions.values,
    levels=levels,
    colors="xkcd:charcoal",
    zorder=1
)

ax.scatter(
    eps_IGW_IDEMIX_df["lon"],
    eps_IGW_IDEMIX_df["rounded mab"],
    c=eps_IGW_IDEMIX_df["eps_IGW"],
    cmap=cmap,
    norm=norm,
    edgecolor='darkgrey',
    marker=MarkerStyle("o"),
    s=275,
    zorder=5
)

# ------------------------------------------------------------------------------------------
# for the legend
# eps_IGW IDEMIX icon
ax.scatter(
    eps_IGW_IDEMIX_df["lon"],
    eps_IGW_IDEMIX_df["rounded mab"],
    #c=energy_levels["eps"],
    color = "tab:gray",
    edgecolor="black",
    marker=MarkerStyle("o"),
    s=175,
    zorder=-10,
    label=r"$\varepsilon_{\mathrm{IGW, IDEMIX}}$",
)

# artificial eps_IGW strain icon
ax.scatter(
    eps_IGW_IDEMIX_df["lon"],
    eps_IGW_IDEMIX_df["rounded mab"],
    color = "tab:gray",
    edgecolor="black",
    marker=MarkerStyle("s"),
    s=70,
    zorder=-10,
    label=r"$\varepsilon_{\mathrm{IGW, fine}}$"
)


water_mass_boundaries = [28.26, 28.40]  # + 28.00 boundary
# gravity_current_boundary = [28.40]  # from Garabato et al 2002
CS = ax.contour(
    binned_thorpe_gamma_n_df.columns,
    binned_thorpe_gamma_n_df.index,
    binned_thorpe_gamma_n_df,
    levels=water_mass_boundaries,
    linestyles=["dashed", "solid"],
    colors="k",
    linewidths=2,
    zorder=2
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
    colors=["white", "black"]
)


# velocity isolines
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
    linewidths=2,
    alpha=0.8
)

#ax.set_ylim(-10, 500)
ax.set_ylim(-10,600)
ax.legend(loc = "upper left", ncol=3, columnspacing=1)
#ax.annotate('gravity current\nboundary', xy=(-48.8, 130), xytext=(-48.5, 270), #fontsize=9,
#            arrowprops = dict(facecolor='black', width = 2, shrink=0.05), ha = "center", va = "center", color = "white", bbox=dict(facecolor='black', alpha = 0.8, edgecolor='black', boxstyle='round, pad = 0.5'))

ax.set_facecolor('lightgrey')
ax.set_ylabel("Meters above bottom")
ax.set_xlabel("Longitude (°)")
#ax[0].set_title(r"Dissipation rate $\varepsilon$ across the slope")

#fig.tight_layout()
fig.savefig("./eps_IGW_comparison.pdf")
# fig.savefig("./eps_IGW_comparison.png", dpi = 400, bbox_inches = "tight")
plt.show()