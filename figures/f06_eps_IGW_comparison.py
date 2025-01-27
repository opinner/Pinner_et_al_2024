import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import pandas as pd
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

binned_neutral_density = pd.read_csv("../derived_data/binned_neutral_density.csv", index_col = 0)
binned_neutral_density.columns = binned_neutral_density.columns.astype("float") #convert column names from strings to floats

#-------------------------------------------------------------------
# read eps_IGW results from strain-based finestructure analysis
eps_IGW_strain_df = pd.read_csv("../scripts/finestructure/method_results/binned_strain_eps.csv")
eps_IGW_strain_df.set_index('Unnamed: 0', inplace=True)

#-------------------------------------------------------------------
# read eps_IGW results from IDEMIX method
eps_IGW_IDEMIX_df = pd.read_csv("../scripts/IDEMIX_parameterization/method_results/eps_IGW_IDEMIX_results.csv")

fig,ax = plt.subplots(1, figsize=(TWO_COLUMN_WIDTH*cm, 0.8*TWO_COLUMN_WIDTH*cm), layout="constrained")
######################################################################## 
######################################################################## 
##################### Axis 0 ###########################################
########################################################################
norm = mcolors.LogNorm(vmin=1e-10, vmax=1e-8)
mpp = ax.pcolormesh(
    -eps_IGW_strain_df.columns.astype(float),
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
levels = [2.5, 3.5]  # Border between IL and BL
plt.rcParams['hatch.color'] = 'xkcd:charcoal'
ax.contourf(
    -binned_regions.columns,
    binned_regions.index,
    binned_regions.values,
    levels=levels,
    hatches=["xx"],
    colors="None",
    zorder=10
)
ax.contour(
    -binned_regions.columns,
    binned_regions.index,
    binned_regions.values,
    levels=levels,
    colors="xkcd:charcoal",
    zorder=10
)

ax.scatter(
    -eps_IGW_IDEMIX_df["lon"],
    eps_IGW_IDEMIX_df["rounded mab"],
    c=eps_IGW_IDEMIX_df["eps_IGW"],
    cmap=cmap,
    norm=norm,
    edgecolor='darkgrey',
    marker=MarkerStyle("o"),
    s=275,
    zorder=50
)

# ------------------------------------------------------------------------------------------
# for the legend
# eps_IGW IDEMIX icon
ax.scatter(
    -eps_IGW_IDEMIX_df["lon"],
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
    -eps_IGW_IDEMIX_df["lon"],
    eps_IGW_IDEMIX_df["rounded mab"],
    color="tab:gray",
    edgecolor="black",
    marker=MarkerStyle("s"),
    s=70,
    zorder=-10,
    label=r"$\varepsilon_{\mathrm{IGW, fine}}$"
)


water_mass_boundaries = [28.26, 28.40]  # + 28.00 boundary
# gravity_current_boundary = [28.40]  # from Garabato et al 2002
CS = ax.contour(
    -binned_neutral_density.columns,
    binned_neutral_density.index,
    binned_neutral_density,
    levels=water_mass_boundaries,
    linestyles=["dashed", "solid"],
    colors=["white", "black"],
    linewidths=2,
    zorder=2
)

# to be shifted in postprocessing
strs = ['WSDW', 'WSBW', "IL", "BL"]
colors = ["white", "black", "black", "xkcd:charcoal"]
for ix, (s,color) in enumerate(zip(strs,colors)):
    ax.text(0.9, 0.9-0.05*ix,s, color=color, fontsize=8, transform=ax.transAxes)

# fmt = {}
# strs = ['WSDW', 'WSBW']
# for l, s in zip(CS.levels, strs):
#     fmt[l] = s
#
# # Label every other level using strings
# ax.clabel(
#     CS,
#     CS.levels,
#     inline=False,
#     fmt=fmt,
#     fontsize=8,
#     colors=["white", "black"]
# )

# velocity isolines
with np.load("./flowfield.npz") as data:
    xi = data['xi']
    yi = data['yi']
    avrg_v = data['avrg_v']

levels = [0.1, 0.2]
CS = ax.contour(
    -xi, yi, avrg_v,
    levels=levels,
    colors='yellow',
    linestyles=["dashed", "solid"],
    linewidths=2,
    alpha=0.8
)

# to be shifted in postprocessing
for ix, s in enumerate(levels):
    label = f"{s:.1f}" +r"$\,$m$\,$s$^{-1}$"
    ax.text(0+0.05*ix,0+0.05*ix,label, color="yellow", fontsize=7, transform=ax.transAxes)

#ax.set_ylim(-10, 500)
ax.set_ylim(-10,600)
ax.legend(loc="upper left", ncol=3, columnspacing=1)
#ax.annotate('gravity current\nboundary', xy=(-48.8, 130), xytext=(-48.5, 270), #fontsize=9,
#            arrowprops = dict(facecolor='black', width = 2, shrink=0.05), ha = "center", va = "center", color = "white", bbox=dict(facecolor='black', alpha = 0.8, edgecolor='black', boxstyle='round, pad = 0.5'))

ax.set_facecolor('lightgrey')
ax.set_ylabel("Meters above bottom")
ax.set_xlabel(r"Longitude (°$\,$W)")
ax.invert_xaxis()
#ax[0].set_title(r"Dissipation rate $\varepsilon$ across the slope")

#fig.tight_layout()
fig.savefig("./eps_IGW_comparison.pdf")
fig.savefig("./eps_IGW_comparison.svg")
plt.show()