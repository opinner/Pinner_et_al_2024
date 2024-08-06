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

# depth-level-wise (row-wise) arithmetic averaging
rows = []
for index, row in thorpe_gamma_n_df.iterrows():
    values = row.to_numpy()
    bin_means = ss.binned_statistic(x=thorpe_lons, values=values, statistic=np.nanmean, bins=BIN_EDGES)[0]
    new_row = pd.DataFrame([bin_means], columns=BIN_EDGES[:-1])
    rows.append(new_row)

binned_thorpe_gamma_n_df = pd.concat(rows, sort=False).reset_index(drop=True)


#-------------------------------------------------------------------
# read eps_IGW results from strain-based finestructure analysis
eps_IGW_strain_df = pd.read_csv("../scripts/shear_strain_parametrization/method_results/binned_strain_eps.csv")
eps_IGW_strain_df.set_index('Unnamed: 0', inplace=True)

#-------------------------------------------------------------------
# read eps_IGW results from IDEMIX method
eps_IGW_IDEMIX_df = pd.read_csv("../scripts/IDEMIX_parametrization/method_results/eps_IGW_IDEMIX_results.csv")



fig,ax = plt.subplots(1, figsize=(TWO_COLUMN_WIDTH*cm, 0.8*TWO_COLUMN_WIDTH*cm))


######################################################################## 
######################################################################## 
##################### Axis 0 ###########################################
########################################################################
norm = mcolors.LogNorm(vmin=1e-10, vmax=1e-8)
mpp = ax.pcolormesh(eps_IGW_strain_df.columns.astype(float), eps_IGW_strain_df.index, eps_IGW_strain_df,
                      cmap=cmap,
                      norm=norm,
                      shading="nearest"
                     )
      
cb = plt.colorbar(mpp, ax=ax, location="top", extend="max")
cb.set_label(r"Wave-induced dissipation rate $\varepsilon_{\mathrm{IGW}}\,$(W$\,$kg$^{-1}$)")

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
    colors="white"
)


ax.scatter(
    eps_IGW_IDEMIX_df["lon"],
    eps_IGW_IDEMIX_df["rounded mab"],
    c=eps_IGW_IDEMIX_df["eps_IGW"],
    cmap=cmap,
    norm=norm,
    edgecolor='darkgrey',
    marker=MarkerStyle("o"),
    s=300,
    zorder=10
)

ax.set_facecolor('lightgrey')
ax.set_ylabel("Meters above bottom")
ax.set_xlabel("Longitude (°)")
#ax[0].set_title(r"Dissipation rate $\varepsilon$ across the slope")


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
    s = 200,
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
    s=80,
    zorder=-10,
    label=r"$\varepsilon_{\mathrm{IGW, fine}}$"
)


ax.set_ylim(0,600)
ax.legend(loc = "upper left", ncol=3, columnspacing=1)
#ax.annotate('gravity current\nboundary', xy=(-48.8, 130), xytext=(-48.5, 270), #fontsize=9,
#            arrowprops = dict(facecolor='black', width = 2, shrink=0.05), ha = "center", va = "center", color = "white", bbox=dict(facecolor='black', alpha = 0.8, edgecolor='black', boxstyle='round, pad = 0.5'))

            
fig.tight_layout()
fig.savefig("./eps_IGW_comparison.pdf")
# fig.savefig("./eps_IGW_comparison.png", dpi = 400, bbox_inches = "tight")
plt.show()