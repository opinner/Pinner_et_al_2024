import pandas as pd
import matplotlib.pyplot as plt
import cmocean as cmocean
import warnings
# Suppress specific RuntimeWarning related to mean of empty slice
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*Mean of empty slice.*")

ONE_COLUMN_WIDTH = 8.3
TWO_COLUMN_WIDTH = 12
GOLDEN_RATIO = 1.61
cm = 1/2.54  # centimeters in inches

#plt.style.use('./paper.mplstyle')
plt.rcParams.update({
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "font.size": 9,
    "svg.fonttype":'path'  # TrueType font
})

# binned_neutral_density = pd.read_csv("../derived_data/binned_neutral_density.csv", index_col=0)
# binned_neutral_density.columns = binned_neutral_density.columns.astype("float") #convert column names from strings to floats
# binned_regions = pd.read_csv("../scripts/preprocessing/method_results/binned_regions.csv", index_col=0)
# binned_regions.columns = binned_regions.columns.astype("float") #convert column names from strings to floats
# binned_regions = binned_regions.iloc[0:600]

neutral_density = pd.read_csv("../scripts/preprocessing/method_results/PS129_transect.csv",
                               index_col = 0)
neutral_density.columns = neutral_density.columns.astype("float") #convert column names from strings to floats


fig, ax = plt.subplots(1, figsize=(TWO_COLUMN_WIDTH*cm, 0.5*TWO_COLUMN_WIDTH*cm), layout="constrained")
mpp = ax.pcolormesh(
    -neutral_density.columns,
    neutral_density.index,
    neutral_density.values,
    cmap=cmocean.cm.dense,
    rasterized=True,  # optimize the drawing for vector graphics
    vmin=27.8
)

#cb = plt.colorbar(mpp, ax=ax, extend="min", location="top")  # pad=0.02, aspect=12
cb = fig.colorbar(mpp, extend="min", ax=ax, aspect=8)

# AASW/WDW, WDW/WSDW and WSDW/WSBW interfaces, from Garabato et al 2002
water_mass_boundaries = [28.0, 28.26, 28.40]
CS = ax.contour(
    -neutral_density.columns,
    neutral_density.index,
    neutral_density.values,
    levels=water_mass_boundaries,
    linestyles=["dotted", "dashed", "solid"],
    colors="white",
    #linewidths=1,
)

# levels = [2.5, 3.5]  # Border between IL and BL
# ax.contour(
#     binned_regions.columns,
#     binned_regions.index,
#     binned_regions.values,
#     levels=levels,
#     colors="white",
# )

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
#     fontsize=10,
#     colors="white"
# )

# to be shifted in postprocessing
strs = ['WDW', 'WSDW', 'WSBW',]
for ix, s in enumerate(strs):
    ax.text(0.9, 0.9-0.05*ix,s, color="white", fontsize=8, transform=ax.transAxes)


# ax.annotate('bottom\ncurrent', xy=(-51.69, 184), xytext=(-51.8, 230),
#              arrowprops=dict(facecolor='black', width=2, shrink=0.05), ha="center", color="white",
#              bbox=dict(facecolor='black', alpha=0.8, edgecolor='black', boxstyle='round'))
# ax.annotate('gravity current\nboundary', xy=(-48.8, 130), xytext=(-48, 270), #fontsize=9,
#             arrowprops = dict(facecolor='black', width = 2, shrink=0.05), ha = "center", va = "center", color = "white", bbox=dict(facecolor='black', alpha = 0.8, edgecolor='black', boxstyle='round, pad = 0.5'))

mooring_info = pd.read_csv("../scripts/IDEMIX_parameterization/method_results/eps_IGW_IDEMIX_results.csv")
moorings_depths = mooring_info["rounded depth"]
moorings_lons = mooring_info["lon"]

# draw measurement positions
ax.plot(-moorings_lons, moorings_depths,
        ls="None",
        marker="o",
        label="moored\nrotor current meters",
        color="tab:red",
        markersize=6,
        markeredgecolor="k",
        zorder=6)

# ax.set_ylim(-10, 500)
# xlim = ax.get_xlim()
# ax.set_xlim((xlim[0] - 0.2, xlim[1] + 0.1))
ax.set_xlabel("Longitude (° W)")
ax.set_ylabel("Depth (m)")
cb.set_label(r"Neutral Density $\gamma^\text{n}\,$(kg$\,$m$^{-3}$)")
cb.ax.plot([0, 1], [water_mass_boundaries[0], water_mass_boundaries[0]], ':', color="white", lw=2)
cb.ax.plot([0, 1], [water_mass_boundaries[1], water_mass_boundaries[1]], '--', color="white", lw=2)
cb.ax.plot([0, 1], [water_mass_boundaries[2], water_mass_boundaries[2]], ls="solid", color="white",lw=2, )
cb.ax.invert_yaxis()

# draw sea floor
# x = list(ax.get_xlim())
# y1 = [0, 0]
# y2 = 2 * list(ax.get_ylim())[0]
# ax.axhline(0, c="k")
# ax.fill_between(x, y1, y2, facecolor="xkcd:charcoal grey", zorder=5)  # , hatch="///")
ax.legend(loc="lower left")  #,facecolor='k', framealpha=0.8, edgecolor = "black", labelcolor = "white")

#ax.grid()
ax.invert_yaxis()
ax.invert_xaxis()
ax.set_facecolor("lightgrey")
#fig.savefig("./stratification.pdf")
fig.savefig("./stratification.svg") #, dpi = 400, bbox_inches = "tight")
plt.show()
