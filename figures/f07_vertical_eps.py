import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

plt.rcParams.update({
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "font.size": 9,
    "svg.fonttype":'path'  # TrueType font
})

# define regions
regions_key = {"air": -1, "shelf": 0, "open ocean": 1, "IL": 2, "BL": 3}
binned_regions = pd.read_csv("../derived_data/gravity_current_regions.csv", index_col=0)
# convert column names from strings to floats
binned_regions.columns = binned_regions.columns.astype("float")
binned_regions = binned_regions.iloc[:600]

binned_neutral_density = pd.read_csv("../derived_data/binned_neutral_density.csv", index_col=0)
binned_neutral_density.columns = binned_neutral_density.columns.astype("float")
binned_neutral_density = binned_neutral_density.iloc[:600]

binned_thorpe_dissipation = pd.read_csv("../derived_data/binned_thorpe_dissipation.csv", index_col=0)
binned_thorpe_dissipation.columns = binned_thorpe_dissipation.columns.astype("float")
binned_thorpe_dissipation = binned_thorpe_dissipation.iloc[0:600]

binned_finestructure_dissipation = pd.read_csv("../derived_data/binned_finestructure_dissipation.csv", index_col=0)
binned_finestructure_dissipation.columns = binned_finestructure_dissipation.columns.astype("float")
binned_finestructure_dissipation = binned_finestructure_dissipation.iloc[:6]
wave_energy_dissipation = pd.read_csv("../derived_data/wave_energy_dissipation.csv", index_col=0)

mooring_lons = wave_energy_dissipation["lon"].unique()
desired_moorings = 0, 1, 4  # moorings A, B, E
desired_lons = [mooring_lons[i] for i in desired_moorings]
closest_lons = []
for lon in desired_lons:
    closest_index = np.argmin(np.abs(np.array(binned_neutral_density.columns) - lon))
    closest_lon = binned_neutral_density.columns[closest_index]
    print(f"mooring at {lon:.2f}, bin at {closest_lon:.2f}")
    closest_lons.append(closest_lon)

ONE_COLUMN_WIDTH = 8.3
TWO_COLUMN_WIDTH = 12
GOLDEN_RATIO = 1.61
cm = 1 / 2.54  # centimeters in inches


def plot_finestructure(ax, eps, mab, label, BL=None):
    for x, z in zip(eps, mab):
        spacing = np.abs(np.mean(np.diff(mab)))
        assert spacing == 125
        upper = z + spacing / 2
        lower = z - spacing / 2
        ax.vlines(x, lower, upper, lw=3, colors="tab:blue", label=label)
        #label = None  # use only the label from the first plot instance
        ax.fill_betweenx(
            [lower, upper], x / 5, x * 5,
            color="tab:blue", alpha=0.3,
            edgecolor=None)
        # draw in hatched BL
        if lower < 0 and BL is not None:
            ax.fill_betweenx(
                [lower, BL], x / 5, x * 5,
                color="none", edgecolor="xkcd:charcoal", alpha=0.8,
                hatch="xx"
            )


fig, axis = plt.subplots(
    ncols=3,
    figsize=(TWO_COLUMN_WIDTH * cm, TWO_COLUMN_WIDTH / GOLDEN_RATIO * cm),
    layout="constrained",
    sharex=True,
    sharey=True
)

xmin, xmax = (27.75, 28.7)
ymin, ymax = (-10, 500)
axis2 = []
for ax, mooring_lon, closest_lon in zip(axis, desired_lons, closest_lons):
    density = binned_neutral_density[closest_lon]
    regions = binned_regions[closest_lon]
    layers = np.squeeze(regions.diff().to_numpy().nonzero())

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.tick_params(axis='x', which='major', labelsize=6)
    ax.tick_params(axis='x', which='minor', labelsize=6)
    ax.plot(density.values, density.index,
            label="neutral density",
            color="k", ls="--", alpha=0.5, lw=1.5)

    #ax.plot(layers, layers, "s")
    layer_labels = ("Seafloor", "BL", "IL", "open ocean")
    if np.any(layers) != 0:
        for i, (layer, label) in enumerate(zip(layers, layer_labels[:-1])):
            if layer <= 0: continue
            ax.axhline(layer, xmin=0.85, xmax=1, color="k")
            ax.text(x=xmax-0.02, y=layer-16, s=label, fontsize=7, color="k", ha="right", va="top")
        ax.text(x=xmax-0.02, y=layer+16, s="open ocean", fontsize=7, rotation="vertical", ha="right", va="bottom", color="k")

    ax2 = ax.twiny()
    axis2.append(ax2)
    ax.xaxis.tick_top()
    ax2.xaxis.tick_bottom()
    ax2.set_xscale("log")
    ax2.set_xlim(7e-11, 2e-6)

    label = "$\\langle\\varepsilon_{\\mathrm{total, Thorpe}}\\rangle$"
    thorpe_eps = binned_thorpe_dissipation[closest_lon]
    ax2.plot(thorpe_eps, thorpe_eps.index, c="tab:red", label=label)
    ax2.fill_betweenx(
        thorpe_eps.index,
        thorpe_eps / 5, thorpe_eps * 5,
        color="tab:red", edgecolor=None, alpha=0.5,
    )

    try:
        BL = layers[1]
    except IndexError:
        BL = None

    label = "$\\langle\\varepsilon_{\\mathrm{IGW, fine}}\\rangle$"
    plot_finestructure(ax2,
                       binned_finestructure_dissipation[closest_lon],
                       binned_finestructure_dissipation.index,
                       label=label,
                       BL=BL)
    # ax.semilogx()
    # ax.semilogx(IGWs, mabs, "o", c="k", label="$\\langle\\varepsilon_{\\mathrm{IGW, IDEMIX}}\\rangle$", zorder=10)

    current_mooring = wave_energy_dissipation[wave_energy_dissipation["lon"] == mooring_lon]

    ax2.semilogx(current_mooring["eps_IGW"],
                 current_mooring["rounded mab"],
                 "o", c="k",
                 label="$\\langle\\varepsilon_{\\mathrm{IGW, IDEMIX}}\\rangle$",
                 zorder=10
                 )
    # Plot multiplicative errors
    for IGW, IGW_error, mab in zip(current_mooring["eps_IGW"], current_mooring["eps_IGW_mult_error"],
                                   current_mooring["rounded mab"]):
        ax2.plot([IGW / 5, IGW * 5], [mab, mab], lw=3, c="xkcd:dark grey", alpha=0.6)
        ax2.plot([IGW / IGW_error, IGW * IGW_error], [mab, mab], lw=3, c="k", alpha=1)


# from https://stackoverflow.com/questions/73915456/how-to-remove-duplicate-labels-in-subplot-legend
lines_labels = [ax2.get_legend_handles_labels() for ax2 in axis2]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
# grab unique labels
unique_labels = set(labels)
# assign labels and legends in dict
legend_dict = dict(zip(labels, lines))
# query dict based on unique labels
unique_lines = [legend_dict[x] for x in unique_labels]
new_symbols = []
for line, label in zip(unique_lines, unique_labels):
    if "Thorpe" in label:
        error_symbol = mpatches.Patch(color="tab:red", edgecolor=None, alpha=0.5, label='Thorpe_error')
        new_symbols.append((error_symbol, line))
    if "IDEMIX" in label:
        error_symbol = mlines.Line2D([], [], color='k', lw=3, label='IDEMIX_error')
        new_symbols.append((error_symbol, line))
    if "fine" in label:
        error_symbol = mpatches.Patch(color="tab:blue", alpha=0.3, edgecolor=None, label='fine_error')
        new_symbols.append((error_symbol, line))

fig.legend(new_symbols, unique_labels, ncols=3, fontsize="7", columnspacing=1, loc="upper center")

# neutral density legend, used in postprocessing
axis[0].legend(fontsize="7")

#ax.set_ylim(-10, 450)
#fig.legend(loc="upper left", framealpha=0.6, labelspacing=0.4, ncols=1, fontsize="9")
#fig.legend(loc="upper left", framealpha=0.6, labelspacing=0.4, ncols=3, fontsize="9")
#ax.set_xlim(3e-11,5e-8)

axis[0].text(x=xmin + 0.04, y=200, s="continental shelf", fontsize=7, rotation="vertical", ha="left", va="center",
        color="k")
fig.supylabel("Height above bottom (m)", fontsize="9")

fig.supxlabel("Dissipation rate $\\varepsilon\\,$(W kg$^{-1})$", fontsize="9")
for ax in axis:
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
axis[1].set_xlabel(r"Neutral Density $\gamma^\text{n}\,$(kg$\,$m$^{-3}$)", fontsize="6")

# make space for legend placement in postprocessing
legend_axis = axis[-1].twinx()
legend_axis.tick_params(
    axis='y',          # changes apply to the y-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelright=False
)
legend_axis.set_ylabel("test space", fontsize=18)

moorings = ["A", "B", "E"]
for ax, mooring in zip(axis, moorings):
    ax.text(
        0.88, 0.91,
        mooring,
        transform=ax.transAxes,
        fontsize="9",
        weight="bold"
    )


fig.savefig("./vertical_eps.svg")
#fig.savefig("./vertical_eps.pdf")
plt.show()
