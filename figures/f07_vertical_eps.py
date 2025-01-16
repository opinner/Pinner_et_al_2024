import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "font.size": 9
})

binned_neutral_density = pd.read_csv("../derived_data/binned_neutral_density.csv", index_col = 0)
binned_neutral_density.columns = binned_neutral_density.columns.astype("float") #convert column names from strings to floats
binned_neutral_density = binned_neutral_density.iloc[:600]

binned_thorpe_dissipation = pd.read_csv("../derived_data/binned_thorpe_dissipation.csv", index_col = 0)
binned_thorpe_dissipation.columns = binned_thorpe_dissipation.columns.astype("float") #convert column names from strings to floats
binned_thorpe_dissipation = binned_thorpe_dissipation.iloc[0:600]

binned_finestructure_dissipation = pd.read_csv("../derived_data/binned_finestructure_dissipation.csv", index_col = 0)
binned_finestructure_dissipation.columns = binned_finestructure_dissipation.columns.astype("float") #convert column names from strings to floats
binned_finestructure_dissipation = binned_finestructure_dissipation.iloc[:6]
wave_energy_dissipation = pd.read_csv("../derived_data/wave_energy_dissipation.csv", index_col = 0)

# define regions
regions_key = {"air": -1, "shelf": 0, "open ocean": 1, "IL": 2, "BL": 3}
binned_regions = pd.read_csv("../derived_data/gravity_current_regions.csv", index_col = 0)
binned_regions.columns = binned_regions.columns.astype("float")

mooring_lons = wave_energy_dissipation["lon"].unique()
desired_moorings = 0,1,4 # moorings A, B, E
desired_lons = [mooring_lons[i] for i in desired_moorings]
closest_lons = []
for lon in desired_lons:
    closest_index = np.argmin(np.abs(np.array(binned_neutral_density.columns) - lon))
    closest_lon = binned_neutral_density.columns[closest_index]
    print(lon, closest_lon)
    closest_lons.append(closest_lon)



ONE_COLUMN_WIDTH = 8.3
TWO_COLUMN_WIDTH = 12
GOLDEN_RATIO = 1.61
cm = 1/2.54  # centimeters in inches

def plot_finestructure(ax, eps, mab, label):
    for x, z in zip(eps, mab):
        spacing = np.abs(np.mean(np.diff(mab)))
        assert spacing == 125
        upper = z + spacing / 2
        lower = z - spacing / 2
        #if lower < 0: lower = 0
        ax.vlines(x, lower, upper, lw=3, colors="tab:blue", label=label)
        label = None  # use only the label from the first plot instance
        ax.fill_betweenx([lower, upper], x / 5, x * 5, color="tab:blue", alpha=0.3)

fig, axis = plt.subplots(
    ncols=3,
    figsize=(TWO_COLUMN_WIDTH * cm, TWO_COLUMN_WIDTH / GOLDEN_RATIO * cm),
    layout="constrained",
    sharex=True,
    sharey=True
)

axis2 = []
for ax, mooring_lon, closest_lon in zip(axis, desired_lons, closest_lons):
    profile = binned_neutral_density[closest_lon]

    ax.set_xlim(27.75, 28.7)
    ax.set_ylim(-10, 500)
    ax.tick_params(axis='x', which='major', labelsize=6)
    ax.tick_params(axis='x', which='minor', labelsize=6)
    ax.plot(profile.values, profile.index, "k", ls="--", alpha=0.5, lw=1.5)

    ax2 = ax.twiny()
    axis2.append(ax2)
    ax.xaxis.tick_top()
    ax2.xaxis.tick_bottom()
    ax2.set_xscale("log")
    ax2.set_xlim(7e-11, 2e-6)

    label = "$\\langle\\varepsilon_{\\mathrm{total, Thorpe}}\\rangle$"
    thorpe_eps = binned_thorpe_dissipation[closest_lon]
    ax2.plot(thorpe_eps, thorpe_eps.index, c="tab:red", label=label)
    ax2.fill_betweenx(thorpe_eps.index, thorpe_eps / 5, thorpe_eps * 5, color="tab:red", alpha=0.5)

    label = "$\\langle\\varepsilon_{\\mathrm{IGW, fine}}\\rangle$"
    plot_finestructure(ax2, binned_finestructure_dissipation[closest_lon], binned_finestructure_dissipation.index,
                       label=label)
    # ax.semilogx()
    # ax.semilogx(IGWs, mabs, "o", c="k", label="$\\langle\\varepsilon_{\\mathrm{IGW, IDEMIX}}\\rangle$", zorder=10)

    current_mooring = wave_energy_dissipation[wave_energy_dissipation["lon"] == mooring_lon]

    ax2.semilogx(current_mooring["eps_IGW"], current_mooring["rounded mab"], "o", c="k",
                 label="$\\langle\\varepsilon_{\\mathrm{IGW, IDEMIX}}\\rangle$", zorder=10)
    # Plot multiplicative errors
    for IGW, IGW_error, mab in zip(current_mooring["eps_IGW"], current_mooring["eps_IGW_mult_error"],
                                   current_mooring["rounded mab"]):
        ax2.plot([IGW / IGW_error, IGW * IGW_error], [mab, mab], lw=3, c="k", alpha=0.5)


# from https://stackoverflow.com/questions/73915456/how-to-remove-duplicate-labels-in-subplot-legend
lines_labels = [ax2.get_legend_handles_labels() for ax2 in axis2]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
# grab unique labels
unique_labels = set(labels)
# assign labels and legends in dict
legend_dict = dict(zip(labels, lines))
# query dict based on unique labels
unique_lines = [legend_dict[x] for x in unique_labels]
#fig.legend(unique_lines, unique_labels, fontsize="9", loc="outside right upper")

fig.supylabel("Height above bottom (m)", fontsize="9")

fig.supxlabel("Dissipation rate $\\varepsilon\\,$(W kg$^{-1})$", fontsize="9")
for ax in axis:
    ax.xaxis.tick_top()
    ax.set_xlabel("$\\gamma\\,$(kg m$^{-3})$", fontsize="6")
    ax.xaxis.set_label_position('top')

moorings = ["A","B","E"]
for ax, mooring in zip(axis,moorings):
    ax.text(
        0.88, 0.91,
        mooring,
        transform=ax.transAxes,
        fontsize="9",
        weight="bold"
    )

#ax.set_ylim(-10, 450)
#ax.legend(loc="upper left", framealpha=0.6, labelspacing=0.4, ncols=1, fontsize="9")
#ax.set_xlim(3e-11,5e-8)
#fig.savefig("./avrg_vertical_eps.png", dpi=400, bbox_inches="tight")
fig.savefig("./vertical_eps.pdf")
plt.show()
