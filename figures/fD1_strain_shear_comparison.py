import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cmocean
import mixsea as mx
import scipy.stats as ss
from pandas import DataFrame
import warnings

warnings.filterwarnings('ignore', category=RuntimeWarning)
import src.read_CTDs

plt.rcParams.update({
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "font.size": 9
})

ONE_COLUMN_WIDTH = 8.3
TWO_COLUMN_WIDTH = 12
GOLDEN_RATIO = 1.61
cm = 1 / 2.54  # centimeters in inches


def get_transect_CTDs():
    data = src.read_CTDs.get_PS129_CTD_data()

    def load_stations(path):
        transect_names, transect_numbers = np.loadtxt(path, dtype=str, delimiter="\t", unpack=True)
        translate = dict(zip(transect_names, transect_numbers))
        # print(translate)
        return translate

    path = "/media/sf_VM_Folder/figures/PS129_Plots/Weddell_Sea_Transect.txt"
    translate = load_stations(path)
    transect_names = translate.keys()
    transect_numbers = translate.values()

    print(data['Event'])

    # Filter the DataFrame
    transect_df = data[data['Event'].str.replace('PS129_', '').isin(transect_names)]
    return transect_df


transect_df = get_transect_CTDs()
assert not transect_df.empty
# load LADCP data
list_of_LADCP_casts, list_of_CTD_casts = src.read_CTDs.get_PS129_CTDs_and_LADCPs()
# only use CTDS from the Weddell Sea transect
list_of_LADCP_casts[:] = [cast for cast in list_of_LADCP_casts if
                          cast.location.lon in transect_df['Longitude'].to_numpy()]
list_of_CTD_casts[:] = [cast for cast in list_of_CTD_casts if cast.location.lon in transect_df['Longitude'].to_numpy()]

# check order
for LADCP_cast, CTD_cast in zip(list_of_LADCP_casts, list_of_CTD_casts):
    if LADCP_cast.name != CTD_cast.name:
        print(LADCP_cast.name, CTD_cast.name)
        print("Wrong order")
        raise AssertionError


def create_fixed_step_array_including_seafloor(start, stop, step, fixed_depth):
    """
    Create a fixed step array, that includes a given fixed depth, in this case the sea floor
    For that the start point is adjusted
    """
    # Adjust the start point to be a multiple of step such that specific_value can be included
    start_adjusted = start + (fixed_depth - start) % step

    # Generate the regularly spaced array starting from the adjusted start
    array = np.arange(start_adjusted, stop + step, step)

    # Ensure the specific value is in the array
    if fixed_depth not in array:
        raise ValueError(
            f"The specific value {fixed_depth} cannot be included with the given start, stop, and step values.")
    return array


# ----------------------------
# Finestructure params
# ----------------------------
shst_params = dict()
# Center points of depth windows. Windows are half overlapping, i.e.
# their size (300m) is double the spacing here (150m).
window_size = 250.0
min_size = 10.0
dz = window_size / 2
shst_params["window_size"] = window_size

# Set up wavenumber vector.
shst_params["m"] = np.arange(
    2 * np.pi / window_size, 2 * np.pi / min_size, 2 * np.pi / window_size
)

# Set up limits for shear and strain variance integrations
mi_sh = np.array([0, 8])
mii_sh = np.array(range(*mi_sh))
mi_st = np.array([0, 20])
mii_st = np.array(range(*mi_st))
shst_params["m_include_sh"] = mii_sh
shst_params["m_include_st"] = mii_st
# Convert indices to more intuitive length scales
m_sh = 2 * np.pi / shst_params["m"][[mi_sh[0], mi_sh[1] - 1]]
m_st = 2 * np.pi / shst_params["m"][[mi_st[0], mi_st[1] - 1]]
print(
    f"Wavenumber indices for integration:\n"
    f"- Shear is integrated from {np.round(m_sh[0])}m to {np.round(m_sh[1])}m scales.\n"
    f"- Strain is integrated from {np.round(m_st[0])}m to {np.round(m_st[1])}m."
)
shst_params["ladcp_is_shear"] = True
shst_params["return_diagnostics"] = True

print("Finestructure method is ongoing...")
eps_list = []
eps_strain_list = []
for ctd, ladcp in zip(list_of_CTD_casts, list_of_LADCP_casts):

    depth = ctd["depth"]
    lowest_segment: int = np.floor(depth.max() - dz / 2)
    t = ctd["t"]
    SP = ctd["SP"]
    lon = ctd.location.lon
    lat = ctd.location.lat

    uz = ladcp["uz"]
    vz = ladcp["vz"]
    depth_sh = ladcp["depth"]

    if lowest_segment < dz:
        print(f"profile with {lowest_segment = }m depth, at {lon}, is too shallow")
        continue
    depth_bins = create_fixed_step_array_including_seafloor(start=dz, stop=10000.0, step=dz, fixed_depth=lowest_segment)
    shst_params["depth_bin"] = depth_bins
    try:
        eps, _krho, diag = mx.shearstrain.shearstrain(
            depth, t, SP, lon, lat, uz, vz, depth_sh, **shst_params
        )
    except ValueError:
        print(f"errors at {ctd.name}")
        continue

    depth_bins = diag["depth_bin"]
    # use meters above bottom as y-axis
    mab_bins = np.floor(depth.max()) - depth_bins
    if mab_bins[-1] % 2 != 0:
        mab_bins = mab_bins - 1
    eps_list.append(pd.DataFrame(index=mab_bins, data={lon: eps}))
    eps_strain_list.append(pd.DataFrame(index=mab_bins, data={lon: diag["eps_st"]}))

eps_df: DataFrame = pd.concat(eps_list, axis=1)
eps_df.sort_index(axis=1, inplace=True)  # sort columns
eps_df.sort_index(inplace=True)  # sort rows TODO Why is this necessary?
#eps_df.columns = [el[0] for el in eps_df.columns]  # convert multiindex to single index

eps_strain_df = pd.concat(eps_strain_list, axis=1)
eps_strain_df.sort_index(axis=1, inplace=True)  # sort columns
eps_strain_df.sort_index(inplace=True)  # sort rows TODO Why is this necessary?
#eps_strain_df.columns = [el[0] for el in eps_strain_df.columns]  # convert multiindex to single index

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
eps_strain_df = eps_strain_df * 2.694  #Correction from Rw =3 to Rw = 7
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

comparison_df = eps_strain_df / eps_df

thorpe_gamma_n_df = pd.read_pickle("../scripts/thorpe_scales/method_results/Thorpe_neutral_density_df_with_mab.pkl")
thorpe_mab = thorpe_gamma_n_df.index

thorpe_lons = thorpe_gamma_n_df.columns.to_numpy()
max_lon = max(thorpe_lons)
min_lon = min(thorpe_lons)

# half a degree bins
#BIN_EDGES = np.arange(min_lon - 1e-3 * min_lon, 0.5 + max_lon + 1e-3 * max_lon, 0.5)
BIN_EDGES = np.arange(-54.25, -46.25, 0.5)
BIN_CENTER = BIN_EDGES[:-1] + 0.25

# depth-level-wise (row-wise) arithmetic averaging
rows = []
for index, row in thorpe_gamma_n_df.iterrows():
    values = row.to_numpy()
    bin_means = ss.binned_statistic(x=thorpe_lons, values=values, statistic=np.nanmean, bins=BIN_EDGES)[0]
    new_row = pd.DataFrame([bin_means], columns=BIN_CENTER)
    rows.append(new_row)

binned_thorpe_gamma_n_df = pd.concat(rows, sort=False).reset_index(drop=True)

fig, ax = plt.subplots(1,
                       figsize=(TWO_COLUMN_WIDTH * cm, 0.5 * TWO_COLUMN_WIDTH * cm),
                       layout="constrained")
ax.set_ylim(0, 1000)
ax.set_xlim(-53.2, -47)
factor = 20

mapp = ax.pcolormesh(
    comparison_df.columns,
    comparison_df.index,
    comparison_df,
    cmap=cmocean.cm.balance,
    norm=mcolors.LogNorm(vmin=1/factor, vmax=1*factor)
)

cb = plt.colorbar(mapp, ax=ax, aspect=15, extend="both")
cb.set_label(r"$\varepsilon_\mathrm{IGW, strain}$ / $\varepsilon_\mathrm{IGW, shear}$")
ax.set_facecolor('lightgrey')

# continental_slope = BIN_EDGES[3]
# deep_sea = BIN_EDGES[-1]
# ax.fill_between(
#     [continental_slope, deep_sea],
#     [0,0],
#     [250/2, 250/2],
#     hatch="xx",
#     facecolor='None',
#     edgecolor='darkgrey',
#     alpha=0.8
# )

water_mass_boundaries = [28.26, 28.40]  # + 28.00 boundary, from Garabato et al 2002

plot_gamma_n_df = binned_thorpe_gamma_n_df[binned_thorpe_gamma_n_df.index < 1050]
CS = ax.contour(
    plot_gamma_n_df.columns,
    plot_gamma_n_df.index,
    plot_gamma_n_df,
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
clabels = ax.clabel(
    CS,
    CS.levels,
    inline=False,
    fmt=fmt,
    fontsize=10,
    colors="black"
)
# adjust bboxes for better readability
[txt.set_bbox(dict(facecolor='lightgrey', alpha=0.8, edgecolor='darkgrey', boxstyle="round", pad=0)) for txt in clabels]

binned_regions = pd.read_csv(
    "../scripts/preprocessing/method_results/binned_regions.csv", index_col=0)
binned_regions.columns = binned_regions.columns.astype("float")  #convert column names from strings to floats
binned_regions = binned_regions.iloc[0:600]
levels = [2.5, 3.5]  # Border between IL and BL
plt.rcParams['hatch.color'] = 'xkcd:charcoal'
ax.contourf(
    binned_regions.columns,
    binned_regions.index,
    binned_regions.values,
    levels=levels,
    hatches=["x"],
    colors="None",
    zorder=10
)
ax.contour(
    binned_regions.columns,
    binned_regions.index,
    binned_regions.values,
    levels=levels,
    colors="xkcd:charcoal",
    zorder=10
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
    linewidths=1.5,
    alpha=0.8
)

ax.set_facecolor('lightgrey')
ax.set_ylabel("Meters above bottom")
ax.set_xlabel("Longitude (°)")

#fig.tight_layout()
fig.savefig("./strain_shear_comparison.pdf")

#print(f"Largest deviations: {comparison_df.max(axis=None):.2e}, {comparison_df.min(axis=None):.2e}")
print("done")
plt.show()
