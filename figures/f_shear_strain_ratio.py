import cmocean
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
plt.style.use('./thesis.mplstyle')
import mixsea as mx
import numpy as np
import pandas as pd

import src.read_CTDs

ONE_COLUMN_WIDTH = 8.3
TWO_COLUMN_WIDTH = 12
GOLDEN_RATIO = 1.61
cm = 1/2.54  # centimeters in inches



def get_transect_CTDs():
    data = src.read_CTDs.get_PS129_CTD_data()

    def load_stations(path):
        transect_names, transect_numbers = np.loadtxt(path, dtype=(str), delimiter="\t", unpack=True)
        translate = dict(zip(transect_names, transect_numbers))
        # print(translate)
        return translate

    path = "/media/sf_VM_Folder/figures/PS129_Plots/Weddell_Sea_Transect.txt"
    translate = load_stations(path)
    transect_names = translate.keys()
    transect_numbers = translate.values()

    # Filter the DataFrame
    return data[data['Event'].str.replace('PS129_', '').isin(transect_names)]


transect_df = get_transect_CTDs()

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

# ----------------------------
# Finestructure params
# ----------------------------
shst_params = dict()

# Center points of depth windows. Windows are half overlapping, i.e.
# their size (300m) is double the spacing here (150m).
window_size = 250.0
min_size = 10.0
dz = window_size / 2
shst_params["depth_bin"] = np.arange(dz, 10000.0, dz)
shst_params["window_size"] = window_size

# Set up wavenumber vector.
shst_params["m"] = np.arange(
    2 * np.pi / window_size, 2 * np.pi / min_size, 2 * np.pi / window_size
)

# Set up limits for shear and strain variance integrations
mi_sh = np.array([0, 8])
mii_sh = np.array(range(*mi_sh))
mi_st = np.array([2, 20])
mii_st = np.array(range(*mi_st))

shst_params["m_include_sh"] = mii_sh
shst_params["m_include_st"] = mii_st

# Convert indices to more intuitive length scales
m_sh = 2 * np.pi / shst_params["m"][[mi_sh[0], mi_sh[1] - 1]]
m_st = 2 * np.pi / shst_params["m"][[mi_st[0], mi_st[1] - 1]]
print(
    f"Wavenumber indices for integration:\n"
    f"- Shear is integrated from {round(m_sh[0])}m to {round(m_sh[1])}m scales.\n"
    f"- Strain is integrated from {round(m_st[0])}m to {round(m_st[1])}m."
)

shst_params["ladcp_is_shear"] = True
shst_params["return_diagnostics"] = True

eps_list = []
Rw_list = []

for ctd, ladcp in zip(list_of_CTD_casts, list_of_LADCP_casts):

    depth = ctd["depth"]
    t = ctd["t"]
    SP = ctd["SP"]
    lon = ctd.location.lon
    lat = ctd.location.lat

    u = ladcp["u"]
    v = ladcp["v"]
    uz = ladcp["uz"]
    vz = ladcp["vz"]
    depth_sh = ladcp["depth"]

    try:
        eps, _krho, diag = mx.shearstrain.shearstrain(
            depth, t, SP, lon, lat, uz, vz, depth_sh, **shst_params
        )
    except ValueError:
        print(f"errors at {ctd.name}")
        continue

    depth_bin = diag["depth_bin"]
    depth_bin_edges = np.concatenate(
        (
            [np.min(depth_bin) - dz / 2],
            0.5 * (depth_bin[1:] + depth_bin[:-1]),
            [np.max(depth_bin) + dz / 2],
        )
    )
    eps_list.append(pd.DataFrame(index=depth_bin, data={lon: eps}))
    Rw_list.append(pd.DataFrame(index=depth_bin, data={lon: diag["Rwcor"]}))

Rw_df = pd.concat(Rw_list, axis=1)
Rw_df.sort_index(axis=1, inplace=True)  # sort columns
Rw_df.sort_index(inplace=True)  # sort rows TODO Why is this necessary?

fine_eps_df = pd.concat(eps_list, axis=1)
fine_eps_df.sort_index(axis=1, inplace=True)  # sort columns
fine_eps_df.sort_index(inplace=True)  # sort rows TODO Why is this necessary?

f, ax = plt.subplots(nrows=1, figsize=(10, 5))
# mab_bin_edges = bin_edges(eps_strain_df.index,dz)
# lon_edges = eps_strain_df.columns - np.diff(eps_strain_df.columns)
mpp = ax.pcolormesh(fine_eps_df.columns, fine_eps_df.index, fine_eps_df,
                    norm=mcolors.LogNorm(vmax=1e-8, vmin=1e-10),
                    shading="nearest",
                    cmap=cmocean.cm.matter
                    )
cb = plt.colorbar(mpp, ax=ax)
cb.set_label(r"$\varepsilon$ / (W kg$^{-1}$)")
ax.invert_yaxis()
ax.set_title(r"Turbulence diagnosed from Finestructure")
f.tight_layout()

f, ax = plt.subplots(nrows=1, figsize=(10, 5))
# mab_bin_edges = bin_edges(eps_strain_df.index,dz)
# lon_edges = eps_strain_df.columns - np.diff(eps_strain_df.columns)
mpp = ax.pcolormesh(Rw_df.columns, Rw_df.index, Rw_df,
                    norm=mcolors.LogNorm(vmax=300, vmin=1),
                    shading="nearest",
                    )
cb = plt.colorbar(mpp, ax=ax)
cb.set_label(r"Shear/Strain ratio Rw")
ax.invert_yaxis()
ax.set_title(r"Shear/Strain ratio Rw")
f.tight_layout()

print(f"Across the Weddell Sea: {np.nanmean(Rw_df.values.flatten()):.1f}±{np.nanstd(Rw_df.values.flatten()):.1f}")
fig, ax = plt.subplots(2, sharex=True, figsize=(0.8 * TWO_COLUMN_WIDTH * cm, 0.6 * TWO_COLUMN_WIDTH * cm), layout="constrained")
ax[0].hist(Rw_df.values.flatten(), bins=np.logspace(np.log10(1), np.log10(300), 20),
           color="k", edgecolor='white', linewidth=1)
ax[0].set_xscale("log")
ax[0].set_title(r"$R_\omega$, across the Weddell Sea transect")
ax[0].axvline(3, color="tab:red", ls="--")
ax[0].axvline(7, color="tab:red")

# Drop all rows that are outside the continental slope
min_longitude = -54
max_longitude = -47
# Filter columns based on the longitude range
filtered_columns = [col for col in Rw_df.columns if min_longitude <= float(col) <= max_longitude]

# Create a new DataFrame with the filtered columns
slope_Rw_df = Rw_df[filtered_columns]
print(
    f"Across continental slope: {np.nanmean(slope_Rw_df.values.flatten()):.1f}±{np.nanstd(slope_Rw_df.values.flatten()):.1f}")
ax[1].hist(slope_Rw_df.values.flatten(), bins=np.logspace(np.log10(1), np.log10(300), 20),
           color="k", edgecolor='white', linewidth=1)
ax[1].set_title(r"$R_\omega$, across the continental slope")
ax[1].axvline(3, color="tab:red", ls="--")
ax[1].axvline(7, color="tab:red")
fig.savefig(f"./Rw_histogram.pdf")
plt.show()
