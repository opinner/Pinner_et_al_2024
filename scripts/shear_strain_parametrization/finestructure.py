import warnings

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import mixsea as mx
import numpy as np
import pandas as pd
# data handling
import scipy.stats as ss

warnings.filterwarnings('ignore', category=RuntimeWarning)

# import my self written functions
from src.read_CTDs import load_Joinville_transect_CTDs

CTDs = load_Joinville_transect_CTDs()
CTDs_grouped = CTDs.groupby("Expedition")
expedition_names = CTDs_grouped.groups.keys()
print(expedition_names)


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
mi_sh = np.array([0, 3])
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


def create_fixed_step_array_includ_seafloor(start, stop, step, fixed_depth):
    # Adjust the start point to be a multiple of step such that specific_value can be included
    start_adjusted = start + (fixed_depth - start) % step

    # Generate the regularly spaced array starting from the adjusted start
    array = np.arange(start_adjusted, stop + step, step)

    # Ensure the specific value is in the array
    if fixed_depth not in array:
        raise ValueError(
            f"The specific value {fixed_depth} cannot be included with the given start, stop, and step values.")

    return array


# empty objects for saving the data later
eps_list = []
eps_strain_list = []

for i, expedition_name in enumerate(expedition_names):
    # print(expedition_name)

    expedition = CTDs_grouped.get_group(expedition_name).reset_index(drop=True)
    expedition = expedition.groupby("Event")
    events = expedition.groups.keys()

    for event in events:
        if expedition_name not in event: continue

        current_profile = expedition.get_group(event).reset_index(drop=True)

        depth = current_profile['Depth water [m]']
        lowest_segment: int = np.floor(depth.max() - dz/2)
        t = current_profile['Temp [°C]']
        SP = current_profile['Sal']
        lon = current_profile["Longitude"].mean(),
        lat = current_profile["Latitude"].mean()

        if lowest_segment < dz:
            print(f"profile with {lowest_segment = }m depth, at {lon}, is too shallow")
            continue
        depth_bins = create_fixed_step_array_includ_seafloor(start=dz, stop=10000.0, step=dz, fixed_depth = lowest_segment)
        shst_params["depth_bin"] = depth_bins
        try:
            eps, krho, diag = mx.shearstrain.nan_shearstrain(
                depth, t, SP, lon, lat, **shst_params
            )
        except ValueError:
            print(f"errors at {expedition_name}, {event}")
            continue

        depth_bins = diag["depth_bin"]
        # use meters above bottom as y axis
        mab_bins = np.floor(depth.max()) - depth_bins
        if mab_bins[-1] % 2 != 0:
            mab_bins = mab_bins-1
        eps_list.append(pd.DataFrame(index=mab_bins, data={lon: eps}))
        eps_strain_list.append(pd.DataFrame(index=mab_bins, data={lon: diag["eps_st"]}))

        # use depth as y axis
        # eps_list.append(pd.DataFrame(index =  depth_bin, data = {lon: eps}))
        # eps_strain_list.append(pd.DataFrame(index = depth_bin, data = {lon: diag["eps_st"]}))

eps_strain_df = pd.concat(eps_strain_list, axis=1)
eps_strain_df.sort_index(axis=1, inplace=True)  # sort columns
eps_strain_df.sort_index(inplace=True)  # sort rows TODO Why is this necessary?
eps_strain_df.columns = [el[0] for el in eps_strain_df.columns]  # convert multiindex to single index

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
eps_strain_df = eps_strain_df * 2.694  #Correction from Rw =3 to Rw = 7
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

# def fexp(f):
#     return int(np.floor(np.log10(abs(f)))) if f != 0 else 0
# def fman(f):
#     return f / 10 ** fexp(f)
# def generate_exp_pattern(start, end):
#     num_steps = (end - start) * 2  # Each step has a multiplier of 1 and 5
#     return [multiplier * 10 ** (start + i // 2) for i, multiplier in enumerate([1, 5] * (num_steps + 1)) if
#             start + i // 2 <= end]
#
#
# start_point = -12
# end_point = -8
# bounds = generate_exp_pattern(start_point, end_point)[:-1]
# print(bounds)
# ncolors = len(bounds) - 1
# cmap = plt.cm.get_cmap('magma_r', ncolors)
# norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=ncolors)
#
# f, ax = plt.subplots(nrows=1, figsize=(10, 5))
#
# # mab_bin_edges = bin_edges(eps_strain_df.index,dz)
# # lon_edges = eps_strain_df.columns - np.diff(eps_strain_df.columns)
# mpp = ax.pcolormesh(eps_strain_df.columns, eps_strain_df.index, eps_strain_df,
#                     cmap=cmap,
#                     norm=norm,
#                     shading="nearest"
#                     )
# cb = plt.colorbar(mpp, ax=ax)
# cb.set_label(r"$\varepsilon$ / (W kg$^{-1}$)")
# cb.ax.set_yticklabels([f'{fman(b):.1f}$\\times10^{{{fexp(b):.0f}}}$' for b in bounds])
# ax.set_facecolor('lightgrey')
# # ax.set_ylim(0,500)
#
# ax.set_title(r"Mixing diagnosed from Strain parametrization")
# # helper.Plot.path_as_footnote(fig = f,
# #                             path = "/home/ole/Desktop/CTD/mixsea/Weddell Sea Thorpe Scale.ipynb",
# #                             rot = "vertical")
# f.tight_layout()
# # f.savefig("./ThorpeDissipation_individ_Cruises.png", dpi = 300)

# trim to gravity current core
vertical_eps_df = eps_strain_df.drop(eps_strain_df.columns[eps_strain_df.columns < -51.5], axis="columns")
vertical_eps_df.drop(vertical_eps_df.columns[vertical_eps_df.columns > -48.5], axis="columns", inplace=True)

mean = vertical_eps_df.mean(axis=1)
std = vertical_eps_df.std(axis=1)
np.savez("method_results/Strain_vertical_eps", mab =vertical_eps_df.index, eps=mean)

# fig, ax = plt.subplots(1)

# """
# column_number = len(vertical_eps_df.columns)
# n = 5
# outliers = (vertical_eps_df > np.tile(mean+n*std, (column_number,1)).T) | (vertical_eps_df < np.tile(mean-n*std, (column_number,1)).T)
# # Set outliers to NaN
# _vertical_eps_df = vertical_eps_df.copy(deep = True)
# _vertical_eps_df[outliers] = pd.NA
# # Calculate new mean, ignoring NaNs
# #mean = _vertical_eps_df.mean(axis = 1)
# #std = _vertical_eps_df.std(axis = 1)
# """
# # ax.fill_betweenx(vertical_eps_df.index, mean-std, mean+std, alpha = 0.5 )
# # ax.semilogx(mean, _vertical_eps_df.index)
# ax.semilogx(mean, vertical_eps_df.index)
# # ax.set_ylim(-10,510)
# print(f"{vertical_eps_df.min(axis=None):.2e},{vertical_eps_df.max(axis=None):.2e}")


lons = eps_strain_df.columns.to_numpy()
max_lon = max(lons)
min_lon = min(lons)
LON_BIN_EDGES = np.arange(min_lon - 1e-3 * min_lon, 0.5+max_lon + 1e-3 * max_lon, 0.5)

rows = []
for index, row in eps_strain_df.iterrows():
    values = row.to_numpy()
    bin_means = ss.binned_statistic(x=lons, values=values, statistic=np.nanmean, bins=LON_BIN_EDGES)[0]
    new_eps = bin_means
    new_row = pd.DataFrame([new_eps], columns=LON_BIN_EDGES[:-1])
    rows.append(new_row)
binned_eps_strain_df = pd.concat(rows, sort=False).reset_index(drop=True)
binned_eps_strain_df.index = eps_strain_df.index

print(binned_eps_strain_df.head())
print(binned_eps_strain_df.info())
binned_eps_strain_df.to_csv("./method_data/binned_strain_eps.csv")

plt.show()
