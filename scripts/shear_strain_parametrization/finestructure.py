import matplotlib.pyplot as plt
import mixsea as mx
import numpy as np
import pandas as pd
# data handling
import scipy.stats as ss
import warnings
# import my self written functions
from src.read_CTDs import load_Joinville_transect_CTDs
warnings.filterwarnings('ignore', category=RuntimeWarning)

OUTLIERS = ['PS71/216-1', 'PS40/099-1', 'PS49/015-2', 'PS71/212-3', 'PS71/210-2']

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
#mi_sh = np.array([0, 8])
#mii_sh = np.array(range(*mi_sh))
mi_st = np.array([2, 20])
mii_st = np.array(range(*mi_st))
#shst_params["m_include_sh"] = mii_sh #not used here, but set for consistency with the shear-based finestructure computation
shst_params["m_include_st"] = mii_st
# Convert indices to more intuitive length scales
#m_sh = 2 * np.pi / shst_params["m"][[mi_sh[0], mi_sh[1] - 1]]
m_st = 2 * np.pi / shst_params["m"][[mi_st[0], mi_st[1] - 1]]
print(
    f"Wavenumber indices for integration:\n"
    f"- Strain is integrated from {round(m_st[0])}m to {round(m_st[1])}m.\n"
    #f"- Shear would be integrated from {round(m_sh[0])}m to {round(m_sh[1])}m scales, but not now"

)
#shst_params["ladcp_is_shear"] = True #not used here, but set for consistency with the shear-based finestructure computation
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


# empty list for saving the strain-based dissipation rates
eps_strain_list = []

for i, expedition_name in enumerate(expedition_names):
    # print(expedition_name)

    expedition = CTDs_grouped.get_group(expedition_name).reset_index(drop=True)
    expedition = expedition.groupby("Event")
    events = expedition.groups.keys()

    for event in events:
        if expedition_name not in event: continue
        if event in OUTLIERS:
            continue

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
            _eps, krho, diag = mx.shearstrain.nan_shearstrain(
                depth, t, SP, lon, lat, **shst_params
            )
            assert np.all(np.isnan(_eps)) # shear-based estimation should not be possible here
        except ValueError:
            print(f"errors at {expedition_name}, {event}")
            continue

        depth_bins = diag["depth_bin"]
        # use meters above bottom as y axis
        mab_bins = np.floor(depth.max()) - depth_bins
        # align mab bins to correct earlier inconsistent rounding up or down
        if mab_bins[-1] % 2 != 0:
            mab_bins = mab_bins-1
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


# trim to gravity current core
vertical_eps_df = eps_strain_df.drop(eps_strain_df.columns[eps_strain_df.columns < -51.5], axis="columns")
vertical_eps_df.drop(vertical_eps_df.columns[vertical_eps_df.columns > -48.5], axis="columns", inplace=True)

# take horizontal average to get a single vertical profile
mean = vertical_eps_df.mean(axis=1)
std = vertical_eps_df.std(axis=1)
np.savez("method_results/Strain_vertical_eps.npz", mab =vertical_eps_df.index, eps=mean)

# Bin resulting dissipation rates
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

print(binned_eps_strain_df.head(),"\n")

print(binned_eps_strain_df.info(),"\n")
binned_eps_strain_df.to_csv("./method_results/binned_strain_eps.csv")
print("Done")
