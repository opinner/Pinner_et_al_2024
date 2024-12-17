import numpy as np
import pandas as pd
import mixsea as mx
from src.read_CTDs import load_Joinville_transect_CTDs
from scipy.interpolate import interp1d  # is considered legacy code, will be in the future removed from scipy
import scipy.stats as ss
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
# Suppress specific RuntimeWarning related to mean of empty slice
warnings.filterwarnings(action="ignore", category=RuntimeWarning, message=".*Mean of empty slice.*")

CTDs = load_Joinville_transect_CTDs()
CTDs_grouped = CTDs.groupby("Event")
events = CTDs_grouped.groups.keys()

# define a common axis with grid spacing of 1
new_mab = np.arange(0, 5000, 1)

# objects for saving the data later
N_df = pd.DataFrame()
gamma_n_df = pd.DataFrame()

for event in events:
    current_profile = CTDs_grouped.get_group(event).reset_index(drop=True)

    _eps, N2 = mx.overturn.eps_overturn(
        current_profile['Depth water [m]'],
        current_profile['Temp [°C]'],
        current_profile['Sal'],
        current_profile["Longitude"].mean(),
        current_profile["Latitude"].mean(),
    )

    N = np.sqrt(N2)  # s* 86400 / (2 * np.pi) # Calculate buoyancy frequency in units of cycles per day (cpd).

    # Plot only in the depth range:
    max_depth = current_profile['Depth water [m]'].max()
    depth = current_profile['Depth water [m]']

    if np.all(np.isnan(N)):
        print(f"{event} only produces NaNs")
        continue

    if max_depth < 200:
        print(f"{event} has {max_depth}m and is too shallow")
        continue
    if current_profile['Temp [°C]'].iloc[-1] > 0.2:
        print(f"{event} shows {current_profile['Temp [°C]'].iloc[-1]:.2f} °C at the bottom, which is too high ")
        continue

    # nearest interpolation to the defined axis
    N_func = interp1d(max_depth - depth, N, kind='nearest', bounds_error=False, fill_value=(np.nan, np.nan))
    gamma_n_func = interp1d(max_depth - depth, current_profile['Neutral density [kg m^-3]'], kind='nearest', bounds_error=False,
                      fill_value=(np.nan, np.nan))

    N_df[current_profile["Longitude"].mean()] = N_func(new_mab)
    gamma_n_df[current_profile["Longitude"].mean()] = gamma_n_func(new_mab)

# sort columns after their longitude value
N_df.sort_index(axis=1, inplace=True)
gamma_n_df.sort_index(axis=1, inplace=True)

gamma_n_df.to_pickle("../../data/Neutral_density_df_with_mab.pkl")

lons = gamma_n_df.columns.to_numpy()
mab = gamma_n_df.index
max_lon = max(lons)
min_lon = min(lons)
# half a degree bins
BIN_EDGES = np.arange(min_lon - 1e-3 * min_lon, 0.5+max_lon + 1e-3 * max_lon, 0.5)
BIN_CENTER = BIN_EDGES[:-1]-0.25

rows = []
for index, row in gamma_n_df.iterrows():
    values = row.to_numpy()
    bin_means = ss.binned_statistic(x=lons, values=values, statistic=np.nanmean, bins=BIN_EDGES)[0]
    new_row = pd.DataFrame([bin_means], columns=BIN_CENTER)
    rows.append(new_row)

binned_gamma_n_df = pd.concat(rows, sort=False).reset_index(drop=True)
binned_gamma_n_df.to_csv("./method_results/binned_gamma_n.csv")
binned_gamma_n_df.to_csv("../../derived_data/binned_neutral_density.csv")

print("done")
