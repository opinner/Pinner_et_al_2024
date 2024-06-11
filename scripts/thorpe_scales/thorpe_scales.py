#!/usr/bin/env python
# coding: utf-8

# Thorpe scale approach in the Weddell Sea

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

plt.rcParams.update({
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "figure.figsize": [8, 6]
})
import mixsea as mx

from src.read_CTDs import read_transect_CTDs
import src.helper as helper
from scipy.interpolate import interp1d  # is considered legacy code, will be in the future removed from scipy
import warnings

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

DENSITY_NOISE = 1e-4  # Noise parameter, Default value = 5e-4
ALPHA = 0.8  # Coefficient relating the Thorpe and Ozmidov scales.
BACKGROUND_EPS = np.nan  # Background value of epsilon applied where no overturns are detected.

CTDs = read_transect_CTDs
CTDs_grouped = CTDs.groupby("Event")
events = CTDs_grouped.groups.keys()

# define a common axis with grid spacing of 1
new_mab = np.arange(0, 5000, 1)

# objects for saving the data later    
eps_df = pd.DataFrame()
N_df = pd.DataFrame()
T_df = pd.DataFrame()

for event in events:
    current_profile = CTDs_grouped.get_group(event).reset_index(drop=True)

    eps, N2 = mx.overturn.eps_overturn(
        current_profile['Depth water [m]'],
        current_profile['Temp [°C]'],
        current_profile['Sal'],
        current_profile["Longitude"].mean(),
        current_profile["Latitude"].mean(),
        dnoise=DENSITY_NOISE,
        alpha=ALPHA,
        background_eps=BACKGROUND_EPS,
    )

    N = np.sqrt(N2)  # s* 86400 / (2 * np.pi) # Calculate buoyancy frequency in units of cycles per day (cpd).

    # Plot only in the depth range:
    max_depth = current_profile['Depth water [m]'].max()
    cut = current_profile['Depth water [m]'] > (max_depth - 500)
    depth = current_profile['Depth water [m]'][cut]

    if np.all(np.isnan(N[cut])):
        print(f"{event} only produces NaNs")
        continue

    if max_depth < 200:
        print(f"{event} has {max_depth}m and is too shallow")
        continue
    if current_profile['Temp [°C]'].iloc[-1] > 0.2:
        print(f"{event} shows {current_profile['Temp [°C]'].iloc[-1]:.2f} °C at the bottom, which is too high ")
        continue

    # nearest interpolation to the defined axis
    N_func = interp1d(max_depth - depth, N[cut], kind='nearest', bounds_error=False, fill_value=(np.nan, np.nan))
    eps_func = interp1d(max_depth - depth, eps[cut], kind='nearest', bounds_error=False, fill_value=(np.nan, np.nan))
    T_func = interp1d(max_depth - depth, current_profile['Temp [°C]'][cut], kind='nearest', bounds_error=False,
                      fill_value=(np.nan, np.nan))

    N_df[current_profile["Longitude"].mean()] = N_func(new_mab)
    eps_df[current_profile["Longitude"].mean()] = eps_func(new_mab)
    T_df[current_profile["Longitude"].mean()] = T_func(new_mab)


# sort columns after their longitude value
eps_df.sort_index(axis=1, inplace=True)
N_df.sort_index(axis=1, inplace=True)
T_df.sort_index(axis=1, inplace=True)


# small data cleaning
try:
    assert T_df.iloc[:, 35].count() < 120
    print(f"Profile at {T_df.columns[35]:.1f}°W contains only {T_df.iloc[:, 35].count()} values and will be removed")
    eps_df.drop(T_df.columns[35], axis="columns", inplace=True)
    N_df.drop(T_df.columns[35], axis="columns", inplace=True)
    T_df.drop(T_df.columns[35], axis="columns", inplace=True)
except AssertionError:
    print("Wrong clean up parameters")
    pass

f, ax = plt.subplots(nrows=1, figsize=(10, 5))
# mab_bin_edges = bin_edges(eps_strain_df.index,dz)
# lon_edges = eps_strain_df.columns - np.diff(eps_strain_df.columns)
mpp = ax.pcolormesh(eps_df.columns, eps_df.index, eps_df,
                    norm= mcolors.LogNorm(vmax = 1e-7),
                    shading="nearest"
                    )
cb = plt.colorbar(mpp, ax=ax)
cb.set_label(r"$\varepsilon$ / (W kg$^{-1}$)")
#ax.set_facecolor('lightgrey')
ax.set_ylim(0,500)
ax.set_title(r"Mixing diagnosed from Strain parametrization")
helper.Plot.path_as_footnote(fig = f,
                             path = "Pinner_et_al_2024/scripts/thorpe_scales/thorpe_scales.py",
                             rot = "vertical")
f.tight_layout()

# cut to only cover the core of the gravity current
vertical_eps_df = eps_df.drop(eps_df.columns[eps_df.columns < -51.5], axis="columns")
vertical_eps_df.drop(vertical_eps_df.columns[vertical_eps_df.columns > -48.5], axis="columns", inplace=True)

# Fill NaN with 'assumed_background_dissipation' only where there is temperature data
assumed_background_dissipation = 1e-10
vertical_eps_df.fillna(value= assumed_background_dissipation, inplace=True)
vertical_eps_df.where(cond=~T_df.isna(), other=np.nan, inplace=True)
mean_profile = vertical_eps_df.mean(axis=1)
std_of_mean_profile = vertical_eps_df.std(axis=1)

# save data
eps_df.to_pickle("./method_data/Thorpe_eps_df_with_mab.pkl")
T_df.to_pickle("./method_data/Thorpe_T_df_with_mab.pkl")
np.savez("./method_data/horizontally_averaged_Thorpe_eps", z=vertical_eps_df.index, eps= mean_profile)

print("done")
plt.show()
