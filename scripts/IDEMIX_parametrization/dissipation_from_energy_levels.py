import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({
    "figure.facecolor":  "white",
    "savefig.facecolor": "white",
})
import scipy.signal as signal

import src.helper as helper
from src.mooring import Mooring
from src.location import Location
from src.ctd_cast import CTDCast
import equations as eq

N_table = pd.read_pickle("./method_data/N_values.pkl")
N_error_table = pd.read_pickle("./method_data/N_std.pkl")

data = np.load("./method_data/results_available_energy.npz", allow_pickle = True)
energy_levels = pd.DataFrame(data = {
    "lon":data["lon"],
    "lat":data["lat"],
    "rounded_depth": data["depth"],
    "rounded_mab": data["mab"],
    "barotropic":data["barotropic"],
    "continuum":data["continuum"],
    "available":data["available"],
})
energy_levels["rounded_depth"] = energy_levels["rounded_depth"].astype("int")
print("Loaded wave energy table")

# Assign N and its error to each location, where the wave energy level is know
N_array = []
N_error_array = []
for index, row in energy_levels.iterrows():
    column_name = f"({row['lat']:.2f},{row['lon']:.2f})"
    # print(column_name)
    try:
        N_value = N_table.loc[N_table['mab'] == row["rounded_mab"], column_name].item()
        N_error = N_error_table.loc[N_error_table['mab'] == row["rounded_mab"], column_name].item()
    except ValueError:
        N_value = np.nan
        N_error = np.nan
    # print(N_error)
    N_array.append(N_value)
    N_error_array.append(N_error)

energy_levels["N"] = N_array
energy_levels["N Error"] = N_error_array
print("Added N values")

# get coriolis frequency at the geographic location of every mooring
energy_levels["coriolis frequency"] = energy_levels.apply(lambda row: helper.Constants.get_coriolis_frequency(row["lat"], unit = "rad/s"), axis=1)
print("Added f values")


