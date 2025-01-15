import numpy as np
import pandas as pd
from src.read_CTDs import load_Joinville_transect_CTDs
import scipy.stats as ss
import warnings

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

CTDs = load_Joinville_transect_CTDs()
CTDs_grouped = CTDs.groupby("Event")
events = CTDs_grouped.groups.keys()

lats = []
lons = []

for event in events:
    current_profile = CTDs_grouped.get_group(event).reset_index(drop=True)
    lons.append(current_profile["Longitude"].mean())
    lats.append(current_profile["Latitude"].mean())

max_lon = max(lons)
min_lon = min(lons)
# half a degree bins
LON_BIN_EDGES = np.arange(-53.75, -46.25, 0.5)
#LON_BIN_EDGES = np.arange(min_lon - 1e-3 * min_lon, 0.5+max_lon + 1e-3 * max_lon, 0.5)
bin_lats = ss.binned_statistic(x=lons, values=lats, statistic=np.nanmean, bins=LON_BIN_EDGES)[0]
bin_lons = LON_BIN_EDGES[:-1] + 0.25

print(bin_lons)
print(len(bin_lats),len(bin_lons))

bin_coordinates = pd.DataFrame(data={"bin_lats": bin_lats, "bin_lons":bin_lons})
bin_coordinates.to_csv("./method_results/bin_coordinates.csv", index=False)
bin_coordinates.to_csv("../../derived_data/bin_coordinates.csv", index=False)