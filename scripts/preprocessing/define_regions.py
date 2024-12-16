import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cmocean

# load binned data
binned_neutral_density = pd.read_csv("../../derived_data/binned_neutral_density.csv", index_col=0)
# convert column names from strings to floats
binned_neutral_density.columns = binned_neutral_density.columns.astype("float")

# define regions
regions_key = {"air": -1, "shelf": 0, "open ocean": 1, "IL": 2, "BL": 3}

# interfacial layer from density greater than WSBW density
binned_regions = binned_neutral_density.where(cond=binned_neutral_density < 28.40, other=regions_key["IL"])

# bottom layer from difference to bottom density
bottom_density = binned_neutral_density.iloc[0]
grav_curr_bottom_density = bottom_density.where(cond=bottom_density > 28.40, other=np.nan)
binned_regions = binned_regions.where(cond=(grav_curr_bottom_density - binned_neutral_density) > 0.01,
                                      other=regions_key["BL"])

# open ocean from density smaller than WSBW density
binned_regions = binned_regions.where(cond=binned_neutral_density > 28.40, other=regions_key["open ocean"])


# shelf from longitude along the transect (the column names)
def condition(col_name): return col_name > -52
# condition = lambda col_name: col_name > -52


# Create the mask
mask = pd.DataFrame(
    {col: condition(col) for col in binned_regions.columns},  # Check condition on column names
    index=binned_regions.index  # Retain the same index
)
binned_regions = binned_regions.where(cond=mask, other=regions_key["shelf"])

# air/na where no neutral density is defined
binned_regions = binned_regions.where(cond=~binned_neutral_density.isna(), other=regions_key["air"])

# save regions
binned_regions.to_csv("./method_results/binned_regions.csv")
binned_regions.to_csv("../../derived_data/gravity_current_regions.csv")

# plot regions
binned_regions = binned_regions.iloc[0:600]
levels = np.arange(-1.5, 4.5, 1)
plt.pcolormesh(
    binned_regions.columns,
    binned_regions.index,
    binned_regions.values,
    norm=mcolors.BoundaryNorm(levels, ncolors=256),
    cmap=cmocean.cm.rain
)
cbar = plt.colorbar(ticks=list(regions_key.values()))
cbar.set_ticks(ticks=list(regions_key.values()), labels=list(regions_key.keys()))
cbar.ax.invert_yaxis()
plt.show()
