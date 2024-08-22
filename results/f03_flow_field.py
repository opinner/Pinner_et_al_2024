import warnings

import cmocean
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import scipy.interpolate as si
import src.helper as helper

# Suppress specific RuntimeWarning related to mean of empty slice
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*Mean of empty slice.*")

ONE_COLUMN_WIDTH = 8.3
TWO_COLUMN_WIDTH = 12
GOLDEN_RATIO = 1.61
cm = 1/2.54  # centimeters in inches
plt.style.use('./paper.mplstyle')

# load all 7 moorings as dataframes
list_of_moorings = helper.IO.load_pickle(name="../data/mooring/list_of_moorings.pkl")
data = np.load("../data/max_depth_dict.npz", allow_pickle=True)
max_depth_dict = data["max_depth_dict"].item()
neutral_density = pd.read_csv("../scripts/preprocessing/method_results/binned_gamma_n.csv")
neutral_density.set_index(keys = 'Unnamed: 0', inplace = True)
neutral_density = neutral_density.drop(neutral_density[neutral_density.index > 550].index)

def add_neutral_density_lines(ax, binned_gamma_n_df):
    # mpp = ax.pcolormesh(
    #     binned_thorpe_gamma_n_df.columns,
    #     mab,
    #     binned_thorpe_gamma_n_df,
    #     cmap=cmocean.cm.haline_r,
    #     vmin=27.8,
    #     rasterized=True # optimize the drawing for vector graphics
    # )

    #cb = plt.colorbar(mpp, ax=ax, extend="min", location="top")  # pad=0.02, aspect=12

    water_mass_boundaries = [28.26, 28.40]  # + 28.00 boundary
    # gravity_current_boundary = [28.40]  # from Garabato et al 2002
    CS = ax.contour(
        binned_gamma_n_df.columns,
        binned_gamma_n_df.index,
        binned_gamma_n_df,
        levels=water_mass_boundaries,
        linestyles=["dashed", "solid"],
        colors="k",
        linewidths=3,
        zorder=50
    )

    fmt = {}
    strs = ['WSDW', 'WSBW']
    for l, s in zip(CS.levels, strs):
        fmt[l] = s

    # Label every other level using strings
    ax.clabel(
        CS,
        CS.levels,
        inline=False,
        fmt=fmt,
        fontsize=10,
        colors="white"
    )


def vertical_then_horizontal_interpolation(x, y, z):
    """
    Perform a two-step interpolation: vertical first, then horizontal.

    Parameters:
    - x: 1D array of x-coordinates
    - y: 1D array of y-coordinates
    - z: 1D array of z-values

    Returns:
    - xi, yi, zi: grid arrays for plotting the interpolated surface
    """

    # Ensure input arrays are numpy arrays
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)

    # Sort the data by x and y for structured interpolation
    sorted_indices = np.lexsort((y, x))  # Sort primarily by x, then by y
    x_sorted = x[sorted_indices]
    y_sorted = y[sorted_indices]
    z_sorted = z[sorted_indices]

    # Create unique x and y values
    unique_x = np.unique(x_sorted)
    # unique_y = np.unique(y_sorted)

    # Define a new grid for interpolation
    xi = np.linspace(min(x_sorted), max(x_sorted), 100)
    yi = np.linspace(min(y_sorted), max(y_sorted), 100)

    # Step 1: Interpolate Vertically for each x value
    z_vertical_interp = np.zeros((len(yi), len(unique_x)))

    for i, x_val in enumerate(unique_x):
        # Select points for a particular x
        mask = (x_sorted == x_val)
        y_subset = y_sorted[mask]
        z_subset = z_sorted[mask]

        # Check if we have enough points for interpolation
        if len(y_subset) > 1:
            # Create a vertical interpolation function
            f_vertical = si.interp1d(y_subset, z_subset, kind='linear', fill_value="extrapolate")
            # Interpolate for all y values
            z_vertical_interp[:, i] = f_vertical(yi)
        else:
            # If not enough points, fill with NaN
            z_vertical_interp[:, i] = np.nan

    # Step 2: Interpolate Horizontally for each y in the vertical interpolation result
    zi = np.zeros((len(yi), len(xi)))

    for j, y_val in enumerate(yi):
        # Filter out NaN values from previous interpolation
        valid_indices = ~np.isnan(z_vertical_interp[j, :])
        valid_x = unique_x[valid_indices]
        valid_z = z_vertical_interp[j, valid_indices]

        # Check if we have enough points for interpolation
        if len(valid_x) > 1:
            # Create a horizontal interpolation function
            f_horizontal = si.interp1d(valid_x, valid_z, kind='linear', fill_value="extrapolate")
            # Interpolate for all x values
            zi[j, :] = f_horizontal(xi)
        else:
            # If not enough points, fill with NaN
            zi[j, :] = np.nan

    return xi, yi, zi


fig1, ax1 = plt.subplots(1, figsize=(TWO_COLUMN_WIDTH*cm, 0.5*TWO_COLUMN_WIDTH*cm))
fig2, ax2 = plt.subplots(1, figsize=(TWO_COLUMN_WIDTH*cm, 0.5*TWO_COLUMN_WIDTH*cm))

avrg_velos = []
max_velos = []
lon_velos = []
mab_velos = []
for nr, mooring in enumerate(list_of_moorings):
    for measurement_depth in mooring.columns:
        if measurement_depth == "time":
            continue

        mab_of_measurement = int(max_depth_dict[mooring.location.lon]) - int(measurement_depth)
        if mab_of_measurement > 600:
            continue

        complex_velocity = mooring[measurement_depth]
        avrg_velocity = np.abs(np.nanmean(complex_velocity))
        max_velocity = np.nanmax(np.abs(complex_velocity))
        #print(avrg_velocity, max_velocity)
        #print(type(mooring.location.lon), type(mab_of_measurement),"\n")
        # ax1.plot(mooring.location.lon, mab_of_measurement, "k.",
        #          alpha=0)
        # ax2.plot(mooring.location.lon, mab_of_measurement, "k.",
        #          alpha=0)  #, s= f"{max_velocity:.2f}", color = "white")
        ax1.text(mooring.location.lon, mab_of_measurement,
                 s=f"{max_velocity:.2f}",
                 color="k",
                 zorder=10,
                 ha="center",
                 va="center",
                 fontsize = 7,
                 bbox = dict(facecolor='white', alpha=0.8, boxstyle='round', pad=0.15)
                 )
        ax2.text(mooring.location.lon, mab_of_measurement,
                 s=f"{max_velocity:.2f}",
                 color="k",
                 zorder=10,
                 ha="center",
                 va="center"
                 )

        avrg_velos.append(avrg_velocity)
        max_velos.append(max_velocity)
        lon_velos.append(mooring.location.lon)
        mab_velos.append(mab_of_measurement)

xi1, yi1, zi_avrg = vertical_then_horizontal_interpolation(lon_velos, mab_velos, avrg_velos)
xi2, yi2, zi_max = vertical_then_horizontal_interpolation(lon_velos, mab_velos, max_velos)

levels = np.arange(0, np.max(avrg_velos)+0.025, 0.025)
#print(levels)
mpp = ax1.contourf(xi1, yi1, zi_avrg,
             levels=levels,
             cmap='viridis'
             )
cb = plt.colorbar(mpp, ax = ax1)
cb.set_label(r'Mean Velocity (m$\,$s$^{-1}$)')

ax2.contourf(xi2, yi2, zi_max,
             levels=15,
             cmap='viridis')

water_mass_boundaries = [28.26, 28.40]  # + 28.00 boundary
# gravity_current_boundary = [28.40]  # from Garabato et al 2002
CS = ax1.contour(
    np.array(neutral_density.columns).astype(float),
    np.array(neutral_density.index).astype(float),
    neutral_density,
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
ax1.clabel(
    CS,
    CS.levels,
    inline=False,
    fmt=fmt,
    fontsize=9,
    colors="white"
)

for x, label in zip(np.unique(lon_velos), ["A","B","C","D","E","F","G"]):
    ax1.text(x, 405,
             s=label,
             color="k",
             zorder=10,
             ha="center",
             va="center",
             fontsize=9,
             )

ax1.set_facecolor('lightgrey')
ax1.set_ylabel("Meters above bottom")
ax1.set_xlabel("Longitude (°)")
ax1.set_ylim(-5, 380)
ax2.set_ylim(-5, 380)
ax1.set_xlim(-52.5, -47.3)
ax2.set_xlim(-52.5, -47.3)
fig1.tight_layout()
fig1.savefig("./flowfield.pdf")
plt.show()
