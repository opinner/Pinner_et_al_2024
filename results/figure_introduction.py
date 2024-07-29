import warnings

import cartopy.crs as ccrs
import cmocean as cm
import matplotlib.colors as mcolors
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rioxarray
import scipy.stats as ss
from matplotlib.gridspec import GridSpec

# Suppress specific RuntimeWarning related to mean of empty slice
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*Mean of empty slice.*")

import src.helper as helper
from src.read_CTDs import load_Joinville_transect_CTDs

ONE_COLUMN_WIDTH = 8.3
TWO_COLUMN_WIDTH = 12
GOLDEN_RATIO = 1.61
CM = 1 / 2.54  # centimeters in inches


def setup_curved_map_shape(ax, min_lon, max_lon, min_lat, max_lat):
    """
    Sets up a high-resolution South Polar Stereographic projection map with specified latitude and longitude extents.

    This function creates a plot with a South Polar Stereographic projection. It defines the projection, sets the extent
    of the map, and adjusts the boundary of the map to fit within the specified latitude and longitude limits. The resolution
    for computing great circles is increased for finer details.

    Args:
        ax (matplotlib.axes._subplots.AxesSubplot): The axis object to set up.
        min_lon (float): The minimum longitude of the map extent.
        max_lon (float): The maximum longitude of the map extent.
        min_lat (float): The minimum latitude of the map extent.
        max_lat (float): The maximum latitude of the map extent.

    Returns:
        None
    """

    ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())

    # Define the limits for the x-axis and y-axis of the map
    xlim = [min_lon, max_lon]  # longitude
    ylim = [min_lat, max_lat]  # latitude

    # Define a rectangular path for the map extent using matplotlib's Path object
    rect = mpath.Path(
        [[xlim[0], ylim[0]], [xlim[1], ylim[0]], [xlim[1], ylim[1]], [xlim[0], ylim[1]], [xlim[0], ylim[0]]])

    # Transform the rectangular path to the data coordinate system
    proj_to_data = ccrs.PlateCarree()._as_mpl_transform(ax) - ax.transData
    rect_in_target = proj_to_data.transform_path(rect)

    # Set the boundary of the map using the transformed rectangular path
    ax.set_boundary(rect_in_target)

    return None


def plot_shaded_topography(ax: object, subset: object, param_dict: dict) -> None:
    """
    Plots shaded topography and bathymetry on a map with hill shading and contours.

    This function adds shaded relief and bathymetry to an existing map plot. It uses hill shading to create a 3D effect
    and overlays bathymetric data. Contour lines are added to highlight specific depth levels in the bathymetric data.

    Args:
        ax (matplotlib.axes._subplots.AxesSubplot): The axis object to plot on.
        subset (xarray.DataArray): The subset of bathymetric data to plot.
        param_dict (dict): Dictionary containing parameters for the plot (e.g., "topo_shading", "cb_plot", "cb_location").

    Returns:
        None
    """
    # Define direction of light source
    ls = mcolors.LightSource(azdeg=315, altdeg=45)

    # extract the single data layer
    subset_data = subset[0, :, :]
    left, bottom, right, top = subset_data.rio.bounds()

    if param_dict["topo_shading"] == True:
        # Plot and shade a base layer
        shade_only = np.array(subset_data)
        shade_only = shade_only.astype(float)  # Ensure float type

        # Compute hill shading
        shade = ls.hillshade(shade_only, vert_exag=0.1)

        # Plot the shaded relief using imshow
        # interpolation = 'none' is needed for good-looking exports to svg or pdf
        ax.imshow(shade, extent=(left, right, bottom, top), transform=ccrs.PlateCarree(), origin='upper',
                  cmap=cm.cm.gray,
                  interpolation='none')

    # Use imshow to plot the bathymetry (negative values)
    bathymetry = np.where(subset_data >= 0, np.nan, subset_data)
    im = ax.imshow(bathymetry, extent=(left, right, bottom, top),
                   transform=ccrs.PlateCarree(),
                   origin='upper', cmap=cm.cm.ice,
                   vmin=-5000,
                   interpolation='none')

    # Prepare data for hill shading only below 0 with transparency
    shade_only = np.where(subset_data > 0, np.nan, subset_data)
    shade_only = shade_only.astype(float)  # Ensure float type

    # Compute hill shading with transparency
    shade = ls.hillshade(shade_only, vert_exag=0.01, fraction=0.01)

    # Create a colormap that goes from fully transparent to black with some transparency
    transparent = mcolors.colorConverter.to_rgba('white', alpha=0)
    black = mcolors.colorConverter.to_rgba('black', alpha=0.3)
    custom_shading_cmap = mcolors.LinearSegmentedColormap.from_list('rb_cmap', [transparent, black], 512)

    # Use imshow to shade bathymetry
    ax.imshow(shade, extent=(left, right, bottom, top), transform=ccrs.PlateCarree(), origin='upper',
              cmap=custom_shading_cmap, interpolation='none')

    # Extract coordinates for contour plotting
    lon, lat = subset.x.values, subset.y.values
    lon_grid, lat_grid = np.meshgrid(lon, lat)

    # Create contour lines for specified depths, in the bathymetry map as well as in the colormap
    levels = [-4000, -3000, -2000, -1000]
    contour = ax.contour(lon_grid, lat_grid, subset_data, transform=ccrs.PlateCarree(),
                         levels=levels, colors="white", linestyles="solid", linewidths=0.5)

    if param_dict["cb_plot"] == True:
        # Add colorbar for bathymetry
        cbar = plt.colorbar(im, ax=ax, location=param_dict["cb_location"])
        cbar.set_label('Bathymetry (m)')
        cbar.add_lines(levels, colors=len(levels) * ["w"], linewidths=1)

    return None


def overview_map(ax, ibsco_data):
    """
    Creates an overview map of the Weddell Sea using IBSC0 data.

    Args:
        ax (matplotlib.axes._subplots.AxesSubplot): The axis object to plot on.
        ibsco_data (xarray.DataArray): The IBSC0 bathymetric data.

    Returns:
        None
    """

    # Define the bounds in lat, lon coordinates for the subset
    min_lon = -65
    max_lon = -15
    min_lat = -80
    max_lat = -60

    # Clip the raster using the defined bounds
    subset = ibsco_data.rio.clip_box(minx=min_lon, miny=min_lat, maxx=max_lon, maxy=max_lat)

    setup_curved_map_shape(ax, min_lon, max_lon, min_lat, max_lat)
    cb_dict = {"topo_shading": True, "cb_plot": False, "cb_location": "bottom"}
    plot_shaded_topography(ax, subset, cb_dict)

    gl = ax.gridlines(
        xlocs=np.arange(-60, -10, 10),
        ylocs=np.arange(-80, -60, 4),
        # This is commented out as it does not produce the desired results
        # as not all grid lines will be labelled
        # draw_labels=["top", "left", "x", "y"],
        dms=True,
        x_inline=False,
        y_inline=False
    )

    # draw scale bar
    # gvpy.maps.cartopy_scale_bar(ax0, location = (0.1,0.76), length = 500)

    ax.plot([-55, -45], [-63, -64],
            color='tab:red', lw=2,
            transform=ccrs.Geodetic())

    return None


def joinville_transect(ax, ibsco_data):
    """
    Plots a detailed map of the Joinville transect with bathymetric data and mooring/CTD locations.

    Args:
        ax (matplotlib.axes._subplots.AxesSubplot): The axis object to plot on.
        ibsco_data (xarray.DataArray): The IBSC0 bathymetric data.

    Returns:
        None
    """

    # Define the bounds in lat, lon coordinates for the subset
    min_lon = -56
    max_lon = -44
    min_lat = -67
    max_lat = -61

    # ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
    ax.set_extent([-55, -45, -65, -63], crs=ccrs.PlateCarree())

    # Clip the raster using the defined bounds
    subset = ibsco_data.rio.clip_box(minx=min_lon, miny=min_lat, maxx=max_lon, maxy=max_lat)
    param_dict = {"topo_shading": False, "cb_plot": True, "cb_location": "top"}
    plot_shaded_topography(ax, subset, param_dict)

    gl = ax.gridlines(xlocs=np.arange(min_lon, max_lon, 2), ylocs=np.arange(min_lat, max_lat, 1),
                      draw_labels={"bottom": "x", "right": "y"},
                      rotate_labels=False,
                      y_inline=False)

    # load all 7 moorings as dataframes
    list_of_moorings = helper.IO.load_pickle(name="../data/mooring/list_of_moorings.pkl")
    for mooring in list_of_moorings:
        ax.plot(mooring.location.lon, mooring.location.lat, "D", c="tab:red", markersize=10,
                transform=ccrs.PlateCarree())

    ax.plot(-mooring.location.lon, -mooring.location.lat, "D", c="tab:red", markersize=10,
            transform=ccrs.PlateCarree(), label="Moorings")

    CTDs = load_Joinville_transect_CTDs()
    unique_coords_df = CTDs.drop_duplicates(subset=["Latitude", "Longitude"])

    color = "grey"
    markersize = 10
    ax.plot(unique_coords_df["Longitude"], unique_coords_df["Latitude"],
            ".",
            markersize=markersize,
            color=color,
            alpha=0.8,
            markeredgewidth=0,
            # markeredgecolor = color,
            transform=ccrs.PlateCarree(),
            label="CTD profiles"
            )

    ax.legend()

    # gvpy.maps.cartopy_scale_bar(ax, location=(0.04, 0.6), length=50)

    # Turn off some axis ticks
    # ax1.tick_params(axis='both', which='both', left=True, right=False, bottom=True, top=False, labelbottom=True)


def neutral_density_cross_section(ax):
    """
    Plots neutral density across the slope plus the location of current velocity measurements

    Args:
        ax (matplotlib.axes._subplots.AxesSubplot): The axis object to plot on.

    Returns:
        None
    """
    thorpe_gamma_n_df = pd.read_pickle("../scripts/thorpe_scales/method_results/Thorpe_neutral_density_df_with_mab.pkl")
    lons = thorpe_gamma_n_df.columns.to_numpy()
    mab = thorpe_gamma_n_df.index
    max_lon = max(lons)
    min_lon = min(lons)
    # half a degree bins
    BIN_EDGES = np.arange(min_lon - 1e-3 * min_lon, max_lon + 1e-3 * max_lon, 0.5)

    rows = []
    for index, row in thorpe_gamma_n_df.iterrows():
        values = row.to_numpy()
        bin_means = ss.binned_statistic(x=lons, values=values, statistic=np.nanmean, bins=BIN_EDGES)[0]
        new_row = pd.DataFrame([bin_means], columns=BIN_EDGES[:-1])
        rows.append(new_row)
    binned_thorpe_gamma_n_df = pd.concat(rows, sort=False).reset_index(drop=True)

    mpp = ax.pcolormesh(
        binned_thorpe_gamma_n_df.columns,
        mab,
        binned_thorpe_gamma_n_df,
        cmap=cm.cm.haline_r,
        vmin=27.8)

    cb = plt.colorbar(mpp, ax=ax, extend="min", location="top")  # pad=0.02, aspect=12

    water_mass_boundaries = [28.26, 28.40]  # + 28.00 boundary
    # gravity_current_boundary = [28.40]  # from Garabato et al 2002
    ax.contour(
        binned_thorpe_gamma_n_df.columns,
        binned_thorpe_gamma_n_df.index,
        binned_thorpe_gamma_n_df,
        levels=water_mass_boundaries,
        linestyles=["dashed", "solid"],
        colors="k",
        linewidths=3,
    )

    # ax.annotate('bottom\ncurrent', xy=(-51.69, 184), xytext=(-51.8, 230),
    #              arrowprops=dict(facecolor='black', width=2, shrink=0.05), ha="center", color="white",
    #              bbox=dict(facecolor='black', alpha=0.8, edgecolor='black', boxstyle='round'))
    # ax.annotate('gravity current\nboundary', xy=(-48.8, 130), xytext=(-48, 270),  # fontsize=9,
    #             arrowprops=dict(facecolor='black', width=2, shrink=0.05), ha="center", va="center", color="white",
    #             bbox=dict(facecolor='black', alpha=0.8, edgecolor='black', boxstyle='round, pad = 0.5'))

    mooring_info = pd.read_csv("../scripts/IDEMIX_parametrization/method_results/eps_IGW_IDEMIX_results.csv")
    moorings_mabs = mooring_info["rounded_mab"]
    moorings_lons = mooring_info["lon"]

    # draw measurement positions
    ax.plot(moorings_lons, moorings_mabs,
            "D",
            label="rotor current meters",
            color="tab:red",
            markersize=10,
            markeredgecolor="k",
            zorder=6)

    ax.set_ylim(-10, 500)
    xlim = ax.get_xlim()
    ax.set_xlim((xlim[0] - 0.2, xlim[1] + 0.8))
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Meters above Seafloor")
    cb.set_label(r"Neutral Density $\gamma_n\,$(kg$\,$m$^{-3}$)")
    cb.ax.plot([water_mass_boundaries[0], water_mass_boundaries[0]], [0, 1], 'k--', lw=2)
    cb.ax.plot([water_mass_boundaries[1], water_mass_boundaries[1]], [0, 1], 'k', lw=2, ls="solid")

    # draw sea floor
    x = list(ax.get_xlim())
    y1 = [0, 0]
    y2 = 2 * list(ax.get_ylim())[0]
    ax.axhline(0, c="k")
    ax.fill_between(x, y1, y2, facecolor="xkcd:charcoal grey", zorder=5)  # , hatch="///")
    ax.legend(loc="upper right")  # ,facecolor='k', framealpha=0.8, edgecolor = "black", labelcolor = "white")


if __name__ == "__main__":
    input_file = '../data/bathymetry/IBCSO_v2_ice-surface_WGS84.tif'
    ibsco_data = rioxarray.open_rasterio(input_file)

    # fig = plt.figure(layout="constrained", figsize=(TWO_COLUMN_WIDTH * CM, 0.8 * TWO_COLUMN_WIDTH * CM))
    fig = plt.figure(figsize=(TWO_COLUMN_WIDTH * CM, 0.8 * TWO_COLUMN_WIDTH * CM))
    gs = GridSpec(nrows=2, ncols=2, figure=fig)
    high_res_proj = ccrs.SouthPolarStereo(central_longitude=-40)
    high_res_proj.threshold = 1e3
    ax0 = fig.add_subplot(gs[0, 0], projection=high_res_proj)
    ax1 = fig.add_subplot(gs[0, 1], projection=ccrs.SouthPolarStereo(central_longitude=-50))
    # identical to ax1 = plt.subplot(gs.new_subplotspec((0, 0), colspan=3))
    ax2 = fig.add_subplot(gs[1:, :])

    overview_map(ax0, ibsco_data)
    joinville_transect(ax1, ibsco_data)
    neutral_density_cross_section(ax2)

    fig.tight_layout()

    # plt.savefig(f"./figures/Intro_figure.png", dpi=800, bbox_inches='tight')
    # plt.savefig(f"./figures/Intro_figure.svg", bbox_inches='tight')
    plt.show()
