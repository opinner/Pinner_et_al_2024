import rioxarray
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.path as mpath
import cmocean
import cartopy.crs as ccrs
from matplotlib.gridspec import GridSpec

ONE_COLUMN_WIDTH = 8.3
TWO_COLUMN_WIDTH = 12
GOLDEN_RATIO = 1.61
cm = 1/2.54  # centimeters in inches

def overview_map(fig,ax):
    input_file = '../data/bathymetry/IBCSO_v2_ice-surface_WGS84.tif'
    rds = rioxarray.open_rasterio(input_file)

    # Define the bounds in lat, lon coordinates for the subset
    min_lon = -65
    min_lat = -80
    max_lon = -20
    max_lat = -60

    # Clip the raster using the defined bounds
    subset = rds.rio.clip_box(minx=min_lon, miny=min_lat, maxx=max_lon, maxy=max_lat)
    left, bottom, right, top = subset.rio.bounds()
    subset_data = subset[0, :, :]

    # Create the plot with a South Polar Stereographic projection
    #fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection=ccrs.SouthPolarStereo(central_longitude=-50)))
    ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())

    # Define parameters for hill shading
    ls = mcolors.LightSource(azdeg=315, altdeg=45)

    # Prepare data for hill shading only above 0
    shade_only = np.array(subset_data)
    shade_only = shade_only.astype(float)  # Ensure float type

    # Compute hill shading
    shade = ls.hillshade(shade_only, vert_exag=0.1)

    # Plot the shaded relief using imshow
    ax.imshow(shade, extent=(left, right, bottom, top), transform=ccrs.PlateCarree(), origin='upper', cmap=cmocean.cm.gray)

    # Use imshow to plot the bathymetry (negative values)
    bathymetry = np.where(subset_data >= 0, np.nan, subset_data)
    im = ax.imshow(bathymetry, extent=(left, right, bottom, top), transform=ccrs.PlateCarree(), origin='upper',
                   cmap=cmocean.cm.ice, vmin=-5000)
    cbar = plt.colorbar(im, ax=ax, orientation='vertical')  # , shrink = 0.4);
    cbar.set_label('Bathymetry (m)')

    # Prepare data for hill shading only below 0 with transparency
    shade_only = np.where(subset_data > 0, np.nan, subset_data)
    shade_only = shade_only.astype(float)  # Ensure float type

    # Compute hill shading with transparency
    shade = ls.hillshade(shade_only, vert_exag=0.01, fraction=0.01)

    # Create a colormap that goes from fully transparent to black with only some transparency
    c_white = mcolors.colorConverter.to_rgba('white', alpha=0)
    c_black = mcolors.colorConverter.to_rgba('black', alpha=0.3)
    cmap_rb = mcolors.LinearSegmentedColormap.from_list('rb_cmap', [c_white, c_black], 512)

    # Use imshow to shade bathymetry
    ax.imshow(shade, extent=(left, right, bottom, top), transform=ccrs.PlateCarree(), origin='upper', cmap=cmap_rb)

    # Extract coordinates for contour plotting
    lon, lat = subset.x.values, subset.y.values
    lon_grid, lat_grid = np.meshgrid(lon, lat)

    # Create contour lines for specified depths, in the bathymetry map as well as in the colormap
    levels = [-4000, -3000, -2000, -1000]
    contour = ax.contour(lon_grid, lat_grid, subset_data, transform=ccrs.PlateCarree(),
                         levels=levels, colors="white", linestyles="solid", linewidths=0.5)
    cbar.add_lines(levels, colors=len(levels) * ["w"], linewidths=1)


    # Define the limits for the x-axis and y-axis of the map
    xlim = [min_lon, max_lon]  # longitude
    ylim = [min_lat, max_lat]  # latitude
    # Define a rectangular path for the map extent using matplotlib's Path object
    rect = mpath.Path([[xlim[0], ylim[0]],
                       [xlim[1], ylim[0]],
                       [xlim[1], ylim[1]],
                       [xlim[0], ylim[1]],
                       [xlim[0], ylim[0]],
                       ])

    # Transform the rectangular path to the data coordinate system
    proj_to_data = ccrs.PlateCarree()._as_mpl_transform(ax0) - ax0.transData
    rect_in_target = proj_to_data.transform_path(rect)

    # Set the boundary of the map using the transformed rectangular path
    ax0.set_boundary(rect_in_target)

    """
    # Configure gridlines for the map
    gl = ax0.gridlines(
        draw_labels=["top", "right", "x", "y"],
        rotate_labels=True, x_inline=False)

    # Draw latitude and longitude gridlines (and their labels)
    # only at the specified values
    lathelp = np.arange(-75, -60, 5)
    gl.ylocator = mticker.FixedLocator(lathelp)
    lonhelp = np.arange(-80, -10, 10)
    gl.xlocator = mticker.FixedLocator(lonhelp)

    # Rotate the longitude labels to be more readable
    gl.xlabel_style = {'rotation': 30}
    """

    # To make sure that all of the wanted map extent is visible in the figure,
    # (especially because the chosen map projections differs strongly
    # at the poles from the PlateCarree Projection)
    # this hack of adding a constant is needed and may need be adjusted,
    # dependent on your map extent and figure size
    ax0.set_extent([xlim[0], xlim[1], ylim[0] + 4, ylim[1]])

    # draw scale bar
    # gvpy.maps.cartopy_scale_bar(ax0, location = (0.1,0.76), length = 500)

    ax0.plot([-55, -45], [-63, -64],
             color='tab:red', lw=2,
             transform=ccrs.Geodetic())


def joinville_transect(fig,ax):
    gl = ax1.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlabel_style = {'rotation': 0}

    ax1.set_extent([-55, -45, -63, -64.3], ccrs.PlateCarree())

    # ax1.set_title("Mooring Locations in the Weddell Sea, Antarctica")
    # divnorm = mcolors.TwoSlopeNorm(vcenter=0, vmin = -4000, vmax = 1000)
    bat_max, bat_min = 0, -5000
    bat_z[bat_z > bat_max] = np.nan
    bat_z[bat_z < bat_min] = np.nan

    levels = np.arange(-5000, 500, 500)
    image = ax1.contourf(bat_lon, bat_lat, bat_z, levels, transform=ccrs.PlateCarree(), cmap=cmocean.cm.deep_r,
                         vmax=bat_max, vmin=bat_min)
    plt.colorbar(image, ax=ax1, orientation="vertical", aspect=15)  # , shrink=0.5)
    ax1.contour(bat_lon, bat_lat, bat_z, levels, linewidths=0.5, linestyles='solid', colors=['black'], alpha=0.4,
                vmax=bat_max, vmin=bat_min, transform=ccrs.PlateCarree())

    for mooring in list_of_moorings:
        ax1.plot(mooring.location.lon, mooring.location.lat, "D", c="tab:red", markersize=15,
                 transform=ccrs.PlateCarree())

    ax1.plot(-mooring.location.lon, -mooring.location.lat, "D", c="tab:red", markersize=10,
             transform=ccrs.PlateCarree(), label="Moorings")

    color = "grey"
    markersize = 10

    for event in events:
        profile = CTDs_grouped.get_group(event).reset_index(drop=True)
        ax1.plot(profile["Longitude"].iloc[0],
                 profile["Latitude"].iloc[0],
                 ".",
                 markersize=markersize,
                 color=color,
                 alpha=0.8,
                 # markeredgewidth=1,
                 # markeredgecolor = "black",
                 transform=ccrs.PlateCarree())

    ax1.plot(-profile["Longitude"].iloc[0],
             -profile["Latitude"].iloc[0],
             ".",
             markersize=markersize,
             color=color,
             alpha=0.8,
             # markeredgewidth=1,
             # markeredgecolor = "black",
             transform=ccrs.PlateCarree(), label="CTD profiles")

    ax1.legend()

    gvpy.maps.cartopy_scale_bar(ax1, location=(0.04, 0.6), length=50)

    # Turn off some axis ticks
    # ax1.tick_params(axis='both', which='both', left=True, right=False, bottom=True, top=False, labelbottom=True)


def stratification(fig,ax):
    # load N data
    N_table = pd.read_pickle("/home/ole/Desktop/Paper1_Figures/data/poster_N_values.pkl")
    ax2.set_facecolor('grey')

    # Draw temperature background
    bounds = np.arange(-1.1, -0.2, 0.1)
    tric = ax2.tricontourf(lon_array, z_array, temperature_array, levels=bounds, extend="max", cmap=cmocean.cm.thermal)
    level = -0.7
    ax2.tricontour(lon_array, z_array, temperature_array, levels=[level], colors="black")
    ax2.annotate('bottom\ncurrent', xy=(-51.69, 184), xytext=(-51.8, 230),
                 arrowprops=dict(facecolor='black', width=2, shrink=0.05), ha="center", color="white",
                 bbox=dict(facecolor='black', alpha=0.8, edgecolor='black', boxstyle='round'))

    # draw measurement positions
    ax2.plot(uv_lon, uv_mab + 5, "s", label="velocity\nmeasurements", color="tab:red", markersize=6,
             markeredgecolor="k", zorder=5)
    ax2.set_ylim(-10, 390)
    xlim = ax2.get_xlim()
    ax2.set_xlim((xlim[0] - 0.2, xlim[1] + 0.8))
    ax2.set_xlabel("Longitude (°)")
    ax2.set_ylabel("Meters above Seafloor")
    cb = plt.colorbar(tric, ax=ax2, pad=0.02, aspect=12)
    cb.set_label(r"In-situ temperature (°C)")
    cb.ax.plot([0, 1], [-0.7, -0.7], 'k--', lw=2)

    # draw N profiles
    factor = 600  # exaggeration in the x-direction
    columns = list(N_table.columns)[2:]
    lons = sorted(list(set(uv_lon))[1:])
    for column, l in zip(columns, lons):
        if column == "mab": continue
        # print(l,column)
        mean = factor * N_table[column][:100].mean()
        ax2.plot(factor * N_table[column][:200] + l - factor * 0.0005, N_table["mab"][:200], c="w", lw=2)

        # only once for the legend, not visible, not part of the loop
    ax2.plot(factor * N_table[column][:200] + l - factor * 0.0005, -50 - 1 * N_table["mab"][:200],
             label="buoyancy\nfrequency $N$",
             c="w", lw=2
             )

    # draw sea floor
    x = list(ax2.get_xlim())
    y1 = [0, 0]
    y2 = 2 * list(ax2.get_ylim())[0]
    ax2.axhline(0, c="k")
    ax2.fill_between(x, y1, y2, facecolor="xkcd:charcoal grey", zorder=5)  # , hatch="///")

    """
    #draw second axis with depths
    lons = sorted(list(set(uv_lon)))
    depth_labels = []
    for l in lons:
        local_depth = transect_depth[np.argmin(np.abs(transect_lon - l))]
        depth_labels.append(f"{local_depth:.0f}")


    ax2b = ax2.secondary_xaxis(location="top")
    ax2b = helper.Plot._axstyle(ax2, ticks="out", grid=False)#, spine_offset=40) 
    ax2b.xaxis.set_ticks([])
    ax2b.xaxis.set_ticklabels([])
    ax2b.minorticks_off()
    ax2b.xaxis.set_ticks(lons)
    ax2b.xaxis.set_ticklabels(depth_labels)
    ax2b.set_xlabel("Depth (m)")
    """
    # ax2.set_title("Cross section through the gravity current")
    ax2.legend(loc="upper right", facecolor='k', framealpha=0.8, edgecolor="black", labelcolor="white")

    # fig.tight_layout()


if __name__ == "__main__":

    fig = plt.figure(layout="constrained", figsize=(TWO_COLUMN_WIDTH * cm, 0.8 * TWO_COLUMN_WIDTH * cm))
    gs = GridSpec(nrows=4, ncols=4, figure=fig)
    high_res_proj = ccrs.SouthPolarStereo(central_longitude=-55)
    high_res_proj.threshold = 1e3
    ax0 = fig.add_subplot(gs[:1, :1], projection=high_res_proj)
    ax1 = fig.add_subplot(gs[:1, 1:], projection=ccrs.SouthPolarStereo(central_longitude=-50))
    # identical to ax1 = plt.subplot(gs.new_subplotspec((0, 0), colspan=3))
    ax2 = fig.add_subplot(gs[1:, :])

    # ax0
    overview_map(fig, ax0)

    # ax1
    #joinville_transect(fig,ax1)

    # ax2
    #stratification(fig,ax2)


    # plt.savefig(f"./figures/Intro_figure.png", dpi=800, bbox_inches='tight')
    # plt.savefig(f"./figures/Intro_figure.svg", bbox_inches='tight')
    plt.show()