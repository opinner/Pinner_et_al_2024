import rioxarray
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.path as mpath
import cmocean
import cartopy.crs as ccrs
cm = 1/2.54  # centimeters in inches
TWO_COLUMN_WIDTH = 12

plt.rcParams.update({
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "font.size": 9
})

def setup_curved_map_shape(central_longitude):
    """
    Sets up a high-resolution South Polar Stereographic projection map with specified latitude and longitude extents.

    This function creates a plot with a South Polar Stereographic projection. It defines the projection, sets the extent
    of the map, and adjusts the boundary of the map to fit within the specified latitude and longitude limits. The resolution
    for computing great circles is increased for finer details. Because I was lazy, this docstring was written by ChatGPT.

    Returns:
        tuple: A tuple containing the figure and axis objects (fig, ax) for the created map.

    Example:
        fig, ax = setup_curved_map_shape()
        ax.coastlines()
        plt.show()
    """
    # Create the plot with a South Polar Stereographic projection
    # Define the projection for a high-resolution South Polar Stereographic map
    central_longitude = np.round(np.mean([min_lon, max_lon]), -1)
    projection = ccrs.SouthPolarStereo(central_longitude=central_longitude)

    # Increase the resolution in the computation of great circles
    # see https://stackoverflow.com/questions/60685245/plot-fine-grained-geodesic-with-cartopy/60724892#comment137539434_60724892
    # if we later run ax.set_boundary() this causes the boundary to be smooth without unwanted corners.
    projection.threshold = 1e3

    #fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection=projection))
    fig, ax = plt.subplots(1, 1,
                           figsize=(0.5 * TWO_COLUMN_WIDTH * cm, 0.5 * TWO_COLUMN_WIDTH * cm,),
                           subplot_kw=dict(projection=projection),
                           layout="constrained"
                           )


    ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())

    # Define the limits for the x-axis and y-axis of the map
    xlim = [min_lon, max_lon]  # longitude
    ylim = [min_lat, max_lat]  # latitude

    # Define a rectangular path for the map extent using matplotlib's Path object
    rect = mpath.Path([[xlim[0], ylim[0]], [xlim[1], ylim[0]], [xlim[1], ylim[1]], [xlim[0], ylim[1]], [xlim[0], ylim[0]]])

    # Transform the rectangular path to the data coordinate system
    proj_to_data = ccrs.PlateCarree()._as_mpl_transform(ax) - ax.transData
    rect_in_target = proj_to_data.transform_path(rect)

    # Set the boundary of the map using the transformed rectangular path
    ax.set_boundary(rect_in_target)

    return fig, ax

# Open the raster data
input_file = '../data/bathymetry/IBCSO_v2_ice-surface_WGS84.tif'
rds = rioxarray.open_rasterio(input_file)

# Define the bounds in lat, lon coordinates for the subset
min_lon = -65
max_lon = -15
min_lat = -80
max_lat = -60

# Clip the raster using the defined bounds
subset = rds.rio.clip_box(minx=min_lon, miny=min_lat, maxx=max_lon, maxy=max_lat)
left, bottom, right, top = subset.rio.bounds()
subset_data = subset[0, :, :]

# Create the plot with a South Polar Stereographic projection
fig, ax = setup_curved_map_shape(central_longitude=-40)
# fig, ax = plt.subplots(1, 1, figsize=(6 * cm, 6 * cm), subplot_kw=dict(projection=ccrs.SouthPolarStereo(central_longitude=-50)))
ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
# Define the limits for the x-axis and y-axis of the map
xlim = [min_lon, max_lon]  # longitude
ylim = [min_lat, max_lat]  # latitude

# Define a rectangular path for the map extent using matplotlib's Path object
# rect = mpath.Path([[xlim[0], ylim[0]], [xlim[1], ylim[0]], [xlim[1], ylim[1]], [xlim[0], ylim[1]], [xlim[0], ylim[0]],])

# Transform the rectangular path to the data coordinate system
# proj_to_data = ccrs.PlateCarree()._as_mpl_transform(ax) - ax.transData
# rect_in_target = proj_to_data.transform_path(rect)
# Set the boundary of the map using the transformed rectangular path
# ax.set_boundary(rect_in_target)

# Define parameters for hill shading
ls = mcolors.LightSource(azdeg=315, altdeg=45)

# Prepare data for hill shading only above 0
shade_only = np.array(subset_data)
shade_only = shade_only.astype(float)  # Ensure float type

# Compute hill shading
shade = ls.hillshade(shade_only, vert_exag=0.1)

# Plot the shaded relief using imshow
ax.imshow(
    shade,
    extent=(left, right, bottom, top),
    transform=ccrs.PlateCarree(),
    origin='upper',
    cmap=cmocean.cm.gray,
    interpolation='none'
)

# Use imshow to plot the bathymetry (negative values)
bathymetry = np.where(subset_data >= 0, np.nan, subset_data)
im = ax.imshow(bathymetry,
               vmin=-5000,
               extent=(left, right, bottom, top),
               transform=ccrs.PlateCarree(),
               origin='upper',
               cmap=cmocean.cm.ice,
               interpolation = 'none'
               )
# cbar = plt.colorbar(im, ax = ax, orientation='vertical')#, shrink = 0.4);
# cbar.set_label('Bathymetry (m)')

# Prepare data for hill shading only below 0 with transparency
shade_only = np.where(subset_data > 0, np.nan, subset_data)
shade_only = shade_only.astype(float)  # Ensure float type

# Compute hill shading with transparency
shade = ls.hillshade(shade_only, vert_exag=0.01, fraction= 0.01)

# Create a colormap that goes from fully transparent to black with only some transparency
c_white = mcolors.colorConverter.to_rgba('white', alpha=0)
c_black = mcolors.colorConverter.to_rgba('black', alpha=0.3)
cmap_rb = mcolors.LinearSegmentedColormap.from_list('rb_cmap', [c_white, c_black], 512)

# Use imshow to shade bathymetry
ax.imshow(shade, extent=(left, right, bottom, top),
          transform=ccrs.PlateCarree(),
          origin='upper',
          cmap=cmap_rb,
          interpolation = 'none'
          )

# Extract coordinates for contour plotting
lon, lat = subset.x.values, subset.y.values
lon_grid, lat_grid = np.meshgrid(lon, lat)

# Create contour lines for specified depths, in the bathymetry map as well as in the colormap
levels = [-4000, -3000, -2000, -1000]
contour = ax.contour(lon_grid, lat_grid, subset_data, transform=ccrs.PlateCarree(),
                     levels=levels, colors="lightgrey", linestyles="solid", linewidths=0.5)
# cbar.add_lines(levels, colors=len(levels)*["w"], linewidths=1)


gl = ax.gridlines(
    xlocs=np.arange(-60,-10,10),
    ylocs=np.arange(-80,-60,4),
    # This is commented out as it does not produce the desired results
    # as not all grid lines will be labelled
    draw_labels=["top", "left", "x", "y"],
    dms=True,
    rotate_labels=False,
    x_inline=False,
    y_inline=False
    )

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

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
fig.patch.set_visible(False)
ax.axis('off')

# Show the plot
plt.savefig(f"./f01b_WeddellSea.svg")
plt.show()
