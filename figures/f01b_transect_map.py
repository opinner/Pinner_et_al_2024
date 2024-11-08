import rioxarray
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cmocean as cmocean
import cartopy.crs as ccrs
import src.helper as helper
from src.read_CTDs import load_Joinville_transect_CTDs

ONE_COLUMN_WIDTH = 8.3
TWO_COLUMN_WIDTH = 12
GOLDEN_RATIO = 1.61
cm = 1/2.54  # centimeters in inches

plt.style.use('./paper.mplstyle')

# Open the raster data
input_file = '../data/bathymetry/IBCSO_v2_ice-surface_WGS84.tif'
rds = rioxarray.open_rasterio(input_file)

# Define the bounds in lat, lon coordinates for the subset
min_lon = -56
min_lat = -67
max_lon = -44
max_lat = -61

# Clip the raster using the defined bounds
subset = rds.rio.clip_box(minx=min_lon, miny=min_lat, maxx=max_lon, maxy=max_lat)
left, bottom, right, top = subset.rio.bounds()
subset_data = subset[0, :, :]

# Create the plot with a South Polar Stereographic projection
fig, ax = plt.subplots(1, 1,
                       figsize=(0.5 * TWO_COLUMN_WIDTH * cm +0.5, 0.5 * TWO_COLUMN_WIDTH * cm +0.5),
                       layout="tight",
                       subplot_kw=dict(projection=ccrs.SouthPolarStereo(central_longitude=-50))
                       )
ax.set_extent([-55, -46, -65, -63], crs=ccrs.PlateCarree())

# Define parameters for hill shading
ls = mcolors.LightSource(azdeg=315, altdeg=45)

# Prepare data for hill shading only above 0
shade_only = np.array(subset_data)
shade_only = shade_only.astype(float)  # Ensure float type

# Compute hill shading
shade = ls.hillshade(shade_only, vert_exag=0.1)

# Plot the shaded relief using imshow
ax.imshow(shade, extent=(left, right, bottom, top),
          transform=ccrs.PlateCarree(),
          origin='upper', cmap=cmocean.cm.gray,
          interpolation = 'none'
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

# cbar_fig, cbar_ax = plt.subplots(1)
cbar = plt.colorbar(im, ax = ax, location = "top")#, shrink = 0.4);
cbar.set_label('Bathymetry (m)')

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
ax.imshow(shade, extent=(left, right, bottom, top), transform=ccrs.PlateCarree(), origin='upper', cmap=cmap_rb, interpolation = 'none')

# Extract coordinates for contour plotting
lon, lat = subset.x.values, subset.y.values
lon_grid, lat_grid = np.meshgrid(lon, lat)

# Create contour lines for specified depths, in the bathymetry map as well as in the colormap
levels = [-4000, -3000, -2000, -1000]
contour = ax.contour(lon_grid, lat_grid, subset_data, transform=ccrs.PlateCarree(),
                     levels=levels, colors="lightgrey", linestyles="solid", linewidths=0.5)
cbar.add_lines(levels, colors=len(levels)*["lightgrey"], linewidths=1)


gl = ax.gridlines(xlocs=np.arange(min_lon, max_lon, 2), ylocs=np.arange(min_lat, max_lat, 1),
                  draw_labels={"bottom": "x", "left": "y"},
                  rotate_labels=False,
                  y_inline=False)


markersize = 7

# load all 7 moorings as dataframes
list_of_moorings = helper.IO.load_pickle(name="../data/mooring/list_of_moorings.pkl")
for mooring in list_of_moorings:
    ax.plot(mooring.location.lon, mooring.location.lat, "D", c="tab:red", markersize=markersize,
             transform=ccrs.PlateCarree())

ax.plot(-mooring.location.lon, -mooring.location.lat, "D", c="tab:red", markersize=markersize,
         transform=ccrs.PlateCarree(), label="Moorings")


CTDs = load_Joinville_transect_CTDs()
unique_coords_df = CTDs.drop_duplicates(subset=["Latitude", "Longitude"])

color = "darkgrey"
ax.plot(unique_coords_df["Longitude"], unique_coords_df["Latitude"],
        ".",
        markersize=markersize+1,
        alpha=0.7,
        #markeredgewidth=,
        markeredgecolor=color,
        markerfacecolor='none',
        transform=ccrs.PlateCarree(),
        label="CTD profiles"
        )

ax.legend(fontsize = "small")


# Show the plot
plt.savefig(f"../results/transect_map.svg")
# plt.savefig(f"../results/transect_map.png", dpi = 400, bbox_inches='tight')
plt.show()
