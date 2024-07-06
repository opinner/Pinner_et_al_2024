import rioxarray
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cmocean as cm
import cartopy.crs as ccrs

# Open the raster data
input_file = '../data/bathymetry/IBCSO_v2_ice-surface_WGS84.tif'
rds = rioxarray.open_rasterio(input_file)

# Define the bounds in lat, lon coordinates for the subset
min_lon = -65
min_lat = -80
max_lon = -40
max_lat = -60

# Clip the raster using the defined bounds
subset = rds.rio.clip_box(minx=min_lon, miny=min_lat, maxx=max_lon, maxy=max_lat)
left, bottom, right, top = subset.rio.bounds()
subset_data = subset[0, :, :]

# Create the plot with a South Polar Stereographic projection
fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection=ccrs.SouthPolarStereo(central_longitude=-50)))
ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())

# Define parameters for hill shading
ls = mcolors.LightSource(azdeg=315, altdeg=45)

# Prepare data for hill shading only above 0
shade_only = np.array(subset_data)
shade_only = shade_only.astype(float)  # Ensure float type

# Compute hill shading
shade = ls.hillshade(shade_only, vert_exag=0.1)

# Plot the shaded relief using imshow
ax.imshow(shade, extent=(left, right, bottom, top), transform=ccrs.PlateCarree(), origin='upper', cmap=cm.cm.gray)

# Use imshow to plot the bathymetry (negative values)
bathymetry = np.where(subset_data >= 0, np.nan, subset_data)
im = ax.imshow(bathymetry, extent=(left, right, bottom, top), transform=ccrs.PlateCarree(), origin='upper', cmap=cm.cm.ice, vmin = -5000)
cbar = plt.colorbar(im, ax = ax, orientation='vertical')#, shrink = 0.4);
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
#ax.imshow(shade, extent=(left, right, bottom, top), transform=ccrs.PlateCarree(), origin='upper', cmap=cmap_rb)

# Extract coordinates for contour plotting
lon, lat = subset.x.values, subset.y.values
lon_grid, lat_grid = np.meshgrid(lon, lat)

# Create contour lines for specified depths, in the bathymetry map as well as in the colormap
levels = [-4000, -3000, -2000, -1000]
contour = ax.contour(lon_grid, lat_grid, subset_data, transform=ccrs.PlateCarree(),
                     levels=levels, colors="white", linestyles="solid", linewidths=0.5)
cbar.add_lines(levels, colors=len(levels)*["w"], linewidths=1)

# Show the plot
plt.show()
