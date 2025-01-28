# Import the necessary libraries
import cartopy.util
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as patches
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np

cm = 1/2.54  # centimeters in inches
TWO_COLUMN_WIDTH = 12

plt.rcParams.update({
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "font.size": 9
})


projection=ccrs.SouthPolarStereo()
projection.threshold = 1e3
fig, ax = plt.subplots(1,
                       subplot_kw=dict(projection=projection),
                       figsize=(2, 2),
                       layout="constrained",
                       )
# Limit the map to -60 degrees latitude and below.
ax.set_extent([-180, 180, -90, -50], ccrs.PlateCarree())

# Draw gridlines on the map
# If labels are drawn, can be chosen here
#ax.gridlines(draw_labels=True)

# Compute a circle in axes coordinates, which we can use as a boundary
# for the map. We can pan/zoom as much as we like - the boundary will be
# permanently circular.
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

# Set the circular boundary for the map using the transformed circle path
ax.set_boundary(circle, transform=ax.transAxes)

# Define the limits for the x-axis and y-axis of the map
xlim = [-65, -15] # longitude
ylim = [-60, -80] # latitude

# Define a rectangular path for the map extent using matplotlib's Path object
rect = mpath.Path([[xlim[0], ylim[0]],
                   [xlim[1], ylim[0]],
                   [xlim[1], ylim[1]],
                   [xlim[0], ylim[1]],
                   [xlim[0], ylim[0]],
                   ])

# Transform the rectangular path to the data coordinate system
proj_to_data = ccrs.PlateCarree()._as_mpl_transform(ax) - ax.transData
rect_in_target = proj_to_data.transform_path(rect)
# Transform the rectangular path to the data coordinate system of the inset map
#proj_to_data = ccrs.PlateCarree()._as_mpl_transform(ax) - ax.transData
#boundary_in_target = proj_to_data.transform_path(boundary)

# Transform the boundary path to a patch object with a width and color
# The zorder is set high to ensure the patch is plotted above the gridlines

Weddell_Sea_patch = patches.PathPatch(rect_in_target, edgecolor='k', facecolor='none', lw=2, zorder=50)

# Turn off axis ticks for the inset map
ax.tick_params(axis='both', which='both', left=False, right=False, bottom=False, top=False,
                   labelbottom=False)

# Add the boundary patch to the inset map
ax.add_patch(Weddell_Sea_patch)

land = cfeature.NaturalEarthFeature(
        category='physical',
        name='land',
        scale='50m',
        #facecolor='white',
        #edgecolor='lightgrey'
    )

ocean = cfeature.NaturalEarthFeature(
        category='physical',
        name='ocean',
        scale='50m',
        #facecolor='white',
        #edgecolor='lightgrey'
    )

# from https://github.com/TomLav/snippets/blob/main/Cartopy%20Ice%20Sheet%20and%20Shelves.ipynb
ice_shelves = cfeature.NaturalEarthFeature(
        category='physical',
        name='antarctic_ice_shelves_polys',
        scale='50m',
        #facecolor='white',
        #edgecolor='lightgrey'
    )

#ax.add_feature(ocean, color="xkcd:sea blue")
ax.set_facecolor("xkcd:sea blue")
ax.add_feature(land, color="xkcd:light grey")
ax.add_feature(ice_shelves, color="white")

# Show the plot
plt.savefig(f"./f01a_SouthernOcean.svg")
plt.show()
