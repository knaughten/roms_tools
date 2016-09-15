from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from cartesian_grid_2d import *

# Make a circumpolar plot of the horizontal resolution of the ROMS grid (square
# root of the area of each cell). 
# Input:
# grid_path = path to ROMS grid file
# save = optional boolean indicating to save the figure, rather than display it
#        on the screen
# fig_name = if save=True, filename to save the figure
def grid_res (grid_path, save=False, fig_name=None):

    # Degrees to radians conversion factor
    deg2rad = pi/180

    # Read grid
    id = Dataset(grid_path, 'r')
    lon = id.variables['lon_rho'][:-15,:-2]
    lat = id.variables['lat_rho'][:-15,:-2]
    mask = id.variables['mask_rho'][:-15,:-2]
    id.close()

    # Get differentials
    dx, dy = cartesian_grid_2d(lon, lat)
    # Calculate resolution: square root of the area, converted to km
    res = sqrt(dx*dy)*1e-3
    # Apply land mask
    res = ma.masked_where(mask==0, res)

    # Polar coordinates for plotting
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)

    # Colour levels
    lev = linspace(0, 40, num=50)

    # Plot
    fig = figure(figsize=(16,12))
    fig.add_subplot(1,1,1, aspect='equal')
    contourf(x, y, res, lev)
    cbar = colorbar()
    cbar.ax.tick_params(labelsize=20)
    title('Grid resolution (km)', fontsize=30)
    axis('off')

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    grid_path = raw_input("Path to ROMS grid file: ")
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("Filename for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    grid_res(grid_path, save, fig_name)

    

    
