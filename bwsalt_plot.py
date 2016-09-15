from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *

# Create a circumpolar plot of bottom water salinity, averaged over the last
# year of simulation.
# Input:
# file_path = path to ocean averages file containing at least one year of
#             5-day averages
# save = optional boolean to save the figure to a file, rather than display it
#        on screen
# fig_name = if save=True, filename for figure
def bwsalt_plot (file_path, save=False, fig_name=None):

    # Degrees to radians conversion factor
    deg2rad = pi/180
    # Centre of missing circle in grid
    lon_c = 50
    lat_c = -83
    # Radius of missing circle (play around with this until it works)
    radius = 10.1
    # Bounds on colour scale
    min_scale = 34
    max_scale = 35

    # Read the grid
    id = Dataset(file_path, 'r')
    lon = id.variables['lon_rho'][:-15,:-2]
    lat = id.variables['lat_rho'][:-15,:-2]
    # Read the last year of bottom water salinity (assume 5-day averages here)
    # and average over time
    bwsalt = mean(id.variables['salt'][-73:,0,:-15,:-2], axis=0)
    id.close()

    # Convert grid to spherical coordinates
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)
    # Find centre in spherical coordinates
    x_c = -(lat_c+90)*cos(lon_c*deg2rad+pi/2)
    y_c = (lat_c+90)*sin(lon_c*deg2rad+pi/2)
    # Build a regular x-y grid and select the missing circle
    x_reg, y_reg = meshgrid(linspace(amin(x), amax(x), num=1000), linspace(amin(y), amax(y), num=1000))
    land_circle = zeros(shape(x_reg))
    land_circle = ma.masked_where(sqrt((x_reg-x_c)**2 + (y_reg-y_c)**2) > radius, land_circle)
    
    # Set colour bounds
    lev = linspace(min_scale, max_scale, num=50)

    # Plot
    fig = figure(figsize=(16,12))
    ax = fig.add_subplot(1,1,1,aspect='equal')
    fig.patch.set_facecolor('white')
    # First shade everything in grey
    contourf(x, y, ones(shape(bwsalt)), 1, colors=(('0.6', '0.6', '0.6')))
    # Fill in the missing circle
    contourf(x_reg, y_reg, land_circle, 1, colors=(('0.6', '0.6', '0.6')))
    # Now shade the salinity (land mask will remain grey)
    contourf(x, y, bwsalt, lev, cmap='jet', extend='both')
    cbar = colorbar(ticks=arange(min_scale, max_scale+0.2, 0.2))
    cbar.ax.tick_params(labelsize=20)
    title(r'Bottom water salinity (psu), annual average', fontsize=30)
    axis('off')

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == '__main__':

    file_path = raw_input("Path to ocean averages file, containing at least 1 year of 5-day averages: ")
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    else:
        save = False
        fig_name = None
    # Make the plot
    bwsalt_plot(file_path, save, fig_name)
    
    
