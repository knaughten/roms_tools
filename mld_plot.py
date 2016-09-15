from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *

# Creates a circumpolar Antarctic plot of mixed layer depth in ROMS. Follows
# the same process as circumpolar_plot.py, but masking out the ice shelves
# requires a special case.
# Input:
# file_path = path to ocean history/averages file
# tstep = timestep to plot in file_path
# colour_bounds = optional bounds on colour scale, stored as an array of size
#                 2 with the lower bound first. If colour_bounds = None, then
#                 determine colour scale bounds automatically.
# save = optional boolean flag indicating that the plot should be saved to a 
#        file rather than displayed on the screen
# fig_name = if save=True, filename for figure
def mld_plot (file_path, tstep, colour_bounds, save=False, fig_name=None):

    deg2rad = pi/180

    # Read grid, zice, Hsbl
    id = Dataset(file_path, 'r')
    lon = id.variables['lon_rho'][:-15,:-2]
    lat = id.variables['lat_rho'][:-15,:-2]
    zice = id.variables['zice'][:-15,:-2]
    hsbl = id.variables['Hsbl'][tstep-1,:-15,:-2]
    id.close()
    # Mask out the ice shelves and change the sign
    mld = ma.masked_where(zice!=0, -hsbl)

    # Convert to spherical coordinates
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)

    if colour_bounds is not None:
        # User has set bounds on colour scale
        lev = linspace(colour_bounds[0], colour_bounds[1], num=40)
    else:
        # Determine bounds automatically
        lev = linspace(0, amax(mld), num=40)

    # Plot
    fig = figure(figsize=(16,12))
    fig.add_subplot(1,1,1, aspect='equal')
    contourf(x, y, mld, lev, cmap='jet', extend='both')
    cbar = colorbar()
    cbar.ax.tick_params(labelsize=20)
    title('Mixed layer depth (m)', fontsize=30)
    axis('off')

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to ocean history/averages file: ")
    tstep = int(raw_input("Timestep number (starting at 1): "))
    colour_bounds = None
    get_bounds = raw_input("Set bounds on colour scale (y/n)? ")
    if get_bounds == 'y':
        lower_bound = float(raw_input("Lower bound: "))
        upper_bound = float(raw_input("Upper bound: "))
        colour_bounds = [lower_bound, upper_bound]
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    mld_plot(file_path, tstep, colour_bounds, save, fig_name)
        
        

    
    
