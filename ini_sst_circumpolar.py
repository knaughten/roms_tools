from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *

# Creates a circumpolar Antarctic plot of initial SST from the
# ROMS initialisation file. Follows the same process as
# circumpolar_plot.py, but since the initialisation file is
# set up differently to the ROMS history/averages files, it
# requires a special case.
# Input:
# grid_path = path to ROMS grid file
# file_path = path to ROMS initialisation file
# fig_name = filename for figure
def ini_sst_circumpolar (grid_path, file_path, fig_name):

    # Grid parameters
    theta_s = 4.0
    theta_b = 0.9
    hc = 40
    N = 31
    deg2rad = pi/180

    # Read surface temps
    id = Dataset(file_path, 'r')
    data = id.variables['temp'][0,-1,:-15,:-1]
    id.close()

    # Read lat and lon
    id = Dataset(grid_path, 'r')
    lon = id.variables['lon_rho'][:-15,:-1]
    lat = id.variables['lat_rho'][:-15,:-1]
    mask = id.variables['mask_rho'][:-15,:-1]
    id.close()

    # Mask with land mask
    data = ma.masked_where(mask==0, data)
    lev = linspace(amin(data),amax(data),40)

    # Convert to spherical coordinates
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)

    # Plot
    fig = figure(figsize=(16,12))
    fig.add_subplot(1,1,1, aspect='equal')
    contourf(x, y, data,lev)
    cbar = colorbar()
    cbar.ax.tick_params(labelsize=20)
    title(r'Initial SST ($^{\circ}$C)', fontsize=30)
    axis('off')

    savefig(fig_name)
    #show()


# Command-line interface
if __name__ == "__main__":

    grid_path = raw_input("Path to ROMS grid file: ")
    file_path = raw_input("Path to initial conditions file: ")
    fig_name = raw_input("Filename for figure: ")
    ini_sst_circumpolar(grid_path, file_path, fig_name)        

    
    
