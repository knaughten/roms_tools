from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *

def cice_transport (file_path, save=False, fig_name=None):

    deg2rad = pi/180

    id = Dataset(file_path,'r')
    dice = mean(id.variables['frazil'][-73:,:-15,:] - id.variables['meltb'][-73:,:-15,:], axis=0)
    lon = id.variables['ULON'][:-15,:]
    lat = id.variables['ULAT'][:-15,:]
    id.close()

    # Convert to spherical coordinates
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)

    lev = linspace(-1, 1, num=40)

    # Plot
    fig = figure(figsize=(16,12))
    fig.add_subplot(1,1,1, aspect='equal')
    contourf(x, y, dice, lev, cmap='RdYlBu_r', extend='both')
    cbar = colorbar()
    cbar.ax.tick_params(labelsize=20)
    title('Freezing minus melting (cm/day), annually averaged', fontsize=30)
    axis('off')

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()    


if __name__ == "__main__":

    file_path = raw_input("Path to CICE history file, containing at least one year of 5-day averages: ")
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    cice_transport(file_path, save, fig_name)
    
