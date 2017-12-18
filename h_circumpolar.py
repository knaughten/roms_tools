from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.colors import ListedColormap

# Creates a circumpolar Antarctic plot of bathymetry.
# Follows the same process as circumpolar_plot.py, but since h is not
# time-dependent, it requires a special case.
# Input:
# grid_path = path to ROMS grid file
# fig_name = filename for figure
def h_circumpolar (grid_path, fig_name):

    deg2rad = pi/180
    nbdry = -60 + 90

    # Read data
    id = Dataset(grid_path, 'r')
    data = id.variables['h'][:-15,:-1]
    lon = id.variables['lon_rho'][:-15,:-1]
    lat = id.variables['lat_rho'][:-15,:-1]
    mask = id.variables['mask_rho'][:-15,:-1]
    id.close()

    # Mask with land mask
    data = ma.masked_where(mask==0, data)

    # Convert to spherical coordinates
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)

    # Plot
    fig = figure(figsize=(128,96))
    fig.add_subplot(1,1,1, aspect='equal')
    img = pcolor(x, y, data, vmin=0, vmax=2800)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    cbar = colorbar(img, extend='max')
    cbar.ax.tick_params(labelsize=160)
    title('Bathymetry (m)', fontsize=240)
    axis('off')

    show()
    savefig(fig_name)


# Command-line interface
if __name__ == "__main__":

    grid_path = raw_input("Path to ROMS grid file: ")
    fig_name = raw_input("Filename for figure: ")
    h_circumpolar(grid_path, fig_name)    
