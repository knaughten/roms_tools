from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *

# Creates a circumpolar Antarctic plot of ice shelf draft.
# Input:
# grid_path = path to ROMS grid file
# fig_name = filename for figure
def zice_circumpolar (grid_path, fig_name):

    # Grid parameters
    theta_s = 0.9
    theta_b = 4.0
    hc = 40
    N = 31
    deg2rad = pi/180
    # Northern boundary 63S for plot
    nbdry = -63+90
    # Centre of missing circle in grid
    lon_c = 50
    lat_c = -83
    # Radius of missing circle (play around with this until it works)
    radius = 10.1

    # Read data
    id = Dataset(grid_path, 'r')
    zice = -1*id.variables['zice'][:-15,:-2]
    lon = id.variables['lon_rho'][:-15,:-2]
    lat = id.variables['lat_rho'][:-15,:-2]
    mask_rho = id.variables['mask_rho'][:-15,:-2]
    id.close()

    # Mask the open ocean and land
    zice = ma.masked_where(zice==0, zice)
    # Get land/zice mask
    open_ocn = copy(mask_rho)
    open_ocn[zice!=0] = 0
    land_zice = ma.masked_where(open_ocn==1, open_ocn)

    # Convert to spherical coordinates
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)
    # Find centre in spherical coordinates
    x_c = -(lat_c+90)*cos(lon_c*deg2rad+pi/2)
    y_c = (lat_c+90)*sin(lon_c*deg2rad+pi/2)
    # Build a regular x-y grid and select the missing circle
    x_reg, y_reg = meshgrid(linspace(-nbdry, nbdry, num=1000), linspace(-nbdry, nbdry, num=1000))
    land_circle = zeros(shape(x_reg))
    land_circle = ma.masked_where(sqrt((x_reg-x_c)**2 + (y_reg-y_c)**2) > radius, land_circle)

    lev = linspace(0, amax(zice), num=50)

    # Plot
    fig = figure(figsize=(16,12))
    ax = fig.add_subplot(1,1,1, aspect='equal')
    fig.patch.set_facecolor('white')
    # First shade land and zice in grey (include zice so there are no white
    # patches near the grounding line where contours meet)
    contourf(x, y, land_zice, 1, colors=(('0.6', '0.6', '0.6')))
    # Fill in the missing circle
    contourf(x_reg, y_reg, land_circle, 1, colors=(('0.6', '0.6', '0.6')))
    # Now shade zice
    contourf(x, y, zice, lev)
    cbar = colorbar(ticks=arange(0,amax(zice),250))
    cbar.ax.tick_params(labelsize=20)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    title('Ice shelf draft (m)', fontsize=30)
    axis('off')

    show()
    #fig.savefig(fig_name)


# Command-line interface
if __name__ == "__main__":

    grid_path = raw_input("Path to ROMS grid file: ")
    fig_name = raw_input("Filename for figure: ")
    zice_circumpolar(grid_path, fig_name)    
