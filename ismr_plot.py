from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *

# Make a circumpolar plot of ice shelf melt rates averaged over the last year
# of simulation.
# Input:
# file_path = path to ROMS ocean averages file, containing at least 1 year of
#             5-day averaged output
# save = optional boolean to save the figure to a file, rather than display it
#        on screen
# fig_name = if save=True, filename for figure
def ismr_plot (file_path, save=False, fig_name=None):

    # Degrees to radians conversion factor
    deg2rad = pi/180
    # Northern boundary 63S for plot
    nbdry = -63+90
    # Centre of missing circle in grid
    lon_c = 50
    lat_c = -83
    # Radius of missing circle (play around with this until it works)
    radius = 10.0

    # Read the grid
    id = Dataset(file_path, 'r')
    lon = id.variables['lon_rho'][:,:]
    lat = id.variables['lat_rho'][:,:]
    mask_rho = id.variables['mask_rho'][:,:]
    zice = id.variables['zice'][:,:]
    # Read the last year of ice shelf melt rates (assume 5-day averages here),
    # average over time, and convert from m/s to m/y
    ismr = mean(id.variables['m'][-73:,:,:], axis=0)*60*60*24*365.25
    id.close()
    # Mask the open ocean and land out of the melt rates
    ismr = ma.masked_where(zice==0, ismr)

    # Get land/zice mask
    open_ocn = copy(mask_rho)
    open_ocn[zice!=0] = 0
    land_zice = ma.masked_where(open_ocn==1, open_ocn)

    # Convert grid to spherical coordinates
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)
    # Find centre in spherical coordinates
    x_c = -(lat_c+90)*cos(lon_c*deg2rad+pi/2)
    y_c = (lat_c+90)*sin(lon_c*deg2rad+pi/2)
    # Build a regular x-y grid and select the missing circle
    x_reg, y_reg = meshgrid(linspace(-nbdry, nbdry, num=1000), linspace(-nbdry, nbdry, num=1000))
    land_circle = zeros(shape(x_reg))
    land_circle = ma.masked_where(sqrt((x_reg-x_c)**2 + (y_reg-y_c)**2) > radius, land_circle)

    # Set bounds on colour scale
    lev = linspace(-10, 10, num=40)

    # Set up plot
    fig = figure(figsize=(16,12))
    ax = fig.add_subplot(1,1,1, aspect='equal')
    fig.patch.set_facecolor('white')
    # First shade land and zice in grey (include zice so there are no white
    # patches near the grounding line where contours meet)
    contourf(x, y, land_zice, 1, colors=(('0.6', '0.6', '0.6')))
    # Fill in the missing circle
    contourf(x_reg, y_reg, land_circle, 1, colors=(('0.6', '0.6', '0.6')))
    # Now shade the melt rate
    contourf(x, y, ismr, lev, cmap='RdBu_r', extend='both')
    cbar = colorbar(ticks=arange(-10, 10+2, 2))
    cbar.ax.tick_params(labelsize=20)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    title('Ice shelf melt rate (m/y), annual average', fontsize=30)
    axis('off')

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()

        
# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to ocean averages file, containing at least 1 year of 5-day averages: ")
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    # Make the plot
    ismr_plot(file_path, save, fig_name)
    
