from numpy import *
from netCDF4 import Dataset, num2date
from matplotlib.pyplot import *

# Make a figure comparing annually-averaged sea ice to ocean salt fluxes,
# from Tamura's dataset to CICE output.
# Input:
# cice_file = path to annually averaged CICE file
# year = year of Tamura data to plot
# save = optional boolean to save the figure to a file, rather than displaying
#        it on the screen
# fig_name = if save=True, path to the desired filename for figure
def ssflux_tamura_annual (cice_file, year, save=False, fig_name=None):

    # Path to Tamura file for this year
    tamura_file = '/short/m68/kaa561/tamura_fluxes/Tamura_ssflux_' + str(year) + '_monthly.nc'
    # Degrees to radians conversion
    deg2rad = pi/180.0
    # Density of freshwater (used by ROMS to convert from kg/m^2/s to psu m/s)
    rho_fw = 1000.0
    # Density of seawater (used by CICE to convert from m/s to kg/m^2/s)
    rho_sw = 1026.0
    # Conversion factor: m/s to cm/day
    mps_to_cmpday = 8.64e6

    # Read the CICE grid
    id = Dataset(cice_file, 'r')
    cice_lon_tmp = id.variables['TLON'][:-15,:]
    cice_lat_tmp = id.variables['TLAT'][:-15,:]
    # Wrap the periodic boundary by 1 cell
    cice_lon = ma.empty([size(cice_lon_tmp,0), size(cice_lon_tmp,1)+1])
    cice_lat = ma.empty([size(cice_lat_tmp,0), size(cice_lat_tmp,1)+1])
    cice_lon[:,:-1] = cice_lon_tmp
    cice_lon[:,-1] = cice_lon_tmp[:,0]
    cice_lat[:,:-1] = cice_lat_tmp
    cice_lat[:,-1] = cice_lat_tmp[:,0]

    # Read the fields we need (already annually averaged)
    fresh_ai = id.variables['fresh_ai'][0,:-15,:]
    sss = id.variables['sss'][0,:-15,:]
    rain_ai = id.variables['rain_ai'][0,:-15,:]
    fsalt_ai = id.variables['fsalt_ai'][0,:-15,:]
    # Convert to units of psu m/s (equivalent to kg/m^2/s of salt)
    # Subtract rain from freshwater flux, since Tamura doesn't count precip
    cice_data_tmp = -1/rho_fw*((fresh_ai-rain_ai)*sss*rho_sw/mps_to_cmpday - fsalt_ai*1e3)
    # Multiply by 1e6 so colour bar is easier to read
    cice_data_tmp *= 1e6
    # Wrap periodic boundary
    cice_data = ma.empty([size(cice_data_tmp,0), size(cice_data_tmp,1)+1])
    cice_data[:,:-1] = cice_data_tmp
    cice_data[:,-1] = cice_data_tmp[:,0]

    # Read Tamura data and grid
    id = Dataset(tamura_file, 'r')
    tamura_lon = id.variables['longitude'][:,:]
    tamura_lat = id.variables['latitude'][:,:]
    tamura_data = mean(id.variables['ssflux'][:,:,:], axis=0)
    id.close()
    # Appply land mask and multiply by 1e6
    tamura_data = ma.masked_where(isnan(tamura_data), tamura_data)
    tamura_data *= 1e6

    # Convert both grids to spherical coordinates
    cice_x = -(cice_lat+90)*cos(cice_lon*deg2rad+pi/2)
    cice_y = (cice_lat+90)*sin(cice_lon*deg2rad+pi/2)
    tamura_x = -(tamura_lat+90)*cos(tamura_lon*deg2rad+pi/2)
    tamura_y = (tamura_lat+90)*sin(tamura_lon*deg2rad+pi/2)

    # Choose colour levels
    bound = 10
    lev = linspace(-bound, bound, num=50)
    # Find boundaries for each side of plot based on extent of grids
    bdry1 = max(amin(cice_x), amin(tamura_x))
    bdry2 = min(amax(cice_x), amax(tamura_x))
    bdry3 = max(amin(cice_y), amin(tamura_y))
    bdry4 = min(amax(cice_y), amax(tamura_y))

    # Plot
    fig = figure(figsize=(20,9))
    # Tamura
    fig.add_subplot(1,2,1, aspect='equal')
    contourf(tamura_x, tamura_y, tamura_data, lev, cmap='RdYlBu_r', extend='both')
    title('Tamura', fontsize=24)
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    # CICE
    fig.add_subplot(1,2,2, aspect='equal')
    img = contourf(cice_x, cice_y, cice_data, lev, cmap='RdYlBu_r', extend='both')
    title('CICE', fontsize=24)
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    # Add a horizontal colourbar at the bottom
    cbaxes = fig.add_axes([0.3, 0.04, 0.4, 0.04])
    cbar = colorbar(img, cax=cbaxes, orientation='horizontal', ticks=arange(-20,20+10,10))
    cbar.ax.tick_params(labelsize=18)
    # Add the main title
    suptitle(r'Ice-to-ocean salt flux (10$^{-6}$ kg/m$^2$/s, annually averaged (' + str(year) + ')', fontsize=30)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    cice_file = raw_input("Path to CICE file, annually averaged: ")
    year = int(raw_input("Year: "))
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    # Make the plot
    ssflux_tamura_annual(cice_file, year, save, fig_name)        
