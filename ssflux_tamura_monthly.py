from numpy import *
from netCDF4 import Dataset, num2date
from matplotlib.pyplot import *
from monthly_avg_cice import *

# Make a figure comparing monthly-averaged sea ice to ocean salt fluxes, 
# from Tamura's dataset to CICE output.
# Input:
# cice_file = path to CICE output file with 5-day averages; if it covers more
#             than one instance of the given month, plot the last one
# month = month number (0-indexed) from 0 to 11
# cice_year = year that this month in cice_file corresponds to
# save = optional boolean to save the figure to a file, rather than displaying
#        it on the screen
# fig_name = if save=True, path to the desired filename for figure
def ssflux_tamura_monthly (cice_file, month, cice_year, save=False, fig_name=None):

    # Beginning and end of Tamura file paths
    tamura_head = '/short/m68/kaa561/tamura_fluxes/Tamura_ssflux_'
    tamura_tail = '_monthly.nc'
    # Name of each month, for the title
    month_name = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
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
    num_lon = id.variables['TLON'].shape[1]
    num_lat = id.variables['TLAT'].shape[0]
    id.close()
    # Wrap the periodic boundary by 1 cell
    cice_lon = ma.empty([size(cice_lon_tmp,0), size(cice_lon_tmp,1)+1])
    cice_lat = ma.empty([size(cice_lat_tmp,0), size(cice_lat_tmp,1)+1])
    cice_lon[:,:-1] = cice_lon_tmp
    cice_lon[:,-1] = cice_lon_tmp[:,0]
    cice_lat[:,:-1] = cice_lat_tmp
    cice_lat[:,-1] = cice_lat_tmp[:,0]

    # Get monthly averages of each variable we need
    fresh_ai = monthly_avg_cice(cice_file, 'fresh_ai', [num_lat, num_lon], month)
    sss = monthly_avg_cice(cice_file, 'sss', [num_lat, num_lon], month)
    rain_ai = monthly_avg_cice(cice_file, 'rain_ai', [num_lat, num_lon], month)
    fsalt_ai = monthly_avg_cice(cice_file, 'fsalt_ai', [num_lat, num_lon], month)
    cice_data_tmp = -1/rho_fw*((fresh_ai-rain_ai)*sss*rho_sw/mps_to_cmpday - fsalt_ai*1e3)
    # Chop off northern boundary
    cice_data_tmp = cice_data_tmp[:-15,:]
    # Multiply by 1e6 so colour bar is easier to read
    cice_data_tmp *= 1e6

    # Wrap periodic boundary
    cice_data = ma.empty([size(cice_data_tmp,0), size(cice_data_tmp,1)+1])
    cice_data[:,:-1] = cice_data_tmp
    cice_data[:,-1] = cice_data_tmp[:,0]

    # Read Tamura data and grid
    id = Dataset(tamura_head + str(cice_year) + tamura_tail, 'r')
    tamura_lon = id.variables['longitude'][:,:]
    tamura_lat = id.variables['latitude'][:,:]
    tamura_data = id.variables['ssflux'][month,:,:]
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
    bound = 20.0
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
    contourf(tamura_x, tamura_y, tamura_data, lev, cmap='RdYlBu_r')
    title('Tamura', fontsize=24)
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    # CICE
    fig.add_subplot(1,2,2, aspect='equal')
    img = contourf(cice_x, cice_y, cice_data, lev, cmap='RdYlBu_r')
    title('CICE', fontsize=24)
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    # Add a horizontal colourbar at the bottom
    cbaxes = fig.add_axes([0.3, 0.04, 0.4, 0.04])
    cbar = colorbar(img, orientation='horizontal', ticks=arange(-1,1+0.5, 0.5), cax=cbaxes, extend='both')
    cbar.ax.tick_params(labelsize=18)
    # Add the main title
    suptitle(r'Ice-to-ocean salt flux (10$^{-6}$ kg/m$^2$/s, ' + month_name[month] + ' ' + str(cice_year), fontsize=30)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    cice_file = raw_input("Path to CICE file: ")
    # Convert month from 1-indexed to 0-indexed
    month = int(raw_input("Month number (1-12): "))-1
    cice_year = int(raw_input("Year this month corresponds to in the CICE file: "))
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    # Make the plot
    ssflux_tamura_monthly(cice_file, month, cice_year, save, fig_name)

    while True:
        # Repeat until the user wants to exit
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the parameters until the user is done
                changes = raw_input("Enter a parameter to change: (1) file path, (2) month number, (3) year, (4) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New CICE file
                        cice_file = raw_input("Path to CICE file: ")
                    elif int(changes) == 2:
                        # New month
                        month = int(raw_input("Month number (1-12): "))-1
                    elif int(changes) == 3:
                        cice_year = int(raw_input("Year this month corresponds to in the CICE file: "))
                    elif int(changes) == 4:
                        # Switch from save to display, or vice versa
                        save = not save
            if save:
                # Get a new figure name
                fig_name = raw_input("File name for figure: ")
            # Make the plot
            ssflux_tamura_monthly(cice_file, month, cice_year, save, fig_name)
        else:
            break
            

    

    
        
    
        
