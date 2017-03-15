from numpy import *
from netCDF4 import Dataset
from matplotlib.pyplot import *
from monthly_avg_cice import *

# Make a figure comparing sea ice concentration from NSIDC (1995 data) and
# CICE (latest year of spinup under repeated 1995 forcing), for the given
# month.
# Inupt:
# cice_file = path to CICE output file; if it covers more than once instance of
#             the given month, plot the last one
# month = month number (0-indexed) from 0 to 11
# save = optional boolean to save the figure to a file, rather than displaying
#        it on the screen
# fig_name = if save=True, path to the desired filename for figure
def nsidc_aice_monthly (cice_file, month, save=False, fig_name=None):

    # Strings to construct paths to NSIDC data; other users may need to edit
    nsidc_head = '../nsidc_aice/seaice_conc_monthly_sh'
    nsidc_head_0 = nsidc_head + '_f11_'
    nsidc_head_1 = nsidc_head + '_f13_'
    nsidc_tail = '_v02r00.nc'
    # Name of each month, for the title
    month_name = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
    # Degrees to radians conversion
    deg2rad = pi/180.0

    # Read CICE grid
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

    # Get monthly average of sea ice area
    cice_data_tmp = monthly_avg_cice(cice_file, 'aice', [num_lat, num_lon], month)
    # Chop off the northern boundary
    cice_data_tmp = cice_data_tmp[:-15,:]
    # Wrap the periodic boundary
    cice_data = ma.empty([size(cice_data_tmp,0), size(cice_data_tmp,1)+1])
    cice_data[:,:-1] = cice_data_tmp
    cice_data[:,-1] = cice_data_tmp[:,0]

    # Construct NSIDC file path
    if month+1 < 10:
        nsidc_file = nsidc_head_0 + '19950' + str(month+1) + nsidc_tail
    else:
        nsidc_file = nsidc_head_1 + '1995' + str(month+1) + nsidc_tail

    # Read NSIDC grid and monthly data
    id = Dataset(nsidc_file, 'r')
    nsidc_lon = id.variables['longitude'][:,:]
    nsidc_lat = id.variables['latitude'][:,:]
    nsidc_data_tmp = id.variables['seaice_conc_monthly_cdr'][0,:,:]
    # Read std just for the land mask
    nsidc_mask = id.variables['stdev_of_seaice_conc_monthly_cdr'][0,:,:]
    id.close()

    # Set land mask on NSIDC sea ice concentration
    nsidc_data = ma.empty(shape(nsidc_data_tmp))
    nsidc_data[:,:] = 0.0
    nsidc_data[~nsidc_mask.mask] = nsidc_data_tmp[~nsidc_mask.mask]
    nsidc_data[nsidc_mask.mask] = ma.masked

    # Convert both grids to spherical coordinates
    cice_x = -(cice_lat+90)*cos(cice_lon*deg2rad+pi/2)
    cice_y = (cice_lat+90)*sin(cice_lon*deg2rad+pi/2)
    nsidc_x = -(nsidc_lat+90)*cos(nsidc_lon*deg2rad+pi/2)
    nsidc_y = (nsidc_lat+90)*sin(nsidc_lon*deg2rad+pi/2)

    # Find boundaries for each side of plot based on extent of NSIDC grid
    bdry1 = amax(nsidc_x[:,0])
    bdry2 = amin(nsidc_x[:,-1])
    bdry3 = amin(nsidc_y[:,0])
    bdry4 = amax(nsidc_y[:,-1])

    # Set consistent colour levels
    lev = linspace(0, 1, num=50)

    # Plot
    fig = figure(figsize=(20,9))
    # NSIDC
    fig.add_subplot(1,2,1, aspect='equal')
    contourf(nsidc_x, nsidc_y, nsidc_data, lev)
    title('NSIDC', fontsize=24)
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    # CICE
    fig.add_subplot(1,2,2, aspect='equal')
    img = contourf(cice_x, cice_y, cice_data, lev)
    title('CICE', fontsize=24)
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    # Add a horizontal colorbar at the bottom
    cbaxes = fig.add_axes([0.35, 0.04, 0.3, 0.04])
    cbar = colorbar(img, orientation='horizontal', ticks=arange(0,1+0.25,0.25), cax=cbaxes)
    cbar.ax.tick_params(labelsize=20)
    # Add the main title
    suptitle(month_name[month] + ' sea ice concentration', fontsize=30)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    cice_file = raw_input("Path to CICE file: ")
    # Convert month from 1-indexed to 0-indexed
    month = int(raw_input("Month number (1-12): ")) - 1
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    # Make the plot
    nsidc_aice_monthly(cice_file, month, save, fig_name)

    while True:
        # Repeat until the user wants to exit
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to parameters until the user is done
                changes = raw_input("Enter a parameter to change: (1) file path, (2) month number, (3) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New CICE file
                        cice_file = raw_input("Path to CICE file: ")
                    elif int(changes) == 2:
                        # New month
                        month = int(raw_input("Month number (1-12): ")) - 1
                    elif int(changes) == 3:
                        # Switch from save to display, or vice versa
                        save = not save
            if save:
                # Get a new figure name
                fig_name = raw_input("File name for figure: ")
            # Make the plot
            nsidc_aice_monthly(cice_file, month, save, fig_name)
        else:
            break

        

        

        

        

        
        

    
        

    
        

    
    
