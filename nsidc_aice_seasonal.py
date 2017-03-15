from numpy import *
from netCDF4 import Dataset
from matplotlib.pyplot import *
from seasonal_avg_cice import *

# Make a figure comparing sea ice concentration from NSIDC (1995 data) and CICE
# (latest year of spinup under repeated 1995 forcing) for each season.
# Input:
# cice_file = path to CICE output file with 5-day averages, containing at least
#             one complete Dec-Nov period (if there are multiple such periods, 
#             this script uses the last one)
# save = optional boolean to save the figure to a file, rather than displaying
#        it on the screen
# fig_name = if save=True, path to the desired filename for figure
def nsidc_aice_seasonal (cice_file, save=False, fig_name=None):

    nsidc_head = '../nsidc_aice/seaice_conc_monthly_sh'
    nsidc_head_0 = nsidc_head + '_f11_'
    nsidc_head_1 = nsidc_head + '_f13_'
    nsidc_tail = '_v02r00.nc'
    # Number of days in each month (this is just for NSIDC)
    ndays_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    # Season names for titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    # Degrees to radians conversion
    deg2rad = pi/180.0

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

    # Read seasonally averaged fields
    cice_data_tmp = seasonal_avg_cice(cice_file, 'aice', [num_lat, num_lon])
    # Chop off northern boundary
    cice_data_tmp = cice_data_tmp[:,:-15,:]
    # Wrap the periodic boundary
    cice_data = ma.empty([size(cice_data_tmp,0), size(cice_data_tmp,1), size(cice_data_tmp,2)+1])
    cice_data[:,:,:-1] = cice_data_tmp
    cice_data[:,:,-1] = cice_data_tmp[:,:,0]

    # Read NSIDC grid from the January file
    id = Dataset(nsidc_head_0 + '199501' + nsidc_tail, 'r')
    nsidc_lon = id.variables['longitude'][:,:]
    nsidc_lat = id.variables['latitude'][:,:]
    id.close()

    # Initialise seasonal averages of NSIDC data
    nsidc_data = ma.empty([4, size(nsidc_lon,0), size(nsidc_lon,1)])
    nsidc_data[:,:,:] = 0.0
    # Process one season at a time
    for season in range(4):
        # Figure out which months we care about
        if season == 0:
            # DJF
            nsidc_months = [12, 1, 2]
        elif season == 1:
            # MAM
            nsidc_months = [3, 4, 5]
        elif season == 2:
            # JJA
            nsidc_months = [6, 7, 8]
        elif season == 3:
            # SON
            nsidc_months = [9, 10, 11]
        season_days = 0 # Number of days in season; this will be incremented

        # Process one month at a time
        for month in nsidc_months:
            # Construct NSIDC file path
            if month < 10:
                nsidc_file = nsidc_head_0 + '19950' + str(month) + nsidc_tail
            else:
                nsidc_file = nsidc_head_1 + '1995' + str(month) + nsidc_tail
            # Read concentration data
            id = Dataset(nsidc_file, 'r')
            nsidc_data_raw = id.variables['seaice_conc_monthly_cdr'][0,:,:]
            # Read std just for the mask
            nsidc_mask = id.variables['stdev_of_seaice_conc_monthly_cdr'][0,:,:]
            id.close()
            # Set land mask
            nsidc_data_tmp = ma.empty(shape(nsidc_data_raw))
            nsidc_data_tmp[:,:] = 0.0
            nsidc_data_tmp[~nsidc_mask.mask] = nsidc_data_raw[~nsidc_mask.mask]
            nsidc_data_tmp[nsidc_mask.mask] = ma.masked
            # Accumulate master array, weighted with number of days per month
            nsidc_data[season,:,:] += nsidc_data_tmp*ndays_month[month-1]
            season_days += ndays_month[month-1]

        # Convert from sum to average
        nsidc_data[season,:,:] /= season_days

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
    # Loop over seasons
    for season in range(4):
        # NSIDC
        ax = fig.add_subplot(2, 4, season+1, aspect='equal')
        contourf(nsidc_x, nsidc_y, nsidc_data[season,:,:], lev)
        if season == 0:
            text(-39, 0, 'NSIDC', fontsize=24, ha='right')
        title(season_names[season], fontsize=24)
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
        # CICE
        ax = fig.add_subplot(2, 4, season+5, aspect='equal')
        img = contourf(cice_x, cice_y, cice_data[season,:,:], lev)
        if season == 0:
            text(-39, 0, 'CICE', fontsize=24, ha='right')
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
    # Add a horizontal colorbar at the bottom
    cbaxes = fig.add_axes([0.25, 0.04, 0.5, 0.02])
    cbar = colorbar(img, orientation='horizontal', ticks=arange(0,1+0.25,0.25), cax=cbaxes)
    cbar.ax.tick_params(labelsize=16)
    # Add the main title
    suptitle('Sea ice concentration', fontsize=30)
    # Decrease space between plots
    subplots_adjust(wspace=0.025,hspace=0.025)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    cice_file = raw_input("Path to CICE file, containing at least one complete Dec-Nov period: ")
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    nsidc_aice_seasonal(cice_file, save, fig_name)

        
        
        

        

        

        

        
        

    
        

    
        

    
    
