from numpy import *
from netCDF4 import Dataset
from matplotlib.pyplot import *
from seasonal_avg_cice import *

# Make a figure comparing seasonally-averaged sea ice to ocean salt fluxes,
# from Tamura's dataset to CICE output.
# Input:
# cice_file = path to CICE output file with 5-day averages, containing at least
#             one complete Dec-Nov period (if there are multiple such periods,
#             this script uses the last one)
# cice_year = year that the last Jan-Nov period in CICE corresponds to
# save = optional boolean to save the figure to a file, rather than displaying
#        it on the screen
# fig_name = if save=True, path to the desired filename for figure
def ssflux_tamura_seasonal (cice_file, cice_year, save=False, fig_name=None):

    # Beginning and end of Tamura file paths
    tamura_head = '/short/m68/kaa561/tamura_fluxes/Tamura_ssflux_'
    tamura_tail = '_monthly.nc'
    # Number of days in each month (this is just for Tamura)
    ndays_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    # Season names for titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
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

    # Read seasonally averaged fields
    print 'Reading fresh_ai'
    fresh_ai = seasonal_avg_cice(cice_file, 'fresh_ai', [num_lat, num_lon])
    print 'Reading sss'
    sss = seasonal_avg_cice(cice_file, 'sss', [num_lat, num_lon])
    print 'Reading rain_ai'
    rain_ai = seasonal_avg_cice(cice_file, 'rain_ai', [num_lat, num_lon])
    print 'Reading fsalt_ai'
    fsalt_ai = seasonal_avg_cice(cice_file, 'fsalt_ai', [num_lat, num_lon])
    cice_data_tmp = -1/rho_fw*((fresh_ai-rain_ai)*sss*rho_sw/mps_to_cmpday - fsalt_ai*1e3)
    # Multiply by 1e6 so colour bar is easier to read
    cice_data_tmp *= 1e6
    # Chop off northern boundary
    cice_data_tmp = cice_data_tmp[:,:-15,:]
    # Wrap periodic boundary
    cice_data = ma.empty([size(cice_data_tmp,0), size(cice_data_tmp,1), size(cice_data_tmp,2)+1])
    cice_data[:,:,:-1] = cice_data_tmp
    cice_data[:,:,-1] = cice_data_tmp[:,:,0]

    # Read Tamura grid
    id = Dataset(tamura_head + str(cice_year) + tamura_tail, 'r')
    tamura_lon = id.variables['longitude'][:,:]
    tamura_lat = id.variables['latitude'][:,:]
    id.close()

    # Check for leap years
    leap_year = False
    if mod(cice_year, 4) == 0:
        # Years divisible by 4 are leap years
        leap_year = True
        if mod(cice_year, 100) == 0:
            # Unless they're also divisible by 100, in which case they aren't
            # leap years
            leap_year = False
            if mod(cice_year, 400) == 0:
                # Unless they're also divisible by 400, in which case they are
                # leap years after all
                leap_year = True
    if leap_year:
        # Update last day in February
        ndays_month[1] += 1

    # Set up array for seasonal averages of Tamura data
    tamura_data = ma.empty([4, size(tamura_lon,0), size(tamura_lon,1)])
    tamura_data[:,:,:] = 0.0
    for season in range(4):
        # Work out which months we're interested in
        if season == 0:
            tamura_months = [11, 0, 1]
        elif season == 1:
            tamura_months = [2, 3, 4]
        elif season == 2:
            tamura_months = [5, 6, 7]
        elif season == 3:
            tamura_months = [8, 9, 10]
        season_days = 0

        for month in tamura_months:
            if season == 0 and month == 11:
                # Read December from the previous year
                id = Dataset(tamura_head + str(cice_year-1) + tamura_tail, 'r')
            else:
                # Read the given month from the current year
                id = Dataset(tamura_head + str(cice_year) + tamura_tail, 'r')
            tamura_data_tmp = id.variables['ssflux'][month,:,:]
            id.close()
            # Apply land mask
            tamura_data_tmp = ma.masked_where(isnan(tamura_data_tmp), tamura_data_tmp)
            # Accumulate seasonal average
            tamura_data[season,:,:] += tamura_data_tmp*ndays_month[month]
            season_days += ndays_month[month]
        # Convert from sum to average
        tamura_data[season,:,:] /= season_days

    # Multiply by 1e6 as for CICE
    tamura_data *= 1e6

    # Convert both grids to spherical coordinates
    cice_x = -(cice_lat+90)*cos(cice_lon*deg2rad+pi/2)
    cice_y = (cice_lat+90)*sin(cice_lon*deg2rad+pi/2)
    tamura_x = -(tamura_lat+90)*cos(tamura_lon*deg2rad+pi/2)
    tamura_y = (tamura_lat+90)*sin(tamura_lon*deg2rad+pi/2)

    # Choose colour levels
    lev = linspace(-10, 10, num=50)
    # Bounds for each side of plot
    bdry1 = -35 
    bdry2 = 35
    bdry3 = -33
    bdry4 = 37

    # Plot
    fig = figure(figsize=(20,9))
    # Loop over seasons
    for season in range(4):
        # Tamura
        ax = fig.add_subplot(2, 4, season+1, aspect='equal')
        contourf(tamura_x, tamura_y, tamura_data[season,:,:], lev, cmap='RdYlBu_r', extend='both')
        if season == 0:
            text(-39, 0, 'Tamura', fontsize=24, ha='right')
        title(season_names[season], fontsize=24)
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
        # CICE
        ax = fig.add_subplot(2, 4, season+5, aspect='equal')
        img = contourf(cice_x, cice_y, cice_data[season,:,:], lev, cmap='RdYlBu_r', extend='both')
        if season == 0:
            text(-39, 0, 'CICE', fontsize=24, ha='right')
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
    # Add a horizontal colourbar at the bottom
    cbaxes = fig.add_axes([0.25, 0.04, 0.5, 0.02])
    cbar = colorbar(img, orientation='horizontal', ticks=arange(-10,10+2, 2), cax=cbaxes)
    cbar.ax.tick_params(labelsize=16)
    # Add the main title
    suptitle(r'Ice-to-ocean salt flux (10$^{-6}$ kg/m$^2$/s)', fontsize=30)
    # Make plots closer together
    subplots_adjust(wspace=0.025,hspace=0.025)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    cice_file = raw_input("Path to CICE file, containing at least one complete Dec-Nov period: ")
    cice_year = int(raw_input("Year that the last complete Jan-Nov period in CICE file corresponds to: "))
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    ssflux_tamura_seasonal(cice_file, cice_year, save, fig_name)        
