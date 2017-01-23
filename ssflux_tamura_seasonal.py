from numpy import *
from netCDF4 import Dataset, num2date
from matplotlib.pyplot import *

# Make a figure comparing seasonally-averaged sea ice to ocean salt fluxes,
# from Tamura's dataset to CICE output.
# Input:
# cice_file = path to CICE output file with 5-day averages, containing at least
#             one complete Dec-Nov period (if there are multiple such periods,
#             this script uses the last one)
# save = optional boolean to save the figure to a file, rather than displaying
#        it on the screen
# fig_name = if save=True, path to the desired filename for figure
def ssflux_tamura_seasonal (cice_file, save=False, fig_name=None):

    # Beginning and end of Tamura file paths
    tamura_head = '/short/m68/kaa561/tamura_fluxes/Tamura_ssflux_'
    tamura_tail = '_monthly.nc'
    # Starting and ending months (1-based) for each season
    start_month = [12, 3, 6, 9]
    end_month = [2, 5, 8, 11]
    # Starting and ending days of the month (1-based) for each season
    # Assume no leap years, we'll fix this later if we need
    start_day = [1, 1, 1, 1]
    end_day = [28, 31, 31, 30]
    # Number of days in each season (again, ignore leap years for now)
    ndays_season = [90, 92, 92, 91]
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

    # Read the CICE grid and time values
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
    time_id = id.variables['time']
    # Get the year, month, and day (all 1-indexed) for each output step
    # These are 5-day averages marked with the next day's date.
    cice_time = num2date(time_id[:], units=time_id.units, calendar=time_id.calendar.lower())

    # Loop backwards through time indices to find the last one we care about
    # (which contains 30 November in its averaging period)
    end_t = -1  # Missing value flag
    for t in range(size(cice_time)-1, -1, -1):
        if cice_time[t].month == start_month[0] and cice_time[t].day in range(start_day[0], start_day[0]+5):
            end_t = t
            break
    # Make sure we actually found it
    if end_t == -1:
        print 'Error: ' + cice_file + ' does not contain a complete Dec-Nov period'
        return

    # Continue looping backwards to find the first time index we care about
    # (which contains 1 December the previous year in its averaging period)
    start_t = -1  # Missing value flag
    for t in range(end_t-60, -1, -1):
        if cice_time[t].month == start_month[0] and cice_time[t].day in range(start_day[0]+1, start_day[0]+6):
            start_t = t
            break
    # Make sure we actually found it
    if start_t == -1:
        print 'Error: ' + cice_file + ' does not contain a complete Dec-Nov period'
        return

    # Check for leap years
    leap_year = False
    cice_year = cice_time[end_t].year
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
        end_day[0] += 1
        ndays_season[0] += 1
        ndays_month[1] += 1

    # Initialise seasonal averages of CICE output
    cice_data_tmp = ma.empty([4, size(cice_lon_tmp,0), size(cice_lon_tmp,1)])
    cice_data_tmp[:,:,:] = 0.0

    # Process one season at a time
    for season in range(4):
        season_days = 0  # Number of days in season; this will be incremented
        next_season = mod(season+1, 4)

        # Find starting timestep
        start_t_season = -1
        for t in range(start_t, end_t+1):
            if cice_time[t].month == start_month[season] and cice_time[t].day in range(start_day[season]+1, start_day[season]+6):
                start_t_season = t
                break
        # Make sure we actually found it
        if start_t_season == -1:
            print 'Error: could not find starting timestep for season ' + season_names[season]
            return

        # Find ending timestep
        end_t_season = -1
        for t in range(start_t_season+1, end_t+1):
            if cice_time[t].month == start_month[next_season] and cice_time[t].day in range(start_day[next_season], start_day[next_season]+5):
                end_t_season = t
                break
        # Make sure we actually found it
        if end_t_season == -1:
            print 'Error: could not find ending timestep for season ' + season_names[season]
            return

        # Figure out how many of the 5 days averaged in start_t_season are
        # actually within this season
        if cice_time[start_t_season].month == start_month[season] and cice_time[start_t_season].day == start_day[season] + 5:
            # Starting day is in position 1 of 5; we care about all of them
            start_days = 5
        elif cice_time[start_t_season].month == start_month[season] and cice_time[start_t_season].day == start_day[season] + 4:
            # Starting day is in position 2 of 5; we care about the last 4
            start_days = 4
        elif cice_time[start_t_season].month == start_month[season] and cice_time[start_t_season].day == start_day[season]+ 3:
            # Starting day is in position 3 of 5; we care about the last 3
            start_days = 3
        elif cice_time[start_t_season].month == start_month[season] and cice_time[start_t_season].day == start_day[season] + 2:
            # Starting day is in position 4 of 5; we care about the last 2
            start_days = 2
        elif cice_time[start_t_season].month == start_month[season] and cice_time[start_t_season].day == start_day[season] + 1:
            # Starting day is in position 5 of 5; we care about the last 1
            start_days = 1
        else:
            print 'Error for season ' + season_names[season] + ': starting index is month ' + str(cice_time[start_t_season].month) + ', day ' + str(cice_time[start_t_season].day)
            return

        # Read the fields we need at start_t_season
        fresh_ai = id.variables['fresh_ai'][start_t_season,:-15,:]
        sss = id.variables['sss'][start_t_season,:-15,:]
        rain_ai = id.variables['rain_ai'][start_t_season,:-15,:]
        fsalt_ai = id.variables['fsalt_ai'][start_t_season,:-15,:]
        # Start accumulating data weighted by days
        # Convert to units of psu m/s (equivalent to kg/m^2/s of salt)
        # Subtract rain from freshwater flux, since Tamura doesn't count precip
        cice_data_tmp[season,:,:] += -1/rho_fw*((fresh_ai-rain_ai)*sss*rho_sw/mps_to_cmpday - fsalt_ai*1e3)*start_days
        season_days += start_days

        # Between start_t_season and end_t_season, we want all the days
        for t in range(start_t_season+1, end_t_season):
            fresh_ai = id.variables['fresh_ai'][t,:-15,:]
            sss = id.variables['sss'][t,:-15,:]
            rain_ai = id.variables['rain_ai'][t,:-15,:]
            fsalt_ai = id.variables['fsalt_ai'][t,:-15,:]
            cice_data_tmp[season,:,:] += -1/rho_fw*((fresh_ai-rain_ai)*sss*rho_sw/mps_to_cmpday - fsalt_ai*1e3)*5
            season_days += 5

        # Figure out how many of the 5 days averaged in end_t_season are
        # actually within this season
        if cice_time[end_t_season].month == start_month[next_season] and cice_time[end_t_season].day == start_day[next_season] + 4:
            # Ending day is in position 1 of 5; we care about the first 1
            end_days = 1
        elif cice_time[end_t_season].month == start_month[next_season] and cice_time[end_t_season].day == start_day[next_season] + 3:
            # Ending day is in position 2 of 5; we care about the first 2
            end_days = 2
        elif cice_time[end_t_season].month == start_month[next_season] and cice_time[end_t_season].day == start_day[next_season] + 2:
            # Ending day is in position 3 of 5; we care about the first 3
            end_days = 3
        elif cice_time[end_t_season].month == start_month[next_season] and cice_time[end_t_season].day == start_day[next_season] + 1:
            # Ending day is in position 4 of 5; we care about the first 4
            end_days = 4
        elif cice_time[end_t_season].month == start_month[next_season] and cice_time[end_t_season].day == start_day[next_season]:
            # Ending day is in position 5 of 5; we care about all 5
            end_days = 5
        else:
            print 'Error for season ' + season_names[season] + ': ending index is month ' + str(cice_time[end_t_season].month) + ', day ' + str(cice_time[end_t_season].day)
            return

        fresh_ai = id.variables['fresh_ai'][end_t_season,:-15,:]
        sss = id.variables['sss'][end_t_season,:-15,:]
        rain_ai = id.variables['rain_ai'][end_t_season,:-15,:]
        fsalt_ai = id.variables['fsalt_ai'][end_t_season,:-15,:]
        cice_data_tmp[season,:,:] += -1/rho_fw*((fresh_ai-rain_ai)*sss*rho_sw/mps_to_cmpday - fsalt_ai*1e3)*end_days
        season_days += end_days

        # Check that we got the correct number of days
        if season_days != ndays_season[season]:
            print 'Error: found ' + str(season_days) + ' days instead of ' + str(ndays_season[season])
            return

        # Finished accumulating data, now convert from sum to average
        cice_data_tmp[season,:,:] /= season_days

    id.close()
    # Multiply by 1e6 so colour bar is easier to read
    cice_data_tmp *= 1e6

    # Wrap periodic boundary
    cice_data = ma.empty([size(cice_data_tmp,0), size(cice_data_tmp,1), size(cice_data_tmp,2)+1])
    cice_data[:,:,:-1] = cice_data_tmp
    cice_data[:,:,-1] = cice_data_tmp[:,:,0]

    # Read Tamura grid
    id = Dataset(tamura_head + str(cice_year) + tamura_tail, 'r')
    tamura_lon = id.variables['longitude'][:,:]
    tamura_lat = id.variables['latitude'][:,:]
    id.close()

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
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    ssflux_tamura_seasonal(cice_file, save, fig_name)        
