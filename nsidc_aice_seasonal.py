from numpy import *
from netCDF4 import Dataset, num2date
from matplotlib.pyplot import *

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

    # Starting and ending months (1-based) for each season
    start_month = [12, 3, 6, 9]
    end_month = [2, 5, 8, 11]
    # Starting and ending days of the month (1-based) for each season
    # Assume no leap years, we'll fix this later if we need
    start_day = [1, 1, 1, 1]
    end_day = [28, 31, 31, 30]
    # Number of days in each season (again, ignore leap years for now)
    ndays_season = [90, 92, 92, 91]
    # Number of days in each month (no leap years, this is just for NSIDC)
    ndays_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    # Season names for titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    # Degrees to radians conversion
    deg2rad = pi/180.0

    # Read CICE grid and time values
    id = Dataset(cice_file, 'r')
    cice_lon = id.variables['TLON'][:,:]
    cice_lat = id.variables['TLAT'][:,:]
    time_id = id.variables['time']
    # Get the year, month, and day (all 1-based) for each output step
    # These are 5-day averages marked with the last day's date.
    cice_time = num2date(time_id[:], units=time_id.units, calendar=time_id.calendar.lower())

    # Loop backwards through time indices to find the last one we care about
    # (which contains 30 November in its averaging period)
    end_t = -1  # Missing value flag
    for t in range(size(cice_time)-1, -1, -1):
        if cice_time[t].month == end_month[-1] and cice_time[t].day == end_day[-1]:
            end_t = t
            break
        if cice_time[t].month == start_month[0] and cice_time[t].day in range(start_day[0], start_day[0]+4):
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
        if cice_time[t].month == start_month[0] and cice_time[t].day in range(start_day[0], start_day[0]+5):
            start_t = t
            break
    # Make sure we actually found it
    if start_t == -1:
        print 'Error: ' + cice_file + ' does not contain a complete Dec-Nov period'
        return

    # Check if end_t occurs on a leap year
    leap_year = False
    if mod(cice_time[end_t].year, 4) == 0:
        # Years divisible by 4 are leap years
        leap_year = True
        if mod(cice_time[end_t].year, 100) == 0:
            # Unless they're also divisible by 100, in which case they aren't
            # leap years
            leap_year = False
            if mod(cice_time[end_t].year, 400) == 0:
                # Unless they're also divisible by 400, in which case they are
                # leap years after all
                leap_year = True
    if leap_year:
        # Update last day in February
        end_day[0] += 1
        ndays_season[0] += 1

    # Initialise seasonal averages of CICE output
    cice_data = ma.empty([4, size(cice_lon,0), size(cice_lon,1)])
    cice_data[:,:,:] = 0.0
    # Process one season at a time
    for season in range(4):
        season_days = 0  # Number of days in season; this will be incremented
        next_season = mod(season+1, 4)

        # Find starting timestep
        start_t_season = -1
        for t in range(start_t, end_t+1):
            if cice_time[t].month == start_month[season] and cice_time[t].day in range(start_day[season], start_day[season]+5):
                start_t_season = t
                break
        # Make sure we actually found it
        if start_t_season == -1:
            print 'Error: could not find starting timestep for season ' + season_title[season]
            return

        # Find ending timestep
        end_t_season = -1
        for t in range(start_t_season+1, end_t+1):
            if cice_time[t].month == end_month[season] and cice_time[t].day == end_day[season]:
                end_t_season = t
                break
            if cice_time[t].month == start_month[next_season] and cice_time[t].day in range(start_day[next_season], start_day[next_season]+4):
                end_t_season = t
                break
        # Make sure we actually found it
        if end_t_season == -1:
            print 'Error: could not find ending timestep for season ' + season_title[season]
            return

        # Figure out how many of the 5 days averaged in start_t_season are
        # actually within this season
        if cice_time[start_t_season].month == start_month[season] and cice_time[start_t_season].day == start_day[season] + 4:
            # Starting day is in position 1 of 5; we care about all of them
            start_days = 5
        elif cice_time[start_t_season].month == start_month[season] and cice_time[start_t_season].day == start_day[season] + 3:
            # Starting day is in position 2 of 5; we care about the last 4
            start_days = 4
        elif cice_time[start_t_season].month == start_month[season] and cice_time[start_t_season].day == start_day[season]+ 2:
            # Starting day is in position 3 of 5; we care about the last 3
            start_days = 3
        elif cice_time[start_t_season].month == start_month[season] and cice_time[start_t_season].day == start_day[season] + 1:
            # Starting day is in position 4 of 5; we care about the last 2
            start_days = 2
        elif cice_time[start_t_season].month == start_month[season] and cice_time[start_t_season].day == start_day[season]:
            # Starting day is in position 5 of 5; we care about the last 1
            start_days = 1
        else:
            print 'Error for season ' + season_title[season] + ': starting index is month ' + str(cice_time[start_t_season].month) + ', day ' + str(cice_time[start_t_season].day)
            return

        # Start accumulating data weighted by days
        cice_data[season,:,:] += id.variables['aice'][start_t_season,:,:]*start_days
        season_days += start_days

        # Beween start_t_season and end_t_season, we want all the days
        for t in range(start_t_season+1, end_t_season):
            cice_data[season,:,:] += id.variables['aice'][t,:,:]*5
            season_days += 5

        # Figure out how many of the 5 days averaged in end_t_season are
        # actually within this season
        if cice_time[end_t_season].month == start_month[next_season] and cice_time[end_t_season].day == start_day[next_season] + 3:
            # Ending day is in position 1 of 5; we care about the first 1
            end_days = 1
        elif cice_time[end_t_season].month == start_month[next_season] and cice_time[end_t_season].day == start_day[next_season] + 2:
            # Ending day is in position 2 of 5; we care about the first 2
            end_days = 2
        elif cice_time[end_t_season].month == start_month[next_season] and cice_time[end_t_season].day == start_day[next_season] + 1:
            # Ending day is in position 3 of 5; we care about the first 3
            end_days = 3
        elif cice_time[end_t_season].month == start_month[next_season] and cice_time[end_t_season].day == start_day[next_season]:
            # Ending day is in position 4 of 5; we care about the first 4
            end_days = 4
        elif cice_time[end_t_season].month == end_month[season] and cice_time[end_t_season].day == end_day[season]:
            # Ending day is in position 5 of 5; we care about all 5
            end_days = 5
        else:
            print 'Error for season ' + season_title[season] + ': ending index is month ' + str(cice_time[end_t_season].month) + ', day ' + str(cice_time[end_t_season].day)
            return

        cice_data[season,:,:] += id.variables['aice'][end_t_season,:,:]*end_days
        season_days += end_days

        # Check that we got the correct number of days        
        if season_days != ndays_season[season]:
            print 'Error: found ' + str(num_days) + ' days instead of ' + str(ndays_season[season])
            return

        # Finished accumulating data, now convert from sum to average
        cice_data[season,:,:] /= season_days

    # Finished reading all CICE data
    id.close()

    # Read NSIDC grid from the January file
    id = Dataset(nsidc_head_0 + '199501' + nsidc_tail, 'r')
    nsidc_lon = id.variables['longitude'][:,:]
    nsidc_lat = id.variables['latitude'][:,:]
    id.close()

    # Initialise seasonal averages of NSIDC data
    nsidc_data = ma.empty([4, size(nsidc_lon,0), size(nsidc_lon,1)])
    nsidc_data[:,:] = 0.0
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
    fig = figure(figsize=(40,18))
    # Loop over seasons
    for season in range(4):
        # NSIDC
        fig.add_subplot(2, 4, season+1, aspect='equal')
        contourf(nsidc_x, nsidc_y, nsidc_data[season,:,:], lev)
        title('NSIDC (' + season_names[season] + ')', fontsize=24)
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
        # CICE
        fig.add_subplot(2, 4, season+5, aspect='equal')
        img = contourf(cice_x, cice_y, cice_data[season,:,:], lev)
        title('CICE (' + season_names[season] + ')', fontsize=24)
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
    # Add a horizontal colorbar at the bottom
    cbaxes = fig.add_axes([0.4, 0.04, 0.2, 0.04])
    cbar = colorbar(img, orientation='horizontal', ticks=arange(0,1+0.25,0.25), cax=cbaxes)
    cbar.ax.tick_params(labelsize=20)
    # Add the main title
    suptitle('Sea ice concentration', fontsize=30)

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

        
        
        

        

        

        

        
        

    
        

    
        

    
    
