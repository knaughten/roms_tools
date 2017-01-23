from numpy import *
from netCDF4 import Dataset, num2date
from matplotlib.pyplot import *

# Creates a 4x2 plot of seasonally averaged sea ice concentration (top row) and
# thickness (bottom row) over the last year of simulation.
# Input:
# cice_file = path to CICE output file with 5-day averages, containing at least
#             one complete Dec-Nov period (if there are multiple such periods,
#             this script uses the last one)
# save = optional boolean to save the figure to a file, rather than displaying
#        it on the screen
# fig_name = if save=True, path to the desired filename for figure
def aice_hi_seasonal (cice_file, save=False, fig_name=None):

    # Starting and ending months (1-based) for each season
    start_month = [12, 3, 6, 9]
    end_month = [2, 5, 8, 11]
    # Starting and ending days of the month (1-based) for each season
    # Assume no leap years, we'll fix this later if we need
    start_day = [1, 1, 1, 1]
    end_day = [28, 31, 31, 30]
    # Number of days in each season (again, ignore leap years for now)
    ndays_season = [90, 92, 92, 91]
    # Season names for titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    # Degrees to radians conversion
    deg2rad = pi/180.0

    # Read CICE grid and time values
    id = Dataset(cice_file, 'r')
    lon_tmp = id.variables['TLON'][:-15,:]
    lat_tmp = id.variables['TLAT'][:-15,:]
    # Wrap the periodic boundary by 1 cell
    lon = ma.empty([size(lon_tmp,0), size(lon_tmp,1)+1])
    lat = ma.empty([size(lat_tmp,0), size(lat_tmp,1)+1])
    lon[:,:-1] = lon_tmp
    lon[:,-1] = lon_tmp[:,0]
    lat[:,:-1] = lat_tmp
    lat[:,-1] = lat_tmp[:,0]
    time_id = id.variables['time']
    # Get the year, month, and day (all 1-based) for each output step
    # These are 5-day averages marked with the next day's date.
    time = num2date(time_id[:], units=time_id.units, calendar=time_id.calendar.lower())

    # Loop backwards through time indices to find the last one we care about
    # (which contains 30 November in its averaging period)
    end_t = -1  # Missing value flag
    for t in range(size(time)-1, -1, -1): 
        if time[t].month == start_month[0] and time[t].day in range(start_day[0], start_day[0]+5):
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
        if time[t].month == start_month[0] and time[t].day in range(start_day[0]+1, start_day[0]+6):
            start_t = t
            break
    # Make sure we actually found it
    if start_t == -1:
        print 'Error: ' + cice_file + ' does not contain a complete Dec-Nov period'
        return

    # Check for leap years
    leap_year = False
    if mod(time[end_t].year, 4) == 0:
        # Years divisible by 4 are leap years
        leap_year = True
        if mod(time[end_t].year, 100) == 0:
            # Unless they're also divisible by 100, in which case they aren't
            # leap years
            leap_year = False
            if mod(time[end_t].year, 400) == 0:
                # Unless they're also divisible by 400, in which case they are
                # leap years after all
                leap_year = True
    if leap_year:
        # Update last day in February
        end_day[0] += 1
        ndays_season[0] += 1

    # Initialise seasonal averages of CICE output
    aice_tmp = ma.empty([4, size(lon_tmp,0), size(lon_tmp,1)])
    aice_tmp[:,:,:] = 0.0
    hi_tmp = ma.empty([4, size(lon_tmp,0), size(lon_tmp,1)])
    hi_tmp[:,:,:] = 0.0
    # Process one season at a time
    for season in range(4):
        season_days = 0  # Number of days in season; this will be incremented
        next_season = mod(season+1, 4)

        # Find starting timestep
        start_t_season = -1
        for t in range(start_t, end_t+1):
            if time[t].month == start_month[season] and time[t].day in range(start_day[season]+1, start_day[season]+6):
                start_t_season = t
                break
        # Make sure we actually found it
        if start_t_season == -1:
            print 'Error: could not find starting timestep for season ' + season_names[season]
            return

        # Find ending timestep
        end_t_season = -1
        for t in range(start_t_season+1, end_t+1):
            if time[t].month == start_month[next_season] and time[t].day in range(start_day[next_season], start_day[next_season]+5):
                end_t_season = t
                break
        # Make sure we actually found it
        if end_t_season == -1:
            print 'Error: could not find ending timestep for season ' + season_names[season]
            return

        # Figure out how many of the 5 days averaged in start_t_season are
        # actually within this season
        if time[start_t_season].month == start_month[season] and time[start_t_season].day == start_day[season] + 5:
            # Starting day is in position 1 of 5; we care about all of them
            start_days = 5
        elif time[start_t_season].month == start_month[season] and time[start_t_season].day == start_day[season] + 4:
            # Starting day is in position 2 of 5; we care about the last 4
            start_days = 4
        elif time[start_t_season].month == start_month[season] and time[start_t_season].day == start_day[season]+ 3:
            # Starting day is in position 3 of 5; we care about the last 3
            start_days = 3
        elif time[start_t_season].month == start_month[season] and time[start_t_season].day == start_day[season] + 2:
            # Starting day is in position 4 of 5; we care about the last 2
            start_days = 2
        elif time[start_t_season].month == start_month[season] and time[start_t_season].day == start_day[season] + 1:
            # Starting day is in position 5 of 5; we care about the last 1
            start_days = 1
        else:
            print 'Error for season ' + season_names[season] + ': starting index is month ' + str(time[start_t_season].month) + ', day ' + str(time[start_t_season].day)
            return

        # Start accumulating data weighted by days
        aice_tmp[season,:,:] += id.variables['aice'][start_t_season,:-15,:]*start_days
        hi_tmp[season,:,:] += id.variables['hi'][start_t_season,:-15,:]*start_days
        season_days += start_days

        # Between start_t_season and end_t_season, we want all the days
        for t in range(start_t_season+1, end_t_season):
            aice_tmp[season,:,:] += id.variables['aice'][t,:-15,:]*5
            hi_tmp[season,:,:] += id.variables['hi'][t,:-15,:]*5
            season_days += 5

        # Figure out how many of the 5 days averaged in end_t_season are
        # actually within this season
        if time[end_t_season].month == start_month[next_season] and time[end_t_season].day == start_day[next_season] + 4:
            # Ending day is in position 1 of 5; we care about the first 1
            end_days = 1
        elif time[end_t_season].month == start_month[next_season] and time[end_t_season].day == start_day[next_season] + 3:
            # Ending day is in position 2 of 5; we care about the first 2
            end_days = 2
        elif time[end_t_season].month == start_month[next_season] and time[end_t_season].day == start_day[next_season] + 2:
            # Ending day is in position 3 of 5; we care about the first 3
            end_days = 3
        elif time[end_t_season].month == start_month[next_season] and time[end_t_season].day == start_day[next_season] + 1:
            # Ending day is in position 4 of 5; we care about the first 4
            end_days = 4
        elif time[end_t_season].month == start_month[next_season] and time[end_t_season].day == start_day[next_season]:
            # Ending day is in position 5 of 5; we care about all 5
            end_days = 5
        else:
            print 'Error for season ' + season_names[season] + ': ending index is month ' + str(time[end_t_season].month) + ', day ' + str(time[end_t_season].day)
            return

        aice_tmp[season,:,:] += id.variables['aice'][end_t_season,:-15,:]*end_days
        hi_tmp[season,:,:] += id.variables['hi'][end_t_season,:-15,:]*end_days
        season_days += end_days

        # Check that we got the correct number of days   
        if season_days != ndays_season[season]:
            print 'Error: found ' + str(season_days) + ' days instead of ' + str(ndays_season[season])
            return

        # Finished accumulating data, now convert from sum to average
        aice_tmp[season,:,:] /= season_days
        hi_tmp[season,:,:] /= season_days

    # Finished reading all CICE data
    id.close()

    # Wrap the periodic boundary
    aice = ma.empty([size(aice_tmp,0), size(aice_tmp,1), size(aice_tmp,2)+1])
    aice[:,:,:-1] = aice_tmp
    aice[:,:,-1] = aice_tmp[:,:,0]
    hi = ma.empty([size(hi_tmp,0), size(hi_tmp,1), size(hi_tmp,2)+1])
    hi[:,:,:-1] = hi_tmp
    hi[:,:,-1] = hi_tmp[:,:,0]

    # Get circumpolar x and y coordinates for plotting
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)

    # Set boundaries of plot
    bdry1 = -35
    bdry2 = 35
    bdry3 = -33
    bdry4 = 37

    # Set consistent colour levels
    lev1 = linspace(0, 1, num=50)
    lev2 = linspace(0, 2.5, num=50)

    # Plot
    fig = figure(figsize=(20,9))
    # Loop over seasons
    for season in range(4):
        # Concentration
        ax = fig.add_subplot(2, 4, season+1, aspect='equal')
        img = contourf(x, y, aice[season,:,:], lev1)
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
        if season == 0:
            text(-39, 0, 'aice (1)', fontsize=21, ha='right')
        title(season_names[season], fontsize=24)
        if season == 3:
            # Add colorbar
            cbaxes1 = fig.add_axes([0.92, 0.55, 0.01, 0.3])
            cbar1 = colorbar(img, ticks=arange(0,1+0.25,0.25), cax=cbaxes1)
            cbar1.ax.tick_params(labelsize=16)
        # Thickness
        ax = fig.add_subplot(2, 4, season+5, aspect='equal')
        img = contourf(x, y, hi[season,:,:], lev2, extend='both')
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
        if season == 0:
            text(-39, 0, 'hi (m)', fontsize=21, ha='right')
        if season == 3:
            # Add colorbar
            cbaxes2 = fig.add_axes([0.92, 0.15, 0.01, 0.3])
            cbar2 = colorbar(img, ticks=arange(0,2.5+0.5,0.5), cax=cbaxes2)
            cbar2.ax.tick_params(labelsize=16)
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
    aice_hi_seasonal(cice_file, save, fig_name)

        
        
        

        

        

        

        
        

    
        

    
        

    
    
