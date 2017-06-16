from netCDF4 import Dataset, num2date
from numpy import *

# Calculate seasonal averages (DJF, MAM, JJA, SON) of the given variable in the
# given CICE file.
# Input:
# file_path = path to CICE output file containing 5-day averages, including at
#             least one complete December-November period. If there are multiple
#             such instances the last one will be used.
# var = variable name
# shape = vector containing the dimensions (excluding time) of the variable
# Output:
# seasonal_data = array of data averaged over each season, dimension 4 x shape
def seasonal_avg_cice (file_path, var, shape):

    # Starting and ending months (1-based) for each season
    start_month = [12, 3, 6, 9]
    end_month = [2, 5, 8, 11]
    # Starting and ending days of the month (1-based) for each season
    # Assume no leap years, we'll fix this later if we need
    start_day = [1, 1, 1, 1]
    end_day = [28, 31, 31, 30]
    # Number of days in each season (again, ignore leap years for now)
    ndays_season = [90, 92, 92, 91]
    season_names = ['DJF', 'MAM', 'JJA', 'SON']

    # Set up array to hold seasonal data
    new_shape = [4] + shape
    seasonal_data = ma.empty(new_shape)
    seasonal_data[:] = 0.0

    # Read the time values
    id = Dataset(file_path, 'r')
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

        # Start accumulating data weighted by days
        seasonal_data[season,:] += id.variables[var][start_t_season,:]*start_days
        season_days += start_days

        # Between start_t_season and end_t_season, we want all the days
        for t in range(start_t_season+1, end_t_season):
            seasonal_data[season,:] += id.variables[var][t,:]*5
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

        seasonal_data[season,:] += id.variables[var][end_t_season,:]*end_days
        season_days += end_days

        # Check that we got the correct number of days
        if season_days != ndays_season[season]:
            print 'Error: found ' + str(season_days) + ' days instead of ' + str(ndays_season[season])
            return

        # Finished accumulating data, now convert from sum to average
        seasonal_data[season,:,:] /= season_days

    id.close()

    return seasonal_data
