from netCDF4 import Dataset, num2date
from numpy import *

def monthly_avg_cice (file_path, var, shape, month):

    # Starting and ending days in each month
    # Assume no leap years, we'll fix this later if we need
    start_day = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    end_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    month_name = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']

    # Set up array to hold monthly data
    monthly_data = ma.empty(shape)
    monthly_data[:] = 0.0

    # Read the time values
    id = Dataset(file_path, 'r')
    time_id = id.variables['time']
    # Get the year, month, and day (all 1-indexed) for each output step
    # These are 5-day averages marked with the next day's date.
    cice_time = num2date(time_id[:], units=time_id.units, calendar='standard')

    next_month = mod(month+1,12)

    end_t = -1
    for t in range(size(cice_time)-1, -1, -1):
        if cice_time[t].month-1 == next_month and cice_time[t].day in range(start_day[next_month], start_day[next_month]+5):
            end_t = t
            break
    if end_t == -1:
        print 'Error: ' + cice_file + ' does not contain a complete ' + month_name[month]
        return

    start_t = -1
    for t in range(end_t, -1, -1):
        if cice_time[t].month-1 == month and cice_time[t].day in range(start_day[month]+1, start_day[month]+6):
            start_t = t
            break
    if start_t == -1:
        print 'Error: ' + cice_file + ' does not contain a complete ' + month_name[month]
        return

    leap_year = False
    cice_year = cice_time[end_t].year
    if month == 11:
        cice_year = cice_time[start_t].year
    if mod(cice_year, 4) == 0:
        leap_year = True
        if mod(cice_year, 100) == 0:
            leap_year = False
            if mod(cice_year, 400) == 0:
                leap_year = True
    if leap_year:
        end_day[1] = 29

    num_days = 0

    # Figure out how many of the 5 days averaged in start_t are actually within
    # this month
    if cice_time[start_t].month-1 == month and cice_time[start_t].day == start_day[month]+5:
        # Starting day is in position 1 of 5; we care about all of them
        start_days = 5
    elif cice_time[start_t].month-1 == month and cice_time[start_t].day == start_day[month]+4:
        # Starting day is in position 2 of 5; we care about the last 4
        start_days = 4
    elif cice_time[start_t].month-1 == month and cice_time[start_t].day == start_day[month]+3:
        # Starting day is in position 3 of 5; we care about the last 3
        start_days = 3
    elif cice_time[start_t].month-1 == month and cice_time[start_t].day == start_day[month]+2:
        # Starting day is in position 4 of 5; we care about the last 2
        start_days = 2
    elif cice_time[start_t].month-1 == month and cice_time[start_t].day == start_day[month]+1:
        # Starting day is in position 5 of 5; we care about the last 1
        start_days = 1
    else:
        print 'Error: starting index is month ' + str(cice_time[start_t].month) + ', day ' + str(cice_time[start_t].day)
        return

    monthly_data += id.variables[var][start_t,:]*start_days
    num_days += start_days

    for t in range(start_t+1, end_t):
        monthly_data += id.variables[var][t,:]*5
        num_days +=5

    # Figure out how many of the 5 days averaged in end_t are actually within
    # this month
    if cice_time[end_t].month-1 == next_month and cice_time[end_t].day == start_day[next_month]+4:
        # Ending day is in position 1 of 5; we care about the first 1
        end_days = 1
    elif cice_time[end_t].month-1 == next_month and cice_time[end_t].day == start_day[next_month]+3:
        # Ending day is in position 2 of 5; we care about the first 2
        end_days = 2
    elif cice_time[end_t].month-1 == next_month and cice_time[end_t].day == start_day[next_month]+2:
        # Ending day is in position 3 of 5; we care about the first 3
        end_days = 3
    elif cice_time[end_t].month-1 == next_month and cice_time[end_t].day == start_day[next_month]+1:
        # Ending day is in position 4 of 5; we care about the first 4
        end_days = 4
    elif cice_time[end_t].month-1 == next_month and cice_time[end_t].day == start_day[next_month]:
        # Ending day is in position 5 of 5; we care about all 5
        end_days = 5
    else:
        print 'Error: ending index is month ' + str(cice_time[end_t].month) + ', day ' + str(cice_time[end_t].day)
        return

    monthly_data += id.variables[var][end_t,:]*end_days
    num_days += end_days

    if num_days != end_day[month]:
        print 'Error: found ' + str(num_days) + ' days instead of ' + str(end_day[month])
        return

    id.close()
    monthly_data /= num_days

    return monthly_data
