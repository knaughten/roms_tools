from netCDF4 import Dataset, num2date
from numpy import *

def monthly_avg_roms (file_path, var, shape, month):

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
    time_id = id.variables['ocean_time']
    # Get the year, month, and day (all 1-based) for each output step
    # These are 5-day averages marked with the middle day's date
    time = num2date(time_id[:], units=time_id.units, calendar=time_id.calendar.lower())

    next_month = mod(month+1,12)
    prev_month = mod(month-1,12)

    end_t = -1
    for t in range(size(time)-1, -1, -1):
        if time[t].month-1 == month and time[t].day in range(end_day[month]-2, end_day[month]+1):
            end_t = t
            break
        if time[t].month-1 == next_month and time[t].day in range(start_day[next_month], start_day[next_month]+2):
            end_t = t
            break
    if end_t == -1:
        print 'Error: ' + file_path + ' does not contain a complete ' + month_name[month]
        return

    start_t = -1
    for t in range(end_t, -1, -1):
        if time[t].month-1 == prev_month and time[t].day in range(end_day[prev_month]-1, end_day[prev_month]+1):
            start_t = t
            break
        if time[t].month-1 == month and time[t].day in range(start_day[month], start_day[month]+3):
            start_t = t
            break
    if start_t == -1:
        print 'Error: ' + file_path + ' does not contain a complete ' + month_name[month]
        return

    leap_year = False
    roms_year = time[end_t].year
    if month == 11:
        roms_year = time[start_t].year
    if mod(roms_year, 4) == 0:
        leap_year = True
        if mod(roms_year, 100) == 0:
            leap_year = False
            if mod(roms_year, 400) == 0:
                leap_year = True
    if leap_year:
        end_day[1] = 29

    num_days = 0

    if time[start_t].month-1 == month and time[start_t].day == start_day[month]+2:
        start_days = 5
    elif time[start_t].month-1 == month and time[start_t].day == start_day[month]+1:
        start_days = 4
    elif time[start_t].month-1 == month and time[start_t].day == start_day[month]:
        start_days = 3
    elif time[start_t].month-1 == prev_month and time[start_t].day == end_day[prev_month]:
        start_days = 2
    elif time[start_t].month-1 == prev_month and time[start_t].day == end_day[prev_month]-1:
        start_days = 1
    else:
        print 'Error: starting index is month ' + str(time[start_t].month) + ', day ' + str(time[start_t].day)
        return

    monthly_data += id.variables[var][start_t,:]*start_days
    num_days += start_days

    for t in range(start_t+1, end_t):
        monthly_data += id.variables[var][t,:]*5
        num_days += 5

    if time[end_t].month-1 == next_month and time[end_t].day == start_day[next_month]+1:
        end_days = 1
    elif time[end_t].month-1 == next_month and time[end_t].day == start_day[next_month]:
        end_days = 2
    elif time[end_t].month-1 == month and time[end_t].day == end_day[month]:
        end_days = 3
    elif time[end_t].month-1 == month and time[end_t].day == end_day[month]-1:
        end_days = 4
    elif time[end_t].month-1 == month and time[end_t].day == end_day[month]-2:
        end_days = 5
    else:
        print 'Error: ending index is month ' + str(time[end_t].month) + ', day ' + str(time[end_t].day)
        return

    monthly_data += id.variables[var][end_t,:]*end_days
    num_days += end_days

    if num_days != end_day[month]:
        print 'Error: found ' + str(num_days) + ' days instead of ' + str(end_day[month])
        return

    id.close()
    monthly_data /= num_days

    return monthly_data

    
