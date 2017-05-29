from netCDF4 import Dataset, num2date
from numpy import *

# Average the given variable in the given CICE file over the given month.
# Input:
# file_path = path to CICE output file containing 5-day averages, including at
#             least one complete instance of the given month.
# var = variable name
# shape = vector containing the dimensions (excluding time) of the variable
# month = month to average over (0-11)
# instance = optional integer indicating which instance of the given month in
#            this file we should use. For instance=1 use the first instance,
#            etc. If instance=-1 (the default) the last instance is used.
# Output:
# monthly_data = array of data averaged over the given month.
def monthly_avg_cice (file_path, var, shape, month, instance=-1):

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

    if instance == -1:
        # Default: find the last instance of this month
        end_t = -1
        for t in range(size(cice_time)-1, -1, -1):
            check_leapyear(cice_time[t].year, end_day)
            if cice_time[t].month-1 == next_month and cice_time[t].day in range(start_day[next_month], start_day[next_month]+5):
                end_t = t
                break
        if end_t == -1:
            print 'Error: ' + cice_file + ' does not contain a complete ' + month_name[month]
            return
        start_t = -1
        for t in range(end_t, -1, -1):
            check_leapyear(cice_time[t].year, end_day)
            if cice_time[t].month-1 == month and cice_time[t].day in range(start_day[month]+1, start_day[month]+6):
                start_t = t
                break
        if start_t == -1:
            print 'Error: ' + cice_file + ' does not contain a complete ' + month_name[month]
            return
    else:
        # Find the given instance of the month
        count = 0
        start_t = -1
        for t in range(size(cice_time)):
            check_leapyear(cice_time[t].year, end_day)
            if cice_time[t].month-1 == month and cice_time[t].day in range(start_day[month]+1, start_day[month]+6):
                count += 1
                if count == instance:
                    start_t = t
                    break
        if start_t == -1:
            print 'Error: ' + file_path + ' does not contain ' + str(instance) + ' ' + month_name[month] + 's'
            return
        end_t = -1
        for t in range(start_t+1, size(cice_time)):
            check_leapyear(cice_time[t].year, end_day)
            if cice_time[t].month-1 == next_month and cice_time[t].day in range(start_day[next_month], start_day[next_month]+5):
                end_t = t
                break
        if end_t == -1:
            print 'Error: ' + file_path + ' does not contain  ' + str(instance) + ' ' + month_name[month] + 's'
            return

    # Figure out what year it is, based on the year of the last timestep we care
    # about
    cice_year = cice_time[end_t].year
    if month == 11:
        # December: the last timestep might be January next year
        # Use the year of the first timestep we care about
        cice_year = cice_time[start_t].year
    # Check for leap years
    check_leapyear(cice_year, end_day)

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

    # Integrate data weighted by start_days
    monthly_data += id.variables[var][start_t,:]*start_days
    num_days += start_days

    # Between start_t and end_t, we care about all the days
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

    # Integrate data weighted by end_days
    monthly_data += id.variables[var][end_t,:]*end_days
    num_days += end_days

    # Make sure we got the right number of days
    if num_days != end_day[month]:
        print 'Error: found ' + str(num_days) + ' days instead of ' + str(end_day[month])
        return

    id.close()
    # Convert from integral to average
    monthly_data /= num_days

    return monthly_data


def check_leapyear (year, end_day):

    leap_year = False
    if mod(year, 4) == 0:
        leap_year = True
        if mod(year, 100) == 0:
            leap_year = False
            if mod(year, 400) == 0:
                leap_year = True
    if leap_year:
        end_day[1] = 29
    else:
        end_day[1] = 28
