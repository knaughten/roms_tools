from netCDF4 import Dataset, num2date
from numpy import *
from matplotlib.pyplot import *

def ssflux_monthly (roms_file, month, bound=None, save=False, fig_name=None):

    # Month names for plot titles
    month_name = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
    # Number of days per month
    start_day = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    end_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    deg2rad = pi/180

    id = Dataset(roms_file, 'r')
    lon = id.variables['lon_rho'][:-15,:-1]
    lat = id.variables['lat_rho'][:-15,:-1]
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
        print 'Error: ' + roms_file + ' does not contain a complete ' + month_name[month]
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
        print 'Error: ' + roms_file + ' does not contain a complete ' + month_name[month]
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

    ssflux = ma.empty(shape(lon))
    ssflux[:,:] = 0.0
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

    ssflux += id.variables['ssflux'][start_t,:-15,:-1]*start_days
    num_days += start_days

    for t in range(start_t+1, end_t):
        ssflux += id.variables['ssflux'][t,:-15,:-1]*5
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

    ssflux += id.variables['ssflux'][end_t,:-15,:-1]*end_days
    num_days += end_days

    if num_days != end_day[month]:
        print 'Error: found ' + str(num_days) + ' days instead of ' + str(end_day[month])
        return

    id.close()
    ssflux /= num_days
    ssflux *= 1e6

    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)

    if bound is None:
        bound = amax(abs(ssflux))
    lev = linspace(-bound, bound, num=50)
    
    fig = figure(figsize=(16,12))
    fig.add_subplot(1,1,1,aspect='equal')
    img = contourf(x, y, ssflux, lev, cmap='RdYlBu_r', extend='both')
    cbar = colorbar()
    cbar.ax.tick_params(labelsize=20)
    title(month_name[month] + r' surface salinity flux (10$^{-6}$ kg/m$^2$/s)', fontsize=30)
    axis('off')

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


if __name__ == "__main__":

    roms_file = raw_input("Path to ROMS output file: ")
    month = int(raw_input("Month number (1-12): "))-1
    bound = None
    get_bound = raw_input("Set bounds on colour scale (y/n)? ")
    if get_bound == 'y':
        bound = float(raw_input("Maximum absolute value (1e-6 kg/m^2/s): "))
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    ssflux_monthly(roms_file, month, bound, save, fig_name)

    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                changes = raw_input("Enter a parameter to change: (1) file path, (2) month, (3) colour bounds, (4) save/display; or enter to continue: ")
                if len(changes) == 0:
                    break
                else:
                    if int(changes) == 1:
                        roms_file = raw_input("Path to ROMS output file: ")
                    elif int(changes) == 2:
                        month = int(raw_input("Month number (1-12): "))-1
                    elif int(changes) == 3:
                        get_bound = raw_input("Set bounds on colour scale (y/n)? ")
                        if get_bound == 'y':
                            bound = float(raw_input("Maximum absolute value (1e-6 kg/m^2/s): "))
                    elif int(changes) == 4:
                        save = not save
            if save:
                fig_name = raw_input("File name for figure: ")
            ssflux_monthly(roms_file, month, bound, save, fig_name)
        else:
            break
                

    
    



