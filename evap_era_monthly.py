from numpy import *
from netCDF4 import Dataset, num2date
from matplotlib.pyplot import *

def evap_era_monthly (roms_file, month, colour_bounds=None, save=False, fig_name=None):

    era_head = '/short/y99/kaa561/CMIP5_forcing/atmos/climatology/ERA_Interim_monthly/ER_'
    era_tail = '_monthly_orig.nc'
    # Starting and ending days in each month
    # Assume no leap years, we'll fix this later if we need
    start_day = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    end_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    # Name of each month, for the title
    month_name = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
    # Degrees to radians conversion
    deg2rad = pi/180.0
    # Conversion factor: m/12h to kg/m^2/2, opposite signs
    era_conv = -1000./(12.*60.*60.)

    id = Dataset(roms_file, 'r')
    roms_lon = id.variables['lon_rho'][:-15,:-1]
    roms_lat = id.variables['lat_rho'][:-15,:-1]
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

    roms_evap = ma.empty(shape(roms_lon))
    roms_evap[:,:] = 0.0
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

    roms_evap += id.variables['evaporation'][start_t,:-15,:-1]*start_days
    num_days += start_days

    for t in range(start_t+1, end_t):
        roms_evap += id.variables['evaporation'][t,:-15,:-1]*5
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

    roms_evap += id.variables['evaporation'][end_t,:-15,:-1]*end_days
    num_days += end_days

    if num_days != end_day[month]:
        print 'Error: found ' + str(num_days) + ' days instead of ' + str(end_day[month])
        return

    id.close()
    roms_evap /= num_days
    roms_evap *= 1e6

    id = Dataset(era_head + str(roms_year) + era_tail, 'r')
    era_lon_1d = id.variables['longitude'][:]
    era_lat_1d = id.variables['latitude'][:]
    era_lon, era_lat = meshgrid(era_lon_1d, era_lat_1d)
    era_evap = id.variables['e'][month,:,:]*era_conv
    id.close()
    era_evap *= 1e6

    roms_x = -(roms_lat+90)*cos(roms_lon*deg2rad+pi/2)
    roms_y = (roms_lat+90)*sin(roms_lon*deg2rad+pi/2)
    era_x = -(era_lat+90)*cos(era_lon*deg2rad+pi/2)
    era_y = (era_lat+90)*sin(era_lon*deg2rad+pi/2)

    if colour_bounds is not None:
        lev = linspace(colour_bounds[0], colour_bounds[1], num=50)
    else:
        min_bound = min(amin(roms_evap), amin(era_evap))
        max_bound = max(amax(roms_evap), amax(era_evap))
        lev = linspace(min_bound, max_bound, num=50)
    bdry = -30+90

    fig = figure(figsize=(20,9))
    fig.add_subplot(1,2,1, aspect='equal')
    contourf(era_x, era_y, era_evap, lev, extend='both')
    title('ERA-Interim', fontsize=24)
    xlim([-bdry, bdry])
    ylim([-bdry, bdry])
    axis('off')
    fig.add_subplot(1,2,2, aspect='equal')
    img = contourf(roms_x, roms_y, roms_evap, lev, extend='both')
    title('ROMS', fontsize=24)
    xlim([-bdry, bdry])
    ylim([-bdry, bdry])
    axis('off')
    cbaxes = fig.add_axes([0.3, 0.04, 0.4, 0.04])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes)
    suptitle(r'Evaporation (10$^{-6}$ kg/m$^2$/s), ' + month_name[month] + ' ' + str(roms_year), fontsize=30)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


if __name__ == "__main__":

    roms_file = raw_input("Path to ROMS file: ")
    month = int(raw_input("Month number (1-12): "))-1
    colour_bounds = None
    get_bounds = raw_input("Set bounds on colour scale (y/n)? ")
    if get_bounds == 'y':
        lower_bound = float(raw_input("Lower bound: "))
        upper_bound = float(raw_input("Upper bound: "))
        colour_bounds = [lower_bound, upper_bound]
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    evap_era_monthly(roms_file, month, colour_bounds, save, fig_name)

    while True:
        # Repeat until the user wants to exit
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the parameters until the user is done
                changes = raw_input("Enter a parameter to change: (1) file path, (2) month number, (3) colour bounds, (4) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New ROMS file
                        roms_file = raw_input("Path to ROMS file: ")
                    elif int(changes) == 2:
                        # New month
                        month = int(raw_input("Month number (1-12): "))-1
                    elif int(changes) == 3:
                        colour_bounds = None
                        get_bounds = raw_input("Set bounds on colour scale (y/n)? ")
                        if get_bounds == 'y':
                            lower_bound = float(raw_input("Lower bound: "))
                            upper_bound = float(raw_input("Upper bound: "))
                            colour_bounds = [lower_bound, upper_bound]
                    elif int(changes) == 4:
                        # Switch from save to display, or vice versa
                        save = not save
            if save:
                # Get a new figure name
                fig_name = raw_input("File name for figure: ")
            # Make the plot
            evap_era_monthly(roms_file, month, colour_bounds, save, fig_name)
        else:
            break
