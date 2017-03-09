from numpy import *
from netCDF4 import Dataset, num2date
from matplotlib.pyplot import *

def cice_evap_era_monthly (cice_file, month, colour_bounds=None, save=False, fig_name=None):

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
    # Conversion factor: m/12h to cm/day
    era_conv = 100*2

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

    cice_evap_tmp = ma.empty(shape(cice_lon_tmp))
    cice_evap_tmp[:,:] = 0.0
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

    cice_evap_tmp += id.variables['evap_ai'][start_t,:-15,:]*start_days
    num_days += start_days

    for t in range(start_t+1, end_t):
        cice_evap_tmp += id.variables['evap_ai'][t,:-15,:]*5
        num_days += 5

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

    cice_evap_tmp += id.variables['evap_ai'][end_t,:-15,:]*end_days
    num_days += end_days

    if num_days != end_day[month]:
        print 'Error: found ' + str(num_days) + ' days instead of ' + str(end_day[month])
        return

    id.close()    
    cice_evap_tmp = -1*cice_evap_tmp/num_days

    cice_evap = ma.empty([size(cice_evap_tmp,0), size(cice_evap_tmp,1)+1])
    cice_evap[:,:-1] = cice_evap_tmp
    cice_evap[:,-1] = cice_evap_tmp[:,0]

    id = Dataset(era_head + str(cice_year) + era_tail, 'r')
    era_lon_1d = id.variables['longitude'][:]
    era_lat_1d = id.variables['latitude'][:]
    era_lon, era_lat = meshgrid(era_lon_1d, era_lat_1d)
    era_evap = id.variables['e'][month,:,:]*era_conv
    id.close()
    era_evap = -1*era_evap

    cice_x = -(cice_lat+90)*cos(cice_lon*deg2rad+pi/2)
    cice_y = (cice_lat+90)*sin(cice_lon*deg2rad+pi/2)
    era_x = -(era_lat+90)*cos(era_lon*deg2rad+pi/2)
    era_y = (era_lat+90)*sin(era_lon*deg2rad+pi/2)

    if colour_bounds is not None:
        lev = linspace(colour_bounds[0], colour_bounds[1], num=50)
    else:
        min_bound = min(amin(cice_evap), amin(era_evap))
        max_bound = max(amax(cice_evap), amax(era_evap))
        lev = linspace(min_bound, max_bound, num=50)
    bdry = -50+90

    fig = figure(figsize=(20,9))
    fig.add_subplot(1,2,1, aspect='equal')
    contourf(era_x, era_y, era_evap, lev, extend='both')
    title('ERA-Interim', fontsize=24)
    xlim([-bdry, bdry])
    ylim([-bdry, bdry])
    axis('off')
    fig.add_subplot(1,2,2, aspect='equal')
    img = contourf(cice_x, cice_y, cice_evap, lev, extend='both')
    title('CICE', fontsize=24)
    xlim([-bdry, bdry])
    ylim([-bdry, bdry])
    axis('off')
    cbaxes = fig.add_axes([0.3, 0.04, 0.4, 0.04])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes)
    suptitle(r'Evaporation (cm/day), ' + month_name[month] + ' ' + str(cice_year), fontsize=30)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


if __name__ == "__main__":

    cice_file = raw_input("Path to CICE file: ")
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
    cice_evap_era_monthly(cice_file, month, colour_bounds, save, fig_name)

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
                        # New CICE file
                        cice_file = raw_input("Path to CICE file: ")
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
            cice_evap_era_monthly(cice_file, month, colour_bounds, save, fig_name)
        else:
            break
