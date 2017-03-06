from numpy import *
from netCDF4 import Dataset, num2date
from matplotlib.pyplot import *
from rotate_vector_cice import *

def compare_nic_seasonal (cice_file, var_name, colour_bounds=None, save=False, fig_name=None):

    nic_dir_head = '/g/data/gh5/access_om_025-CORE_NYF/output'
    output_number = 137
    nic_dir_tail = '/ice/HISTORY/'
    nic_file_head = 'iceh.0'
    nic_year_number = 133
    max_j = 300

    id = Dataset(cice_file, 'r')
    units = id.variables[var_name].units
    if var_name in ['uvel', 'vvel', 'strairx', 'strairy', 'strocnx', 'strocny']:
        rotate = True
        angle = id.variables['ANGLE'][:-15,:]
        if var_name == 'uvel':
            cmp_flag = 'x'
            other_name = 'vvel'
        elif var_name == 'vvel':
            cmp_flag = 'y'
            other_name = 'uvel'
        elif var_name in ['strairx', 'strocnx']:
            cmp_flag = 'x'
            other_name = var_name.replace('x', 'y')
        elif var_name in ['strairy', 'strocny']:
            cmp_flag = 'y'
            other_name = var_name.replace('x', 'y')            
    else:
        rotate = False
    grid_string = id.variables[var_name].coordinates
    if grid_string.startswith('ULON'):
        lon_name = 'ULON'
        lat_name = 'ULAT'
    else:
        lon_name = 'TLON'
        lat_name = 'TLAT'
    id.close()

    # Starting and ending months (1-based) for each season
    start_month = [12, 3, 6, 9]
    end_month = [2, 5, 8, 11]
    # Starting and ending days of the month (1-based) for each season
    # Assume no leap years, we'll fix this later if we need
    start_day = [1, 1, 1, 1]
    end_day = [28, 31, 31, 30]
    # Number of days in each season (again, ignore leap years for now)
    ndays_season = [90, 92, 92, 91]
    # Number of days in each month (this is just for Nic's output)
    ndays_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    # Season names for titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    # Degrees to radians conversion
    deg2rad = pi/180.0

    # Read CICE grid and time values
    id = Dataset(cice_file, 'r')
    cice_lon_tmp = id.variables[lon_name][:-15,:]
    cice_lat_tmp = id.variables[lat_name][:-15,:]
    # Wrap the periodic boundary by 1 cell
    cice_lon = ma.empty([size(cice_lon_tmp,0), size(cice_lon_tmp,1)+1])
    cice_lat = ma.empty([size(cice_lat_tmp,0), size(cice_lat_tmp,1)+1])
    cice_lon[:,:-1] = cice_lon_tmp
    cice_lon[:,-1] = cice_lon_tmp[:,0]
    cice_lat[:,:-1] = cice_lat_tmp
    cice_lat[:,-1] = cice_lat_tmp[:,0]
    time_id = id.variables['time']
    # Get the year, month, and day (all 1-based) for each output step
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
        # Don't update ndays_month because Nic's setup has no leap years

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

        # Start accumulating data weighted by days
        if rotate:
            this_cmp = id.variables[var_name][start_t_season,:-15,:]
            other_cmp = id.variables[other_name][start_t_season,:-15,:]
            if cmp_flag == 'x':
                data_tmp, other_tmp = rotate_vector_cice(this_cmp, other_cmp, angle)
            elif cmp_flag == 'y':
                other_tmp, data_tmp = rotate_vector_cice(other_cmp, this_cmp, angle)
            cice_data_tmp[season,:,:] += data_tmp*start_days
        else:
            cice_data_tmp[season,:,:] += id.variables[var_name][start_t_season,:-15,:]*start_days
        season_days += start_days

        # Between start_t_season and end_t_season, we want all the days
        for t in range(start_t_season+1, end_t_season):
            if rotate:
                this_cmp = id.variables[var_name][t,:-15,:]
                other_cmp = id.variables[other_name][t,:-15,:]
                if cmp_flag == 'x':
                    data_tmp, other_tmp = rotate_vector_cice(this_cmp, other_cmp, angle)
                elif cmp_flag == 'y':
                    other_tmp, data_tmp = rotate_vector_cice(other_cmp, this_cmp, angle)
                cice_data_tmp[season,:,:] += data_tmp*5
            else:
                cice_data_tmp[season,:,:] += id.variables[var_name][t,:-15,:]*5
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

        if rotate:
            this_cmp = id.variables[var_name][end_t_season,:-15,:]
            other_cmp = id.variables[other_name][end_t_season,:-15,:]
            if cmp_flag == 'x':
                data_tmp, other_tmp = rotate_vector_cice(this_cmp, other_cmp, angle)
            elif cmp_flag == 'y':
                other_tmp, data_tmp = rotate_vector_cice(other_cmp, this_cmp, angle)
            cice_data_tmp[season,:,:] += data_tmp*end_days
        else:
            cice_data_tmp[season,:,:] += id.variables[var_name][end_t_season,:-15,:]*end_days
        season_days += end_days

        # Check that we got the correct number of days        
        if season_days != ndays_season[season]:
            print 'Error: found ' + str(season_days) + ' days instead of ' + str(ndays_season[season])
            return

        # Finished accumulating data, now convert from sum to average
        cice_data_tmp[season,:,:] /= season_days

    # Finished reading all CICE data
    id.close()
    # Wrap the periodic boundary
    cice_data = ma.empty([size(cice_data_tmp,0), size(cice_data_tmp,1), size(cice_data_tmp,2)+1])
    cice_data[:,:,:-1] = cice_data_tmp
    cice_data[:,:,-1] = cice_data_tmp[:,:,0]

    # Read Nic's grid from the January file
    id = Dataset(nic_dir_head + str(output_number) + nic_dir_tail + nic_file_head + str(nic_year_number) + '-01.nc', 'r')
    nic_lon = id.variables[lon_name][:max_j,:]
    nic_lat = id.variables[lat_name][:max_j,:]
    id.close()

    nic_data = ma.empty([4, size(nic_lon,0), size(nic_lon,1)])
    nic_data[:,:,:] = 0.0
    for season in range(4):
        if season == 0:
            months = [12, 1, 2]
        elif season == 1:
            months = [3, 4, 5]
        elif season == 2:
            months = [6, 7, 8]
        elif season == 3:
            months = [9, 10, 11]
        season_days = 0

        for month in months:
            if month == 12:
                filename = nic_dir_head + str(output_number-1) + nic_dir_tail + nic_file_head + str(nic_year_number-1) + '-' + str(month) + '.nc'
            else:
                if month < 10:
                    filename = nic_dir_head + str(output_number) + nic_dir_tail + nic_file_head + str(nic_year_number) + '-0' + str(month) + '.nc'
                else:
                    filename = nic_dir_head + str(output_number) + nic_dir_tail + nic_file_head + str(nic_year_number) + '-' + str(month) + '.nc'
            id = Dataset(filename, 'r')
            nic_data[season,:,:] += id.variables[var_name][0,:max_j,:]*ndays_month[month-1]
            season_days += ndays_month[month-1]
        nic_data[season,:,:] /= season_days

    # Convert both grids to spherical coordinates
    cice_x = -(cice_lat+90)*cos(cice_lon*deg2rad+pi/2)
    cice_y = (cice_lat+90)*sin(cice_lon*deg2rad+pi/2)
    nic_x = -(nic_lat+90)*cos(nic_lon*deg2rad+pi/2)
    nic_y = (nic_lat+90)*sin(nic_lon*deg2rad+pi/2)

    bdry1 = -35
    bdry2 = 39
    bdry3 = -35
    bdry4 = 39

    if colour_bounds is not None:
        lev = linspace(colour_bounds[0], colour_bounds[1], num=50)
        if colour_bounds[0] == -colour_bounds[1]:
            colour_map = 'RdYlBu_r'
        else:
            colour_map = 'jet'
    else:
        if var_name in ['uvel', 'vvel', 'strairx', 'strairy', 'strocnx', 'strocny']:
            max_val = max(amax(abs(nic_data)), amax(abs(cice_data)))
            lev = linspace(-max_val, max_val, num=50)
            colour_map = 'RdYlBu_r'
        else:
            min_val = min(amin(nic_data), amin(cice_data))
            max_val = max(amax(nic_data), amax(cice_data))
            lev = linspace(min_val, max_val, num=50)
            colour_map = 'jet'

    fig = figure(figsize=(20,9))
    for season in range(4):
        ax = fig.add_subplot(2, 4, season+1, aspect='equal')
        contourf(nic_x, nic_y, nic_data[season,:,:], lev, cmap=colour_map, extend='both')
        if season == 0:
            text(-39, 0, 'Nic', fontsize=24, ha='right')
        title(season_names[season], fontsize=24)
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
        ax = fig.add_subplot(2, 4, season+5, aspect='equal')
        img = contourf(cice_x, cice_y, cice_data[season,:,:], lev, cmap=colour_map, extend='both')
        if season == 0:
            text(-39, 0, 'Me', fontsize=24, ha='right')
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
    cbaxes = fig.add_axes([0.25, 0.04, 0.5, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes)
    cbar.ax.tick_params(labelsize=16)
    suptitle(var_name + ' (' + units + ')', fontsize=30)
    subplots_adjust(wspace=0.025,hspace=0.025)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


if __name__ == "__main__":

    cice_file = raw_input("Path to CICE file, containing at least one complete Dec-Nov period: ")
    var_name = raw_input("Variable name: ")
    colour_bounds = None
    get_bounds = raw_input("Set bounds on colour scale (y/n)? ")
    if get_bounds == 'y':
        lower_bound = float(raw_input("Lower bound: "))
        upper_bound = float(raw_input("Upper bound: "))
        colour_bounds = [lower_bound, upper_bound]
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    compare_nic_seasonal(cice_file, var_name, colour_bounds, save, fig_name)

    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                changes = raw_input("Enter a parameter to change: (1) file path, (2) variable name, (3) colour bounds, (4) save/display; or enter to continue: ")
                if len(changes) == 0:
                    break
                else:
                    if int(changes) == 1:
                        cice_file = raw_input("Path to CICE file, containing at least one complete Dec-Nov period: ")
                    elif int(changes) == 2:
                        var_name = raw_input("Variable name: ")
                    elif int(changes) == 3:
                        colour_bounds = None
                        get_bounds = raw_input("Set bounds on colour scale (y/n)? ")
                        if get_bounds == 'y':
                            lower_bound = float(raw_input("Lower bound: "))
                            upper_bound = float(raw_input("Upper bound: "))
                            colour_bounds = [lower_bound, upper_bound]
                    elif int(changes) == 4:
                        save = not save
            if save:
                fig_name = raw_input("File name for figure: ")
            compare_nic_seasonal(cice_file, var_name, colour_bounds, save, fig_name)
        else:
            break


    
