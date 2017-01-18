from numpy import *
from netCDF4 import Dataset, num2date
from matplotlib.pyplot import *

# Make a figure comparing sea ice concentration from NSIDC (1995 data) and
# CICE (latest year of spinup under repeated 1995 forcing), for the given
# month.
# Inupt:
# cice_file = path to CICE output file; if it covers more than once instance of
#             the given month, plot the last one
# month = month number (0-indexed) from 0 to 11
# save = optional boolean to save the figure to a file, rather than displaying
#        it on the screen
# fig_name = if save=True, path to the desired filename for figure
def nsidc_aice_monthly (cice_file, month, save=False, fig_name=None):

    # Strings to construct paths to NSIDC data; other users may need to edit
    nsidc_head = '../nsidc_aice/seaice_conc_monthly_sh'
    nsidc_head_0 = nsidc_head + '_f11_'
    nsidc_head_1 = nsidc_head + '_f13_'
    nsidc_tail = '_v02r00.nc'
    # Starting and ending days in each month
    # Assume no leap years, we'll fix this later if we need
    start_day = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    end_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    # Name of each month, for the title
    month_name = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
    # Degrees to radians conversion
    deg2rad = pi/180.0

    # Read CICE grid and time values
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
    # Get the year, month, and day (all 1-based) for each output step
    # These are 5-day averages marked with the next day's date.
    cice_time = num2date(time_id[:], units=time_id.units, calendar=time_id.calendar.lower())

    # Loop backwards through time indices to find the last one we care about
    # (which contains the last day of the month in its averaging period)
    end_t = -1  # Missing value flag
    for t in range(size(cice_time)-1, -1, -1):
        # Note that cice_time[t].month is converted to 0-indexed
        if cice_time[t].month-1 == next_month and cice_time[t].day in range(start_day[next_month], start_day[next_month]+5):
            end_t = t
            break
    # Make sure we actually found it
    if end_t == -1:
        print 'Error: ' + cice_file + ' does not contain a complete ' + month_name[month]
        return

    # Continue looping backwards to find the first time index we care about
    start_t = -1  # Missing value flag
    for t in range(end_t, -1, -1):
        if cice_time[t].month-1 == month and cice_time[t].day in range(start_day[month]+1, start_day[month]+6):
            start_t = t
            break
    # Make sure we actually found it
    if start_t == -1:
        print 'Error: ' + cice_file + ' does not contain a complete ' + month_name[month]
        return

    # Check if this is a leap year
    leap_year = False
    cice_year = cice_time[end_t].year
    if month == 11:
        cice_year = cice_time[start_t].year
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
        end_day[1] = 29

    # Calculate monthly average of CICE output
    cice_data_tmp = ma.empty(shape(cice_lon_tmp))
    cice_data_tmp[:,:] = 0.0
    num_days = 0
    
    # Figure out how many of the 5 days averaged in start_t are actually within
    # this month
    if cice_time[start_t].month-1 == month and cice_time[start_t].day == start_day[month] + 5:
        # Starting day is in position 1 of 5; we care about all of them
        start_days = 5
    elif cice_time[start_t].month-1 == month and cice_time[start_t].day == start_day[month] + 4:
        # Starting day is in position 2 of 5; we care about the last 4
        start_days = 4
    elif cice_time[start_t].month-1 == month and cice_time[start_t].day == start_day[month]+ 3:
        # Starting day is in position 3 of 5; we care about the last 3
        start_days = 3
    elif cice_time[start_t].month-1 == month and cice_time[start_t].day == start_day[month] + 2:
        # Starting day is in position 4 of 5; we care about the last 2
        start_days = 2
    elif cice_time[start_t].month-1 == month and cice_time[start_t].day == start_day[month] + 1:
        # Starting day is in position 5 of 5; we care about the last 1
        start_days = 1
    else:
        print 'Error: starting index is month ' + str(cice_time[start_t].month) + ', day ' + str(cice_time[start_t].day)
        return

    # Start accumulating data weighted by days
    cice_data_tmp[:,:] += id.variables['aice'][start_t,:-15,:]*start_days
    num_days += start_days

    # Beween start_t and end_t, we want all the days
    for t in range(start_t+1, end_t):
        cice_data_tmp[:,:] += id.variables['aice'][t,:-15,:]*5
        num_days += 5

    # Figure out how many of the 5 days averaged in end_t are actually within
    # this month
    next_month = mod(month+1, 12)
    if cice_time[end_t].month-1 == next_month and cice_time[end_t].day == start_day[next_month] + 4:
        # Ending day is in position 1 of 5; we care about the first 1
        end_days = 1
    elif cice_time[end_t].month-1 == next_month and cice_time[end_t].day == start_day[next_month] + 3:
        # Ending day is in position 2 of 5; we care about the first 2
        end_days = 2
    elif cice_time[end_t].month-1 == next_month and cice_time[end_t].day == start_day[next_month] + 2:
        # Ending day is in position 3 of 5; we care about the first 3
        end_days = 3
    elif cice_time[end_t].month-1 == next_month and cice_time[end_t].day == start_day[next_month] + 1:
        # Ending day is in position 4 of 5; we care about the first 4
        end_days = 4
    elif cice_time[end_t].month-1 == next_month and cice_time[end_t].day == start_day[next_month]:
        # Ending day is in position 5 of 5; we care about all 5
        end_days = 5
    else:
        print 'Error: ending index is month ' + str(cice_time[end_t].month) + ', day ' + str(cice_time[end_t].day)
        return

    cice_data_tmp[:,:] += id.variables['aice'][end_t,:-15,:]*end_days
    num_days += end_days

    # Check that we got the correct number of days
    if num_days != end_day[month]:
        print 'Error: found ' + str(num_days) + ' days instead of ' + str(end_day[month])
        return

    # Finished accumulating data
    id.close()
    # Convert from sum to average
    cice_data_tmp[:,:] /= num_days

    # Wrap the periodic boundary
    cice_data = ma.empty([size(cice_data_tmp,0), size(cice_data_tmp,1)+1])
    cice_data[:,:-1] = cice_data_tmp
    cice_data[:,-1] = cice_data_tmp[:,0]

    # Construct NSIDC file path
    if month+1 < 10:
        nsidc_file = nsidc_head_0 + '19950' + str(month+1) + nsidc_tail
    else:
        nsidc_file = nsidc_head_1 + '1995' + str(month+1) + nsidc_tail

    # Read NSIDC grid and monthly data
    id = Dataset(nsidc_file, 'r')
    nsidc_lon = id.variables['longitude'][:,:]
    nsidc_lat = id.variables['latitude'][:,:]
    nsidc_data_tmp = id.variables['seaice_conc_monthly_cdr'][0,:,:]
    # Read std just for the land mask
    nsidc_mask = id.variables['stdev_of_seaice_conc_monthly_cdr'][0,:,:]
    id.close()

    # Set land mask on NSIDC sea ice concentration
    nsidc_data = ma.empty(shape(nsidc_data_tmp))
    nsidc_data[:,:] = 0.0
    nsidc_data[~nsidc_mask.mask] = nsidc_data_tmp[~nsidc_mask.mask]
    nsidc_data[nsidc_mask.mask] = ma.masked

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
    fig = figure(figsize=(20,9))
    # NSIDC
    fig.add_subplot(1,2,1, aspect='equal')
    contourf(nsidc_x, nsidc_y, nsidc_data, lev)
    title('NSIDC', fontsize=24)
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    # CICE
    fig.add_subplot(1,2,2, aspect='equal')
    img = contourf(cice_x, cice_y, cice_data, lev)
    title('CICE', fontsize=24)
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    # Add a horizontal colorbar at the bottom
    cbaxes = fig.add_axes([0.35, 0.04, 0.3, 0.04])
    cbar = colorbar(img, orientation='horizontal', ticks=arange(0,1+0.25,0.25), cax=cbaxes)
    cbar.ax.tick_params(labelsize=20)
    # Add the main title
    suptitle(month_name[month] + ' sea ice concentration', fontsize=30)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    cice_file = raw_input("Path to CICE file: ")
    # Convert month from 1-indexed to 0-indexed
    month = int(raw_input("Month number (1-12): ")) - 1
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    # Make the plot
    nsidc_aice_monthly(cice_file, month, save, fig_name)

    while True:
        # Repeat until the user wants to exit
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to parameters until the user is done
                changes = raw_input("Enter a parameter to change: (1) file path, (2) month number, (3) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New CICE file
                        cice_file = raw_input("Path to CICE file: ")
                    elif int(changes) == 2:
                        # New month
                        month = int(raw_input("Month number (1-12): ")) - 1
                    elif int(changes) == 3:
                        # Switch from save to display, or vice versa
                        save = not save
            if save:
                # Get a new figure name
                fig_name = raw_input("File name for figure: ")
            # Make the plot
            nsidc_aice_monthly(cice_file, month, save, fig_name)
        else:
            break

        

        

        

        

        
        

    
        

    
        

    
    
