nfrom numpy import *
from netCDF4 import Dataset, num2date
from matplotlib.pyplot import *

# Make a figure comparing monthly-averaged sea ice to ocean salt fluxes, 
# from Tamura's dataset to CICE output.
# Input:
# cice_file = path to CICE output file with 5-day averages; if it covers more
#             than one instance of the given month, plot the last one
# month = month number (0-indexed) from 0 to 11
# save = optional boolean to save the figure to a file, rather than displaying
#        it on the screen
# fig_name = if save=True, path to the desired filename for figure
def ssflux_tamura_monthly (cice_file, month, save=False, fig_name=None):

    # Beginning and end of Tamura file paths
    tamura_head = '/short/m68/kaa561/tamura_fluxes/Tamura_ssflux_'
    tamura_tail = '_monthly.nc'
    # Starting and ending days in each month
    # Assume no leap years, we'll fix this later if we need
    start_day = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    end_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    # Name of each month, for the title
    month_name = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
    # Degrees to radians conversion
    deg2rad = pi/180.0
    # Density of freshwater (used by ROMS to convert from kg/m^2/s to psu m/s)
    rho_fw = 1000.0
    # Density of seawater (used by CICE to convert from m/s to kg/m^2/s)
    rho_sw = 1026.0
    # Conversion factor: m/s to cm/day
    mps_to_cmpday = 8.64e6

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

    # Get index of the next month
    next_month = mod(month+1, 12)

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
        # Timestep end_t is likely in the next year, so find the year of start_t
        cice_year = cice_time[start_t].year
    if mod(cice_year, 4) == 0:
        # Years divisible by 4 are leap years
        leap_year = True
        if mod(cice_year, 100) == 0:
            # Unless they're also divisible by 100, in which case they aren't
            # leap years
            leap_year = False
            if mod(cice_.year, 400) == 0:
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

    # Read the fields we need at start_t
    fresh_ai = id.variables['fresh_ai'][start_t,:-15,:]
    sss = id.variables['sss'][start_t,:-15,:]
    rain_ai = id.variables['rain_ai'][start_t,:-15,:]
    fsalt_ai = id.variables['fsalt_ai'][start_t,:-15,:]
    # Start accumulating data weighted by days
    # Convert to units of psu m/s (equivalent to kg/m^2/s of salt)
    # Subtract rain from freshwater flux, since Tamura doesn't count precip
    cice_data_tmp += -1/rho_fw*((fresh_ai-rain_ai)*sss*rho_sw/mps_to_cmpday - fsalt_ai*1e3)*start_days
    num_days += start_days

    # Between start_t and end_t, we want all the days
    for t in range(start_t+1, end_t):
        fresh_ai = id.variables['fresh_ai'][t,:-15,:]
        sss = id.variables['sss'][t,:-15,:]
        rain_ai = id.variables['rain_ai'][t,:-15,:]
        fsalt_ai = id.variables['fsalt_ai'][t,:-15,:]
        cice_data_tmp += -1/rho_fw*((fresh_ai-rain_ai)*sss*rho_sw/mps_to_cmpday - fsalt_ai*1e3)*5
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

    fresh_ai = id.variables['fresh_ai'][end_t,:-15,:]
    sss = id.variables['sss'][end_t,:-15,:]
    rain_ai = id.variables['rain_ai'][end_t,:-15,:]
    fsalt_ai = id.variables['fsalt_ai'][end_t,:-15,:]
    cice_data_tmp += -1/rho_fw*((fresh_ai-rain_ai)*sss*rho_sw/mps_to_cmpday - fsalt_ai*1e3)*end_days
    num_days += end_days

    # Check that we got the correct number of days
    if num_days != end_day[month]:
        print 'Error: found ' + str(num_days) + ' days instead of ' + str(end_day[month])
        return

    # Finished accumulating data
    id.close()
    # Convert from sum to average
    cice_data_tmp /= num_days
    # Multiply by 1e6 so colour bar is easier to read
    cice_data_tmp *= 1e6

    # Wrap periodic boundary
    cice_data = ma.empty([size(cice_data_tmp,0), size(cice_data_tmp,1)+1])
    cice_data[:,:-1] = cice_data_tmp
    cice_data[:,-1] = cice_data_tmp[:,0]

    # Read Tamura data and grid
    id = Dataset(tamura_head + str(cice_year) + tamura_tail, 'r')
    tamura_lon = id.variables['longitude'][:,:]
    tamura_lat = id.variables['latitude'][:,:]
    tamura_data = id.variables['ssflux'][month,:,:]
    id.close()
    # Appply land mask and multiply by 1e6
    tamura_data = ma.masked_where(isnan(tamura_data), tamura_data)
    tamura_data *= 1e6

    # Convert both grids to spherical coordinates
    cice_x = -(cice_lat+90)*cos(cice_lon*deg2rad+pi/2)
    cice_y = (cice_lat+90)*sin(cice_lon*deg2rad+pi/2)
    tamura_x = -(tamura_lat+90)*cos(tamura_lon*deg2rad+pi/2)
    tamura_y = (tamura_lat+90)*sin(tamura_lon*deg2rad+pi/2)

    # Choose colour levels
    bound = 20.0
    lev = linspace(-bound, bound, num=50)
    # Find boundaries for each side of plot based on extent of grids
    bdry1 = max(amin(cice_x), amin(tamura_x))
    bdry2 = min(amax(cice_x), amax(tamura_x))
    bdry3 = max(amin(cice_y), amin(tamura_y))
    bdry4 = min(amax(cice_y), amax(tamura_y))

    # Plot
    fig = figure(figsize=(20,9))
    # Tamura
    fig.add_subplot(1,2,1, aspect='equal')
    contourf(tamura_x, tamura_y, tamura_data, lev, cmap='RdYlBu_r')
    title('Tamura', fontsize=24)
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    # CICE
    fig.add_subplot(1,2,2, aspect='equal')
    img = contourf(cice_x, cice_y, cice_data, lev, cmap='RdYlBu_r')
    title('CICE', fontsize=24)
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    # Add a horizontal colourbar at the bottom
    cbaxes = fig.add_axes([0.3, 0.04, 0.4, 0.04])
    cbar = colorbar(img, orientation='horizontal', ticks=arange(-1,1+0.5, 0.5), cax=cbaxes, extend='both')
    cbar.ax.tick_params(labelsize=18)
    # Add the main title
    suptitle(r'Ice-to-ocean salt flux (10$^{-6}$ kg/m$^2$/s, ' + month_name[month] + ' ' + str(cice_year), fontsize=30)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    cice_file = raw_input("Path to CICE file: ")
    # Convert month from 1-indexed to 0-indexed
    month = int(raw_input("Month number (1-12): "))-1
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    # Make the plot
    ssflux_tamura_monthly(cice_file, month, save, fig_name)

    while True:
        # Repeat until the user wants to exit
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the parameters until the user is done
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
                        month = int(raw_input("Month number (1-12): "))-1
                    elif int(changes) == 3:
                        # Switch from save to display, or vice versa
                        save = not save
            if save:
                # Get a new figure name
                fig_name = raw_input("File name for figure: ")
            # Make the plot
            ssflux_tamura_monthly(cice_file, month, save, fig_name)
        else:
            break
            

    

    
        
    
        
