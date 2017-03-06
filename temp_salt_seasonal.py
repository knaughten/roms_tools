from netCDF4 import Dataset, num2date
from numpy import *
from matplotlib.pyplot import *
from calc_z import *

# Make a 4x2 plot showing lat vs. depth slices of seasonally averaged 
# temperature (top) and salinity (bottom) at the given longitude, over the
# last year of simulation.
# Input:
# file_path = path to ROMS output file, containing at least one complete 
#             Dec-Nov period (if there are multiple such periods, the last one
#             will be used for seasonal averages)
# lon0 = the specific longitude to plot (between -180 and 180)
# depth_bdry = deepest depth to plot (negative, in m)
# save = optional boolean flag; if True, the figure will be saved with file name
#        fig_name; if False, the figure will display on the screen
# fig_name = optional string containing filename for figure, if save=True
def temp_salt_seasonal (file_path, lon0, depth_bdry, save=False, fig_name=None):

    # Path to SOSE seasonal climatology file
    sose_file = '../SOSE_seasonal_climatology.nc'

    # Grid parameters
    theta_s = 4.0
    theta_b = 0.9
    hc = 40
    N = 31

    # Starting and ending months (1-based) for each season
    start_month = [12, 3, 6, 9]
    end_month = [2, 5, 8, 11]
    # Starting and ending days of the month (1-based) for each season
    # Assume no leap years, we'll fix this later if we need
    start_day = [1, 1, 1, 1]
    end_day = [28, 31, 31, 30]
    # Number of days in each season (again, ignore leap years for now)
    ndays_season = [90, 92, 92, 91]
    # Season names for titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']

    # Bounds on colour scale
    temp_min = -2.5
    temp_max = 7.5
    temp_ticks = 2
    salt_min = 33.8
    salt_max = 34.8
    salt_ticks = 0.2

    # Choose what to write on the title about longitude
    if lon0 < 0:
        lon_string = 'T/S slices at ' + str(int(round(-lon0))) + r'$^{\circ}$W'
    else:
        lon_string = 'T/S slices at ' + str(int(round(lon0))) + r'$^{\circ}$E'
    # Edit longitude bounds to be from 0 to 360, to fit with ROMS convention
    if lon0 < 0:
        lon0 += 360

    print 'Processing ROMS data'

    # Read grid and time values
    id = Dataset(file_path, 'r')
    h = id.variables['h'][:-15,:]
    zice = id.variables['zice'][:-15,:]
    lon_roms_2d = id.variables['lon_rho'][:-15,:]
    lat_roms_2d = id.variables['lat_rho'][:-15,:]
    time_id = id.variables['ocean_time']
    # Get the year, month, and day (all 1-based) for each output step
    # These are 5-day averages marked with the middle day's date
    time = num2date(time_id[:], units=time_id.units, calendar=time_id.calendar.lower())

    # Loop backwards through time indices to find the last one we care about
    # (which contains 30 November in its averaging period)
    end_t = -1  # Missing value flag
    for t in range(size(time)-1, -1, -1):
        if time[t].month == end_month[-1] and time[t].day in range(end_day[-1]-2, end_day[-1]+1):
            end_t = t
            break
        if time[t].month == start_month[0] and time[t].day in range(start_day[0], start_day[0]+2):
            end_t = t
            break
    # Make sure we actually found it
    if end_t == -1:
        print 'Error: ' + file_path + ' does not contain a complete Dec-Nov period'
        return

    # Continue looping backwards to find the first time index we care about
    # (which contains 1 December the previous year in its averaging period)
    start_t = -1  # Missing value flag
    for t in range(end_t-60, -1, -1):
        if time[t].month == end_month[-1] and time[t].day in range(end_day[-1]-1, end_day[-1]+1):
            start_t = t
            break
        if time[t].month == start_month[0] and time[t].day in range(start_day[0], start_day[0]+3):
            start_t = t
            break
    # Make sure we actually found it
    if start_t == -1:
        print 'Error: ' + file_path + ' does not contain a complete Dec-Nov period'
       return

    # Check if end_t occurs on a leap year
    leap_year = False
    if mod(time[end_t].year, 4) == 0:
        # Years divisible by 4 are leap years
        leap_year = True
        if mod(time[end_t].year, 100) == 0:
            # Unless they're also divisible by 100, in which case they aren't
            # leap years
            leap_year = False
            if mod(time[end_t].year, 400) == 0:
                # Unless they're also divisible by 400, in which case they are
                # leap years after all
                leap_year = True
    if leap_year:
        # Update last day in February
        end_day[0] += 1
        ndays_season[0] += 1

    # Initialise seasonal averages
    temp_3d_roms = ma.empty([4, N, size(lon_roms_2d,0), size(lon_roms_2d,1)])
    temp_3d_roms[:,:,:,:] = 0.0
    salt_3d_roms = ma.empty([4, N, size(lon_roms_2d,0), size(lon_roms_2d,1)])
    salt_3d_roms[:,:,:,:] = 0.0
    # Process one season at a time
    for season in range(4):
        print 'Calculating seasonal averages for ' + season_names[season]
        season_days = 0  # Number of days in season; this will be incremented
        next_season = mod(season+1, 4)

        # Find starting timestep
        start_t_season = -1
        for t in range(start_t, end_t+1):
            if time[t].month == end_month[season-1] and time[t].day in range(end_day[season-1]-1, end_day[season-1]+1):
                start_t_season = t
                break
            if time[t].month == start_month[season] and time[t].day in range(start_day[season], start_day[season]+3):
                start_t_season = t
                break
        # Make sure we actually found it
        if start_t_season == -1:
            print 'Error: could not find starting timestep for season ' + season_title[season]
            return

        # Find ending timestep
        end_t_season = -1
        for t in range(start_t_season+1, end_t+1):
            if time[t].month == end_month[season] and time[t].day in range(end_day[season]-2, end_day[season]+1):
                end_t_season = t
                break
            if time[t].month == start_month[next_season] and time[t].day in range(start_day[next_season], start_day[next_season]+2):
                end_t_season = t
                break
        # Make sure we actually found it
        if end_t_season == -1:
            print 'Error: could not find ending timestep for season ' + season_title[season]
            return

        # Figure out how many of the 5 days averaged in start_t_season are
        # actually within this season
        if time[start_t_season].month == start_month[season] and time[start_t_season].day == start_day[season]+2:
            # Starting day is in position 1 of 5; we care about all of them
            start_days = 5
        elif time[start_t_season].month == start_month[season] and time[start_t_season].day == start_day[season]+1:
            # Starting day is in position 2 of 5; we care about the last 4
            start_days = 4
        elif time[start_t_season].month == start_month[season] and time[start_t_season].day == start_day[season]:
            # Starting day is in position 3 of 5; we care about the last 3
            start_days = 3
        elif time[start_t_season].month == end_month[season-1] and time[start_t_season].day == end_day[season-1]:
            # Starting day is in position 4 of 5; we care about the last 2
            start_days = 2
        elif time[start_t_season].month == end_month[season-1] and time[start_t_season].day == end_day[season-1]-1:
            # Starting day is in position 5 of 5; we care about the last 1
            start_days = 1
        else:
            print 'Error for season ' + season_title[season] + ': starting index is month ' + str(time[start_t_season].month) + ', day ' + str(time[start_t_season].day)
            return

        # Start accumulating data weighted by days
        temp_3d_roms[season,:,:,:] += id.variables['temp'][start_t_season,:,:-15,:]*start_days
        salt_3d_roms[season,:,:,:] += id.variables['salt'][start_t_season,:,:-15,:]*start_days
        season_days += start_days

        # Between start_t_season and end_t_season, we want all the days
        for t in range(start_t_season+1, end_t_season):
            temp_3d_roms[season,:,:,:] += id.variables['temp'][t,:,:-15,:]*5
            salt_3d_roms[season,:,:,:] += id.variables['salt'][t,:,:-15,:]*5
            season_days += 5

        # Figure out how many of the 5 days averaged in end_t_season are
        # actually within this season
        if time[end_t_season].month == start_month[next_season] and time[end_t_season].day == start_day[next_season]+1:
            # Ending day is in position 1 of 5; we care about the first 1
            end_days = 1
        elif time[end_t_season].month == start_month[next_season] and time[end_t_season].day == start_day[next_season]:
            # Ending day is in position 2 of 5; we care about the first 2
            end_days = 2
        elif time[end_t_season].month == end_month[season] and time[end_t_season].day == end_day[season]:
            # Ending day is in position 3 of 5; we care about the first 3
            end_days = 3
        elif time[end_t_season].month == end_month[season] and time[end_t_season].day == end_day[season]-1:
            # Ending day is in position 4 of 5; we care about the first 4
            end_days = 4
        elif time[end_t_season].month == end_month[season] and time[end_t_season].day == end_day[season]-2:
            # Ending day is in position 5 of 5; we care about all 5
            end_days = 5
        else:
            print 'Error for season ' + season_title[season] + ': ending index is month ' + str(time[end_t_season].month) + ', day ' + str(time[end_t_season].day)
            return

        temp_3d_roms[season,:,:,:] += id.variables['temp'][end_t_season,:,:-15,:]*end_days
        salt_3d_roms[season,:,:,:] += id.variables['salt'][end_t_season,:,:-15,:]*end_days
        season_days += end_days

        # Check that we got the correct number of days
        if season_days != ndays_season[season]:
            print 'Error: found ' + str(season_days) + ' days instead of ' + str(ndays_season[season])
            return

        # Finished accumulating data, now convert from sum to average
        temp_3d_roms[season,:,:,:] /= season_days
        salt_3d_roms[season,:,:,:] /= season_days

    # Finished reading data
    id.close()

    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    z_roms_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N)

    # Calculate zonal slices for each season
    temp_roms = ma.empty([4, N, size(lat_roms_2d,0)])
    temp_roms[:,:,:] = 0.0
    salt_roms = ma.empty([4, N, size(lat_roms_2d,0)])
    salt_roms[:,:,:] = 0.0
    for season in range(4):
        print 'Calculating zonal slices for ' + season_names[season]
        temp_tmp, z_roms, lat_roms = interp_lon_roms(temp_3d_roms[season,:,:,:], z_roms_3d, lat_roms_2d, lon_roms_2d, lon0)
        temp_roms[season,:,:] = temp_tmp
        salt_tmp, z_roms, lat_roms = interp_lon_roms(salt_3d_roms[season,:,:,:], z_roms_3d, lat_roms_2d, lon_roms_2d, lon0)
        salt_roms[season,:,:] = salt_tmp

    # Set colour levels
    lev1 = linspace(temp_min, temp_max, num=50)
    lev2 = linspace(salt_min, salt_max, num=50)

    # Choose boundaries based on extent of ROMS grid
    sbdry = amin(lat_roms)
    nbdry = amax(lat_roms)

    # Plot
    print 'Plotting'
    fig = figure(figsize=(20,9))
    # Loop over seasons
    for season in range(4):
        # Temperature
        fig.add_subplot(2, 4, season+1)
        img = contourf(lat_roms, z_roms, temp_roms[season,:,:], lev1, cmap='jet', extend='both')
        xlim([sbdry, nbdry])
        ylim([depth_bdry, 0])
        title('Temperature (' + season_names[season] + ')', fontsize=24)
        if season == 0:
            ylabel('depth (m)', fontsize=18)
        if season == 3:
            cbaxes1 = fig.add_axes([0.92, 0.55, 0.01, 0.3])
            cbar1 = colorbar(img, cax=cbaxes1, ticks=arange(temp_min, temp_max+temp_ticks, temp_ticks))
            cbar1.ax.tick_params(labelsize=16)
        # Salinity
        fig.add_subplot(2, 4, season+5)
        img = contourf(lat_roms, z_roms, salt_roms[season,:,:], lev2, cmap='jet', extend='both')
        xlim([sbdry, nbdry])
        ylim([depth_bdry, 0])
        title('Salinity (' + season_names[season] + ')', fontsize=24)
        if season == 0:
            ylabel('depth (m)', fontsize=18)
        xlabel('Latitude', fontsize=18)
        if season == 3:
            cbaxes2 = fig.add_axes([0.92, 0.15, 0.01, 0.3])
            cbar2 = colorbar(img, cax=cbaxes2, ticks=arange(salt_min, salt_max+salt_ticks, salt_ticks))
            cbar2.ax.tick_params(labelsize=16)
    # Add the main title
    suptitle(lon_string, fontsize=30)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()    


# Linearly interpolate ROMS data, z, and latitude to the specified longitude.
# Input:
# data_3d = array of data, dimension depth x lat x lon
# z_3d = array of depth values (negative, in metres), dimension depth x lat x lon
# lat_2d = array of latitude values, dimension lat x lon
# lon_2d = array of longitude values, dimension lat x lon (between -180 and 180)
# lon0 = longitude to interpolate to (between -180 and 180)
# Output:
# data = array of data interpolated to lon0, dimension depth x lat
# z = array of depth values interpolated to lon0, dimension depth x lat
# lat = array of latitude values interpolated to lon0, dimension depth x lat
def interp_lon_roms (data_3d, z_3d, lat_2d, lon_2d, lon0):

    # Save dimensions
    num_depth = size(data_3d, 0)
    num_lat = size(data_3d, 1)
    num_lon = size(data_3d, 2)
    # Set up output arrays
    data = ma.empty([num_depth, num_lat])
    z = ma.empty([num_depth, num_lat])
    lat = ma.empty([num_depth, num_lat])

    # Loop over latitudes; can't find a cleaner way to do this
    for j in range(num_lat):
        # Extract the longitude values of this slice
        lon_tmp = lon_2d[j,:]
        # Get indices and coefficients for interpolation
        ie, iw, coeffe, coeffw = interp_lon_roms_helper(lon_tmp, lon0)        
        data[:,j] = coeffe*data_3d[:,j,ie] + coeffw*data_3d[:,j,iw]
        z[:,j] = coeffe*z_3d[:,j,ie] + coeffw*z_3d[:,j,iw]
        lat[:,j] = coeffe*lat_2d[j,ie] + coeffw*lat_2d[j,iw]

    return data, z, lat


# Calculate indices and coefficients for linear interpolation of ROMS longitude.
# This takes care of all the mod 360 nonsense.
# Input:
# lon = 1D array of longitude values (straight out of ROMS i.e. between slightly < 0 and slightly > 360)
# lon0 = longitude to interpolate to (between 0 and 360)
# Output:
# ie, iw, coeffe, coeffw = integers (ie and iw) and coefficients (coeffe and 
#                          coeffw) such that coeffe*lon[ie] + coeffw*lon[iw] = 
#                          lon0, which will also hold for any variable on this 
#                          longitude grid. ie is the index of the nearest point
#                          to the east of lon0; iw the nearest to the west.
def interp_lon_roms_helper (lon, lon0):

    if lon0 < amin(lon) or lon0 > amax(lon):
        # Special case: lon0 on periodic boundary
        # Be careful with mod 360 here

        # Find the periodic boundary
        dlon = lon[1:] - lon[0:-1]
        bdry = argmax(abs(dlon))
        if dlon[bdry] < -300:
            # Jumps from almost 360 to just over 0
            iw = bdry
            ie = bdry + 1
        else:
            # Periodic boundary lines up with the array boundary
            iw = size(lon) - 1
            ie = 0
        # Calculate difference between lon0 and lon[iw], mod 360 if necessary
        dlon_num = lon0 - lon[iw]
        if dlon_num < -300:
            dlon_num += 360
        # Calculate difference between lon[ie] and lon[iw], mod 360
        dlon_den = lon[ie] - lon[iw] + 360

    else:
        # General case

        # Add or subtract 360 from longitude values which wrap around
        # so that longitude increases monotonically from west to east
        i = arange(1, size(lon)+1)
        index1 = nonzero((i > 1200)*(lon < 100))
        lon[index1] = lon[index1] + 360
        index2 = nonzero((i < 200)*(lon > 300))
        lon[index2] = lon[index2] - 360

        # Take mod 360 of lon0 if necessary
        if all(lon < lon0):
            lon0 -= 360
        if all(lon > lon0):
            lon0 += 360
        
        # Find the first index eastward of lon0
        ie = nonzero(lon > lon0)[0][0]
        # The index before it will be the last index westward of lon0
        iw = ie - 1

        dlon_num = lon0 - lon[iw]
        dlon_den = lon[ie] - lon[iw]

    if dlon_num > 5 or dlon_den > 5:
        print 'interp_lon_helper: Problem at periodic boundary'
        return
    coeff1 = dlon_num/dlon_den
    coeff2 = 1 - coeff1

    return ie, iw, coeff1, coeff2


# Linearly interpolate SOSE data to the specified longitude.
# Input:
# data_3d = arary of data, dimension depth x lat x lon
# lon = 1D array of longitude values (between 0 and 360)
# lon0 = longitude to interpolate to (between 0 and 360)
# Output:
# data = array of data interpolated to lon0, dimension depth x lat
def interp_lon_sose (data_3d, lon, lon0):

    if lon0 < lon[0] or lon0 > lon[-1]:
        # Special case: lon0 on periodic boundary
        # Be careful with mod 360 here
        iw = size(lon)-1
        ie = 0
        # Calculate difference between lon0 and lon[iw], mod 360 if necessary
        dlon_num = lon0 - lon[iw]
        if dlon_num < -300:
            dlon_num += 360
        # Calculate difference between lon[ie] and lon[iw], mod 360
        dlon_den = lon[ie] - lon[iw] + 360
    else:
        # General case
        # Find the first index eastwards of lon0
        ie = nonzero(lon > lon0)[0][0]
        # The index before it will be the last index westward of lon0
        iw = ie - 1
        dlon_num = lon0 - lon[iw]
        dlon_den = lon[ie] - lon[iw]
    # Coefficients for interpolation
    coeff1 = dlon_num/dlon_den
    coeff2 = 1 - coeff1

    # Interpolate
    data = coeff1*data_3d[:,:,ie] + coeff2*data_3d[:,:,iw]
    return data


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to ocean averages file, containing at least one complete Dec-Nov period: ")
    lon0 = float(raw_input("Enter longitude (-180 to 180): "))
    depth_bdry = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    temp_salt_seasonal(file_path, lon0, depth_bdry, save, fig_name)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the input parameters; repeat until the user is finished
                changes = raw_input("Enter a parameter to change: (1) file path, (2) longitude, (3) deepest depth, (4) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New file path
                        file_path = raw_input("Path to ocean averages file, containing at least one complete Dec-Nov period: ")
                    elif int(changes) == 2:
                        # New longitude
                        lon0 = float(raw_input("Enter longitude (-180 to 180): "))
                    elif int(changes) == 3:
                        # New depth bound
                        depth_bdry = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
                    elif int(changes) == 4:
                        # Change from save to display, or vice versa
                        save = not save
            if save:
                # Get file name for figure
                fig_name = raw_input("File name for figure: ")

            # Make the plot
            temp_salt_seasonal(file_path, lon0, depth_bdry, save, fig_name)
        else:
            break
        
                

    
