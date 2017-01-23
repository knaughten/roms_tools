from numpy import *
from netCDF4 import Dataset, num2date
from matplotlib.pyplot import *
from rotate_vector_cice import *

# Create a 2x2 plot replicating Figure 1 of Holland & Kimura 2016, showing
# sea ice velocity vectors overlaying sea ice concentration for the seasonal
# averages FMA, MMJ, ASO, NDJ.
# Input:
# cice_file = path to CICE output file with 5-day averages, containing at least
#             one complete Feb-Jan period (if there are multiple such periods,
#             this script uses the last one)
# save = optional boolean to save the figure to a file, rather than displaying
#        it on the screen
# fig_name = if save=True, path to the desired filename for figure
def ice_drift_seasonal (cice_file, save=False, fig_name=None):

    # Starting and ending months (1-based) for each season
    start_month = [2, 5, 8, 11]
    end_month = [4, 7, 10, 1]
    # Starting and ending days of the month (1-based) for each season
    start_day = [1, 1, 1, 1]
    end_day = [30, 31, 31, 31]
    # Number of days in each season
    # Assume no leap years, we'll fix this later if needed    
    ndays_season = [89, 92, 92, 92]
    # Season names for titles
    season_names = ['FMA', 'MMJ', 'ASO', 'NDJ']
    # Order of figures (clockwise)
    figure_order = [1, 2, 4, 3]
    # Degrees to radians conversion
    deg2rad = pi/180.0
    # Side length of blocks to average vectors over (can't plot vector at
    # every single point or the plot will be way too crowded)
    block = 15

    # Read CICE grid (including angle) and time values
    id = Dataset(cice_file, 'r')
    lon_tmp = id.variables['TLON'][:-15,:]
    lat_tmp = id.variables['TLAT'][:-15,:]
    angle_tmp = id.variables['ANGLET'][:-15,:]
    # Wrap the periodic boundary by 1 cell
    lon = ma.empty([size(lon_tmp,0), size(lon_tmp,1)+1])
    lat = ma.empty([size(lat_tmp,0), size(lat_tmp,1)+1])
    angle = ma.empty([size(angle_tmp,0), size(angle_tmp,1)+1])
    lon[:,:-1] = lon_tmp
    lon[:,-1] = lon_tmp[:,0]
    lat[:,:-1] = lat_tmp
    lat[:,-1] = lat_tmp[:,0]
    angle[:,:-1] = angle_tmp
    angle[:,-1] = angle_tmp[:,0]
    time_id = id.variables['time']
    # Get the year, month, and day (all 1-based) for each output step
    # These are 5-day averages marked with the next day's date.
    time = num2date(time_id[:], units=time_id.units, calendar=time_id.calendar.lower())

    # Loop backwards through time indices to find the last one we care about
    # (which contains 31 Jan in its averaging period)
    end_t = -1  # Missing value flag
    for t in range(size(time)-1, -1, -1):
        if time[t].month == start_month[0] and time[t].day in range(start_day[0], start_day[0]+5):
            end_t = t
            break
    # Make sure we actually found it
    if end_t == -1:
        print 'Error: ' + cice_file + ' does not contain a complete Feb-Jan period'
        return

    # Continue looping backwards to find the first time index we care about
    # (which contains 1 Feb the previous year in its averaging period)
    start_t = -1  # Missing value flag
    for t in range(end_t-60, -1, -1):
        if time[t].month == start_month[0] and time[t].day in range(start_day[0]+1, start_day[0]+6):
            start_t = t
            break
    # Make sure we actually found it
    if start_t == -1:
        print 'Error: ' + cice_file + ' does not contain a complete Feb-Jan period'
        return

    # Check for leap years
    leap_year = False
    if mod(time[start_t].year, 4) == 0:
        # Years divisible by 4 are leap years
        leap_year = True
        if mod(time[start_t].year, 100) == 0:
            # Unless they're also divisible by 100, in which case they aren't
            # leap years
            leap_year = False
            if mod(time[start_t].year, 400) == 0:
                # Unless they're also divisible by 400, in which case they are
                # leap years after all
                leap_year = True
    if leap_year:
        # Update last day in February
        ndays_season[0] += 1

    # Initialise seasonal averages of CICE output
    aice_tmp = ma.empty([4, size(lon_tmp,0), size(lon_tmp,1)])
    aice_tmp[:,:,:] = 0.0
    uxy_tmp = ma.empty([4, size(lon_tmp,0), size(lon_tmp,1)])
    uxy_tmp[:,:,:] = 0.0
    vxy_tmp = ma.empty([4, size(lon_tmp,0), size(lon_tmp,1)])
    vxy_tmp[:,:,:] = 0.0
    # Process one season at a time
    for season in range(4):
        season_days = 0  # Number of days in season; this will be incremented
        next_season = mod(season+1, 4)

        # Find starting timestep
        start_t_season = -1
        for t in range(start_t, end_t+1):
            if time[t].month == start_month[season] and time[t].day in range(start_day[season]+1, start_day[season]+6):
                start_t_season = t
                break
        # Make sure we actually found it
        if start_t_season == -1:
            print 'Error: could not find starting timestep for season ' + season_names[season]
            return

        # Find ending timestep
        end_t_season = -1
        for t in range(start_t_season+1, end_t+1):
            if time[t].month == start_month[next_season] and time[t].day in range(start_day[next_season], start_day[next_season]+5):
                end_t_season = t
                break
        # Make sure we actually found it
        if end_t_season == -1:
            print 'Error: could not find ending timestep for season ' + season_names[season]
            return

        # Figure out how many of the 5 days averaged in start_t_season are
        # actually within this season
        if time[start_t_season].month == start_month[season] and time[start_t_season].day == start_day[season] + 5:
            # Starting day is in position 1 of 5; we care about all of them
            start_days = 5
        elif time[start_t_season].month == start_month[season] and time[start_t_season].day == start_day[season] + 4:
            # Starting day is in position 2 of 5; we care about the last 4
            start_days = 4
        elif time[start_t_season].month == start_month[season] and time[start_t_season].day == start_day[season]+ 3:
            # Starting day is in position 3 of 5; we care about the last 3
            start_days = 3
        elif time[start_t_season].month == start_month[season] and time[start_t_season].day == start_day[season] + 2:
            # Starting day is in position 4 of 5; we care about the last 2
            start_days = 2
        elif time[start_t_season].month == start_month[season] and time[start_t_season].day == start_day[season] + 1:
            # Starting day is in position 5 of 5; we care about the last 1
            start_days = 1
        else:
            print 'Error for season ' + season_names[season] + ': starting index is month ' + str(time[start_t_season].month) + ', day ' + str(time[start_t_season].day)
            return

        # Start accumulating data weighted by days
        aice_tmp[season,:,:] += id.variables['aice'][start_t_season,:-15,:]*start_days
        uxy_tmp[season,:,:] += id.variables['uvel'][start_t_season,:-15,:]*start_days
        vxy_tmp[season,:,:] += id.variables['vvel'][start_t_season,:-15,:]*start_days
        season_days += start_days

        # Between start_t_season and end_t_season, we want all the days
        for t in range(start_t_season+1, end_t_season):
            aice_tmp[season,:,:] += id.variables['aice'][t,:-15,:]*5
            uxy_tmp[season,:,:] += id.variables['uvel'][t,:-15,:]*5
            vxy_tmp[season,:,:] += id.variables['vvel'][t,:-15,:]*5
            season_days += 5

        # Figure out how many of the 5 days averaged in end_t_season are
        # actually within this season
        if time[end_t_season].month == start_month[next_season] and time[end_t_season].day == start_day[next_season] + 4:
            # Ending day is in position 1 of 5; we care about the first 1
            end_days = 1
        elif time[end_t_season].month == start_month[next_season] and time[end_t_season].day == start_day[next_season] + 3:
            # Ending day is in position 2 of 5; we care about the first 2
            end_days = 2
        elif time[end_t_season].month == start_month[next_season] and time[end_t_season].day == start_day[next_season] + 2:
            # Ending day is in position 3 of 5; we care about the first 3
            end_days = 3
        elif time[end_t_season].month == start_month[next_season] and time[end_t_season].day == start_day[next_season] + 1:
            # Ending day is in position 4 of 5; we care about the first 4
            end_days = 4
        elif time[end_t_season].month == start_month[next_season] and time[end_t_season].day == start_day[next_season]:
            # Ending day is in position 5 of 5; we care about all 5
            end_days = 5
        else:
            print 'Error for season ' + season_names[season] + ': ending index is month ' + str(time[end_t_season].month) + ', day ' + str(time[end_t_season].day)
            return

        aice_tmp[season,:,:] += id.variables['aice'][end_t_season,:-15,:]*end_days
        uxy_tmp[season,:,:] += id.variables['uvel'][end_t_season,:-15,:]*end_days
        vxy_tmp[season,:,:] += id.variables['vvel'][end_t_season,:-15,:]*end_days
        season_days += end_days

        # Check that we got the correct number of days   
        if season_days != ndays_season[season]:
            print 'Error: found ' + str(season_days) + ' days instead of ' + str(ndays_season[season])
            return

        # Finished accumulating data, now convert from sum to average
        aice_tmp[season,:,:] /= season_days
        uxy_tmp[season,:,:] /= season_days
        vxy_tmp[season,:,:] /= season_days

    # Finished reading all CICE data
    id.close()

    # Wrap the periodic boundary
    aice = ma.empty([size(aice_tmp,0), size(aice_tmp,1), size(aice_tmp,2)+1])
    aice[:,:,:-1] = aice_tmp
    aice[:,:,-1] = aice_tmp[:,:,0]
    u_xy = ma.empty([size(uxy_tmp,0), size(uxy_tmp,1), size(uxy_tmp,2)+1])
    u_xy[:,:,:-1] = uxy_tmp
    u_xy[:,:,-1] = uxy_tmp[:,:,0]
    v_xy = ma.empty([size(vxy_tmp,0), size(vxy_tmp,1), size(vxy_tmp,2)+1])
    v_xy[:,:,:-1] = vxy_tmp
    v_xy[:,:,-1] = vxy_tmp[:,:,0]

    # Rotate from local x-y space to lon-lat space
    u, v = rotate_vector_cice(u_xy, v_xy, angle)
    # Calculate speed
    speed = sqrt(u**2 + v**2)
    # Convert velocity to polar coordinates, rotate to account for longitude in
    # circumpolar projection, and convert back to vector components
    theta = arctan2(v, u)
    theta_circ = theta - lon*deg2rad
    u_circ = speed*cos(theta_circ)
    v_circ = speed*sin(theta_circ)

    # Calculate x and y coordinates for plotting circumpolar projection
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)

    # Average x, y, u_circ, and v_circ over block x block intervals
    # Calculate number of blocks
    size0 = int(ceil(size(x,0)/float(block)))
    size1 = int(ceil((size(x,1)-1)/float(block)))
    # Set up arrays for averaged fields
    x_block = ma.empty([size0, size1])
    y_block = ma.empty([size0, size1])
    u_circ_block = ma.empty([4, size0, size1])
    v_circ_block = ma.empty([4, size0, size1])
    # Loop over seasons
    for season in range(4):
        # Set up arrays containing boundary indices
        posn0 = range(0, size(x,0), block)
        posn0.append(size(x,0))
        posn1 = range(0, size(x,1), block)
        posn1.append(size(x,1))
        # Double loop to average each block (can't find a more efficient way to
        # do this)
        for j in range(size0):
            for i in range(size1):
                start0 = posn0[j]
                end0 = posn0[j+1]
                start1 = posn1[i]
                end1 = posn1[i+1]
                if season == 0:
                    # x_block and y_block are season-independent so just do them
                    # for the first season
                    x_block[j,i] = mean(x[start0:end0, start1:end1])
                    y_block[j,i] = mean(y[start0:end0, start1:end1])
                u_circ_block[season,j,i] = mean(u_circ[season, start0:end0, start1:end1])
                v_circ_block[season,j,i] = mean(v_circ[season, start0:end0, start1:end1])

    # Set up colour levels for aice
    lev = linspace(0, 1, num=50)
    # Set boundaries for each side of plot
    bdry1 = -35
    bdry2 = 35
    bdry3 = -33
    bdry4 = 37

    # Make the plot
    fig = figure(figsize=(16,12))
    # Loop over seasons
    for season in range(4):
        ax = fig.add_subplot(2, 2, figure_order[season], aspect='equal')
        # Contour concentration
        img = contourf(x, y, aice[season,:,:], lev, cmap='jet')
        # Add velocity vectors
        q = quiver(x_block, y_block, u_circ_block[season,:,:], v_circ_block[season,:,:], color='black')
        # Configure plot
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
        title(season_names[season], fontsize=24)
        if season == 0:
            # Add colourbar
            cbaxes = fig.add_axes([0.07, 0.6, 0.02, 0.3])
            cbar = colorbar(img, ticks=arange(0,1+0.25,0.25), cax=cbaxes)
            cbar.ax.tick_params(labelsize=16)
        if season == 3:
            # Add 20 cm/s reference vector
            quiverkey(q, 0.07, 0.3, 0.2, '20 cm/s', coordinates='figure', fontproperties={'size':16})
    # Add main title
    suptitle('Sea ice concentration (1) and velocity (m/s)', fontsize=30)
    # Make plots closer together
    subplots_adjust(wspace=0.025,hspace=0.1)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    cice_file = raw_input("Path to CICE file, containing at least one complete Feb-Jan period: ")
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    ice_drift_seasonal(cice_file, save, fig_name)
