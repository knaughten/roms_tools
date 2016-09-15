from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from cartesian_grid_2d import *

# Calculate and plot timeseries of basal mass loss and area-averaged ice shelf
# melt rates from major ice shelves and from the entire continent during a 
# ROMS simulation.
# Input:
# file_path = path to ocean history/averages file
# log_path = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
def timeseries_massloss (file_path, log_path):

    # Titles and figure names for each ice shelf
    names = ['All Ice Shelves', 'Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    fig_names = ['total_massloss.png', 'larsen_d.png', 'larsen_c.png', 'wilkins_georgevi_stange.png', 'ronne_filchner.png', 'abbot.png', 'pig.png', 'thwaites.png', 'dotson.png', 'getz.png', 'nickerson.png', 'sulzberger.png', 'mertz.png', 'totten_moscowuni.png', 'shackleton.png', 'west.png', 'amery.png', 'princeharald.png', 'baudouin_borchgrevink.png', 'lazarev.png', 'nivl.png', 'fimbul_jelbart_ekstrom.png', 'brunt_riiserlarsen.png', 'ross.png']
    # Limits on longitude and latitude for each ice shelf
    # These depend on the source geometry, in this case RTopo 1.05
    # Note there is one extra index at the end of each array; this is because
    # the Ross region crosses the line 180W and therefore is split into two
    lon_min = [-180, -62.67, -65.5, -79.17, -85, -104.17, -102.5, -108.33, -114.5, -135.67, -149.17, -155, 144, 115, 94.17, 80.83, 65, 33.83, 19, 12.9, 9.33, -10.05, -28.33, -180, 158.33]
    lon_max = [180, -59.33, -60, -66.67, -28.33, -88.83, -99.17, -103.33, -111.5, -114.33, -140, -145, 146.62, 123.33, 102.5, 89.17, 75, 37.67, 33.33, 16.17, 12.88, 7.6, -10.33, -146.67, 180]
    lat_min = [-90, -73.03, -69.35, -74.17, -83.5, -73.28, -75.5, -75.5, -75.33, -74.9, -76.42, -78, -67.83, -67.17, -66.67, -67.83, -73.67, -69.83, -71.67, -70.5, -70.75, -71.83, -76.33, -85, -84.5]
    lat_max = [-30, -69.37, -66.13, -69.5, -74.67, -71.67, -74.17, -74.67, -73.67, -73, -75.17, -76.41, -66.67, -66.5, -64.83, -66.17, -68.33, -68.67, -68.33, -69.33, -69.83, -69.33, -71.5, -77.77, -77]
    # Observed mass loss (Rignot 2013) and uncertainty for each ice shelf, in Gt/y
    obs_massloss = [1325, 1.4, 20.7, 135.4, 155.4, 51.8, 101.2, 97.5, 45.2, 144.9, 4.2, 18.2, 7.9, 90.6, 72.6, 27.2, 35.5, -2, 21.6, 6.3, 3.9, 26.8, 9.7, 47.7]
    obs_massloss_error = [235, 14, 67, 40, 45, 19, 8, 7, 4, 14, 2, 3, 3, 8, 15, 10, 23, 3, 18, 2, 2, 14, 16, 34]
    # Observed ice shelf melt rates and uncertainty
    obs_ismr = [0.85, 0.1, 0.4, 3.1, 0.3, 1.7, 16.2, 17.7, 7.8, 4.3, 0.6, 1.5, 1.4, 7.7, 2.8, 1.7, 0.6, -0.4, 0.4, 0.7, 0.5, 0.5, 0.1, 0.1]
    obs_ismr_error = [0.1, 0.6, 1, 0.8, 0.1, 0.6, 1, 1, 0.6, 0.4, 0.3, 0.3, 0.6, 0.7, 0.6, 0.7, 0.4, 0.6, 0.4, 0.2, 0.2, 0.2, 0.2, 0.1]          
    # Density of ice in kg/m^3
    rho_ice = 916

    old_time = []
    # Check if the log file exists
    if exists(log_path):
        print 'Reading previously calculated values'
        f = open(log_path, 'r')
        # Skip the first line (header for time array)
        f.readline()
        for line in f:
            try:
                old_time.append(float(line))
            except(ValueError):
                # Reached the header for the next variable
                break
        # Set up array for mass loss values at each ice shelf
        old_massloss = empty([len(names), len(old_time)])
        index = 0
        # Loop over ice shelves
        while index < len(names):
            t = 0
            for line in f:
                try:
                    old_massloss[index, t] = float(line)
                    t += 1
                except(ValueError):
                    # Reached the header for the next ice shelf
                    break
            index += 1

    # Calculate dA (masked with ice shelf mask) and lon and lat coordinates
    print 'Analysing grid'
    dA, lon, lat = calc_grid(file_path)

    # Read time data and convert from seconds to years
    id = Dataset(file_path, 'r')
    new_time = id.variables['ocean_time'][:]/(365.25*24*60*60)
    if exists(log_path):
        # Concatenate with time values from log file
        start_t = len(old_time)
        time = old_time
        for t in range(size(new_time)):
            time.append(new_time[t])
        time = array(time)
    else:
        start_t = 0
        time = new_time
    # Set up array of mass loss values
    massloss = empty([len(names), size(time)])
    if exists(log_path):
        # Fill first start_t timesteps with existing values
        massloss[:,0:start_t] = old_massloss[:,:]

    print 'Reading data'
    # Read melt rate and convert from m/s to m/y
    ismr = id.variables['m'][:,:-15,:-3]*365.25*24*60*60
    id.close()

    print 'Setting up arrays'
    # Set up array of masked area values for each ice shelf
    dA_masked = ma.empty([len(names), size(dA,0), size(dA,1)])
    # Set up array of conversion factors from mass loss to area-averaged melt
    # rate for each ice shelf
    factors = empty(len(names))
    for index in range(len(names)):
        # Mask dA for the current ice shelf (note dA is already masked
        # with the global ice shelf mask, so just restrict the indices
        # to isolate the given ice shelf)
        if index == len(names)-1:
            # Ross region is split into two
            dA_tmp = ma.masked_where(((lon < lon_min[index]) + (lon > lon_max[index]) + (lat < lat_min[index]) + (lat > lat_max[index]))*((lon < lon_min[index+1]) + (lon > lon_max[index+1]) + (lat < lat_min[index+1]) + (lat > lat_max[index+1])), dA)
        else:
            dA_tmp = ma.masked_where((lon < lon_min[index]) + (lon > lon_max[index]) + (lat < lat_min[index]) + (lat > lat_max[index]), dA)
        dA_masked[index,:,:] = dA_tmp[:,:]
        # Calculate area of this ice shelf
        area_tmp = sum(dA_masked[index,:,:])
        # Calculate conversion factor
        factors[index] = 1e12/(rho_ice*area_tmp)
        print 'Area of ' + names[index] + ': ' + str(area_tmp) + ' m^2'

    # Build timeseries
    print 'Calculating variables'
    for t in range(start_t, size(time)):
        # Loop over ice shelves
        for index in range(len(names)):
            # Integrate ice shelf melt rate over area to get volume loss
            volumeloss = sum(ismr[t-start_t,:,:]*dA_masked[index,:,:])
            # Convert to mass loss in Gt/y
            massloss[index, t] = 1e-12*rho_ice*volumeloss

    # Plot each timeseries
    print 'Plotting'
    for index in range(len(names)):
        # Calculate the bounds on observed mass loss and melt rate
        massloss_low = obs_massloss[index] - obs_massloss_error[index]
        massloss_high = obs_massloss[index] + obs_massloss_error[index]
        ismr_low = obs_ismr[index] - obs_ismr_error[index]
        ismr_high = obs_ismr[index] + obs_ismr_error[index]
        # Set up plot: mass loss and melt rate are directly proportional (with
        # a different constant of proportionality for each ice shelf depending
        # on its area) so plot one line with two y-axes
        fig, ax1 = subplots()
        ax1.plot(time, massloss[index,:], color='black')
        # In blue, add dashed lines for observed mass loss
        ax1.axhline(massloss_low, color='b', linestyle='dashed')
        ax1.axhline(massloss_high, color='b', linestyle='dashed')
        # Make sure y-limits won't cut off observed melt rate
        ymin = amin([ismr_low/factors[index], massloss_low, amin(massloss[index,:])])
        ymax = amax([ismr_high/factors[index], massloss_high, amax(massloss[index,:])])
        # Adjust y-limits to line up with ticks
        ticks = ax1.get_yticks()
        min_tick = ticks[0]
        max_tick = ticks[-1]
        dtick = ticks[1]-ticks[0]
        while min_tick >= ymin:
            min_tick -= dtick
        while max_tick <= ymax:
            max_tick += dtick
        ax1.set_ylim([min_tick, max_tick])
        # Title and ticks in blue for this side of the plot
        ax1.set_ylabel('Basal Mass Loss (Gt/y)', color='b')
        for t1 in ax1.get_yticklabels():
            t1.set_color('b')
        ax1.set_xlabel('Years')
#        setp(ax1.get_xticklabels(), fontsize=15)
#        setp(ax1.get_yticklabels(), fontsize=15)
        ax1.grid(True)
        # Twin axis for melt rates
        ax2 = ax1.twinx()
        # Make sure the scales line up
        limits = ax1.get_ylim()
        ax2.set_ylim([limits[0]*factors[index], limits[1]*factors[index]])
        # In red, add dashed lines for observed ice shelf melt rates
        ax2.axhline(ismr_low, color='r', linestyle='dashed')
        ax2.axhline(ismr_high, color='r', linestyle='dashed')
        # Title and ticks in red for this side of the plot
        ax2.set_ylabel('Area-Averaged Ice Shelf Melt Rate (m/y)', color='r')
        for t2 in ax2.get_yticklabels():
            t2.set_color('r')
#        setp(ax2.get_yticklabels(), fontsize=15)
        # Name of the ice shelf for the main title
        title(names[index])
        fig.savefig(fig_names[index])
        
    print 'Saving results to log file'
    f = open(log_path, 'w')
    f.write('Time (years):\n')
    for t in range(size(time)):
        f.write(str(time[t]) + '\n')
    for index in range(len(names)):
        f.write(names[index] + ' Basal Mass Loss\n')
        for t in range(size(time)):
            f.write(str(massloss[index, t]) + '\n')
    f.close()


# Given the path to a ROMS grid file, calculate differential of area and
# longitude and latitude.
# Input: file_path = string containing path to ROMS history/averages file
# Output:
# dA = differential of area on the 2D rho-grid, masked with zice
# lon = longitude values on the rho-grid, from -180 to 180
# lat = latitude values on the rho-grid
def calc_grid (file_path):

    # Read grid variables
    id = Dataset(file_path, 'r')
    lon = id.variables['lon_rho'][:-15,:-3]
    lat = id.variables['lat_rho'][:-15,:-3]
    zice = id.variables['zice'][:-15,:-3]
    id.close()

    # Calculate dx and dy in another script
    dx, dy = cartesian_grid_2d(lon, lat)

    # Calculate dA and mask with zice
    dA = ma.masked_where(zice==0, dx*dy)

    # Save dimensions
    num_lat = size(lon, 0)
    num_lon = size(lon, 1)

    # Make longitude values go from -180 to 180, not 0 to 360
    index = lon > 180
    lon[index] = lon[index] - 360

    return dA, lon, lat


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input('Enter path to ocean history/averages file: ')
    log_path = raw_input('Enter path to log file to save values and/or read previously calculated values: ')

    timeseries_massloss(file_path, log_path)

