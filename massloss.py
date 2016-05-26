from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from cartesian_grid_2d import *

# Calculate and plot timeseries of basal mass loss and area-averaged ice shelf
# melt rates from major ice shelves during a ROMS simulation.
# Input:
# file_path = path to ocean history/averages file
# log_path = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
def massloss (file_path, log_path):

    # Titles and figure names for each ice shelf
    names = ['Amery Ice Shelf', 'Ross Ice Shelf', 'Getz Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Abbot Ice Shelf', 'George VI Ice Shelf', 'Larsen C Ice Shelf', 'Ronne-Filchner Ice Shelf', 'Brunt and Riiser-Larsen Ice Shelves', 'Fimbul and Jelbart Ice Shelves']
    fig_names = ['amery_massloss.png', 'ross_massloss.png', 'getz_massloss.png', 'pig_massloss.png', 'abbot_massloss.png', 'george_vi_massloss.png', 'larsen_c_massloss.png', 'ronne_filchner_massloss.png', 'brunt_riiser_larsen_massloss.png', 'fimbul_jelbart_massloss.png']
    # Limits on i- and j- coordinates for each ice shelf; this will vary
    # depending on the grid
    # Note there is one extra index at the end of each array; this is because
    # the Fimbul-Jelbart region crosses the periodic boundary and therefore
    # is split into two
    i_min = [250, 700, 905, 1015, 1000, 1100, 1165, 1060, 1280, 1375, 1]
    i_max = [350, 872, 975, 1030, 1090, 1155, 1190, 1240, 1369, 1443, 12]
    j_min = [1,   20,  150, 140,  160,  150,  187,  1,    65,   80,   100]
    j_max = [125, 123, 175, 160,  185,  200,  220,  135,  116,  150,  120]
    # Observed mass loss (Rignot 2013) and uncertainty for each ice shelf
    obs_massloss = [35.5, 47.7, 144.9, 101.2, 51.8, 89, 20.7, 155.4, 9.7, 22.5]
    obs_massloss_error = [23, 34, 14, 8, 19, 17, 67, 45, 16, 12]
    # Observed ice shelf melt rates and uncertainty
    obs_ismr = [0.6, 0.1, 4.3, 16.2, 1.7, 3.8, 0.4, 0.3, 0.1, 0.5]
    obs_ismr_error = [0.4, 0.1, 0.4, 1, 0.6, 0.7, 1, 0.1, 0.2, 0.2]
          
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

    # Calculate dA (masked with ice shelf mask) and i and j coordinates
    print 'Analysing grid'
    dA, i, j = calc_grid(file_path)

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

    # Set up array of conversion factors from mass loss to area-averaged melt
    # rate for each ice shelf
    factors = empty(len(names))

    # Process each timestep separately to prevent memory overflow
    for t in range(start_t, size(time)):
        print 'Processing timestep ' + str(t-start_t+1) + ' of ' + str(size(time)-start_t)

        # Read ice shelf melt rate, converting to float128 to prevent overflow
        # during integration
        ismr = ma.asarray(id.variables['m'][t-start_t,:,:], dtype=float128)
        # Convert from m/s to m/y
        ismr = ismr*365.25*24*60*60

        # Loop over ice shelves
        for index in range(len(names)):

            # Mask dA for the current ice shelf (note dA is already masked
            # with the global ice shelf mask, so just restrict the i and j
            # to isolate the given ice shelf)
            if index == len(names)-1:
                # Fimbul-Jelbart region is split into two
                dA_tmp = ma.masked_where(((i < i_min[index]) + (i > i_max[index]) + (j < j_min[index]) + (j > j_max[index]))*((i < i_min[index+1]) + (i > i_max[index+1]) + (j < j_min[index+1]) + (j > j_max[index+1])), dA)
            else:
                dA_tmp = ma.masked_where((i < i_min[index]) + (i > i_max[index]) + (j < j_min[index]) + (j > j_max[index]), dA)

            # Integrate ice shelf melt rate over area to get volume loss
            volumeloss = sum(ismr*dA_tmp)
            # Convert to mass loss in Gt/y
            massloss[index, t] = 1e-12*rho_ice*volumeloss

            if t == start_t:
                # Calculate conversion factor on first timestep
                factors[index] = 1e12/(rho_ice*sum(dA_tmp))

    id.close()

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
        ymin = min(0.95*ismr_low/factors[index], ax1.get_ylim()[0])
        ymax = max(1.05*ismr_high/factors[index], ax1.get_ylim()[1])
        ax1.set_ylim([ymin, ymax])
        # Title and ticks in blue for this side of the plot
        ax1.set_ylabel('Basal Mass Loss (Gt/y)', color='b')
        for t1 in ax1.get_yticklabels():
            t1.set_color('b')
        ax1.set_xlabel('Years')
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
# i- and j-indices.
# Input: file_path = string containing path to ROMS history/averages file
# Output:
# dA = differential of area on the 2D rho-grid, masked with zice
# i,j = i- and j- coordinates on the rho-grid, starting at 1
def calc_grid (file_path):

    # Read grid variables
    id = Dataset(file_path, 'r')
    lon = id.variables['lon_rho'][:,:]
    lat = id.variables['lat_rho'][:,:]
    zice = id.variables['zice'][:,:]
    id.close()

    # Calculate dx and dy in another script
    dx, dy = cartesian_grid_2d(lon, lat)

    # Calculate dA and mask with zice
    dA = ma.masked_where(zice==0, dx*dy)

    # Save dimensions
    num_lat = size(lon, 0)
    num_lon = size(lon, 1)

    # Calculate i- and j-coordinates
    i = tile(arange(1, num_lon+1), (num_lat, 1))
    j = transpose(tile(arange(1, num_lat+1), (num_lon, 1)))

    return dA, i, j


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input('Enter path to ocean history/averages file: ')
    log_path = raw_input('Enter path to log file to save values and/or read previously calculated values: ')

    massloss(file_path, log_path)

