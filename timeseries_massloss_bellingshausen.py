from netCDF4 import Dataset
from numpy import *
from timeseries_massloss import calc_grid

# Calculate timeseries of basal mass loss for the Wilkins, Stange, and
# George VI ice shelves separately. Output to a log file.
# Input:
# directory = path to directory containing containing ocean averages files,
#             assuming 5-day averages
# start_index, end_index = integers containing range of files to process. For
#                          example, start_index=1 and end_index=102 will
#                          process files ocean_avg_0001.nc through
#                          ocean_avg_0102.nc.
# log_path = path to desired output log file
def timeseries_massloss_bellingshausen (directory, start_index, end_index, log_path):

    # Name of each ice shelf
    names = ['Wilkins Ice Shelf', 'Stange Ice Shelf', 'George VI Ice Shelf']
    # Limits on longitude and latitude for each ice shelf
    # These depend on the source geometry, in this case RTopo 1.05
    # Note there is one extra index at the end of each array; this is because
    # the George VI region is awkwardly shaped and split into 2 (and hence why
    # all 3 ice shelves are lumped together in timeseries_massloss.py)
    lon_min = [-75, -79, -74.5, -69.5]
    lon_max = [-69, -74.5, -67, -66]
    lat_min = [-71.5, -73.8, -73.8, -72.6]
    lat_max = [-69, -72.6, -72.6, -70]
    # Density of ice in kg/m^3
    rho_ice = 916

    print 'Analysing grid'
    # Using the first file, calculate dA (masked with ice shelf mask) and
    # lon and lat coordinates
    dA, lon, lat = calc_grid(directory + index_to_file(start_index))

    print 'Setting up arrays'
    # Set up array of masked area values for each ice shelf
    dA_masked = ma.empty([len(names), size(dA,0), size(dA,1)])
    for shelf in range(len(names)):
        # Mask dA for the current ice shelf (note dA is already masked
        # with the global ice shelf mask, so just restrict the indices
        # to isolate the given ice shelf)
        if shelf == len(names)-1:
            # George VI region is split into two
            dA_tmp = ma.masked_where(((lon < lon_min[shelf]) + (lon > lon_max[shelf]) + (lat < lat_min[shelf]) + (lat > lat_max[shelf]))*((lon < lon_min[shelf+1]) + (lon > lon_max[shelf+1]) + (lat < lat_min[shelf+1]) + (lat > lat_max[shelf+1])), dA)
        else:
            dA_tmp = ma.masked_where((lon < lon_min[shelf]) + (lon > lon_max[shelf]) + (lat < lat_min[shelf]) + (lat > lat_max[shelf]), dA)
        dA_masked[shelf,:,:] = dA_tmp[:,:]
    # Figure out how many time indices in total there are
    total_time = 0
    for index in range(start_index, end_index+1):
        id = Dataset(directory + index_to_file(index), 'r')
        total_time += id.variables['ocean_time'].size
        id.close()
    # Set up array of mass loss values for each ice shelf
    massloss = empty([len(names), total_time])
    # Also time values
    time = empty(total_time)

    t_posn = 0
    for index in range(start_index, end_index+1):
        file = directory + index_to_file(index)
        print 'Processing ' + file
        id = Dataset(file, 'r')
        # Read melt rate and convert from m/s to m/y
        ismr = id.variables['m'][:,:-15,1:-1]*365.25*24*60*60
        # Read time and convert from seconds to years
        curr_time = id.variables['ocean_time'][:]/(365.25*24*60*60)
        id.close()
        # Save the time values
        num_time = size(curr_time)
        time[t_posn:t_posn+num_time] = curr_time
        for t in range(num_time):
            for shelf in range(len(names)):
                # Integrate ice shelf melt rate over area to get volume loss
                volumeloss = sum(ismr[t,:,:]*dA_masked[shelf,:,:])
                # Convert to mass loss in Gt/y
                massloss[shelf, t_posn+t] = 1e-12*rho_ice*volumeloss
        t_posn += num_time

    print 'Saving results to log file'
    f = open(log_path, 'w')
    f.write('Time (years):\n')
    for t in range(size(time)):
        f.write(str(time[t]) + '\n')
    for shelf in range(len(names)):
        f.write(names[shelf] + ' Basal Mass Loss\n')
        for t in range(size(time)):
            f.write(str(massloss[shelf, t]) + '\n')
    f.close()
    

# Given an integer, return the filename for the corresponding ocean averages
# file. For example, index_to_file(1) = 'ocean_avg_0001.nc', and
# index_to_file(95) = 'ocean_avg_0095.nc'.
def index_to_file (index):

    if index < 10:
        return 'ocean_avg_000' + str(index) + '.nc'
    elif index < 100:
        return 'ocean_avg_00' + str(index) + '.nc'
    elif index < 1000:
        return 'ocean_avg_0' + str(index) + '.nc'
    else:
        return 'ocean_avg_' + str(index) + '.nc'


# Command-line interface
if __name__ == "__main__":

    directory = raw_input("Path to ROMS output directory: ")
    start_index = int(raw_input("Index of first ocean averages file (e.g. ocean_avg_0001.nc has index 1): "))
    end_index = int(raw_input("Index of last ocean averages file: "))
    log_path = raw_input("Path to desired logfile: ")
    timeseries_massloss_bellingshausen(directory, start_index, end_index, log_path)


    
