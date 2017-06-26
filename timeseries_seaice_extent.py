from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from cartesian_grid_2d import *

# Calculate and plot timeseries of sea ice extent (area of ice with
# concentration >= 15%) during a ROMS-CICE simulation.
# Input:
# file_path = path to CICE history file
# log_path = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
def timeseries_seaice_extent (file_path, log_path):

    time = []
    extent = []
    # Check if the log file exists
    if exists(log_path):
        print 'Reading previously calculated values'
        f = open(log_path, 'r')
        # Skip first line (header for time array)
        f.readline()
        for line in f:
            try:
                time.append(float(line))
            except(ValueError):
                # Reached the header for the next variable
                break
        for line in f:
            extent.append(float(line))
        f.close()

    print 'Analysing grid'
    id = Dataset(file_path, 'r')
    lon = id.variables['TLON'][:-15,:]
    lat = id.variables['TLAT'][:-15,:]
    # Calculate area on the tracer grid
    dx, dy = cartesian_grid_2d(lon, lat)
    dA = dx*dy
    # Read time values and convert from days to years
    new_time = id.variables['time'][:]/365.25
    # Concatenate with time values from log file
    for t in range(size(new_time)):
        time.append(new_time[t])

    print 'Reading data'
    # Read sea ice concentration
    # Throw away northern sponge layer
    aice = id.variables['aice'][:,:-15,:]
    id.close()
    # Select cells with concentration >= 15%
    flag = aice >= 0.15

    print 'Building timeseries'
    for t in range(size(new_time)):
        # Integrate extent and convert to million km^2
        extent.append(sum(flag[t,:,:]*dA)*1e-12)

    print 'Plotting'
    clf()
    plot(time, extent)
    xlabel('Years')
    ylabel(r'Sea Ice Extent (million km$^2$)')
    grid(True)
    savefig('seaice_extent.png')

    print 'Saving results to log file'
    f = open(log_path, 'w')
    f.write('Time (years):\n')
    for elm in time:
        f.write(str(elm) + '\n')
    f.write('Sea Ice Extent (million km^2):\n')
    for elm in extent:
        f.write(str(elm) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to CICE history file: ")
    log_path = raw_input("Path to logfile to save values and/or read previously calculated values: ")
    timeseries_seaice_extent(file_path, log_path)

    
