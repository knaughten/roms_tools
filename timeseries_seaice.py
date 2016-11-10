from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from cartesian_grid_2d import *

# Calculate and plot timeseries of total sea ice area and volume during a
# ROMS-CICE simulation.
# Input:
# file_path = path to CICE history file
# log_path = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
def timeseries_seaice (file_path, log_path):

    time = []
    total_area = []
    total_volume = []
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
            try:
                total_area.append(float(line))
            except(ValueError):
                break
        for line in f:
            total_volume.append(float(line))
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
    # Read sea ice concentration and height
    # Throw away overlapping periodic boundary and northern sponge layer
    aice = id.variables['aice'][:,:-15,:]
    hi = id.variables['hi'][:,:-15,:]
    id.close()

    print 'Setting up arrays'
    # Remove masks and fill with zeros (was having weird masking issues here)
    aice_nomask = aice.data
    aice_nomask[aice.mask] = 0.0
    hi_nomask = hi.data
    hi_nomask[hi.mask] = 0.0

    # Build timeseries
    for t in range(size(new_time)):
        # Integrate area and convert to million km^2
        total_area.append(sum(aice_nomask[t,:,:]*dA)*1e-12)
        # Integrate volume and convert to million km^3
        total_volume.append(sum(aice_nomask[t,:,:]*hi_nomask[t,:,:]*dA)*1e-12)

    print 'Plotting total sea ice area'
    clf()
    plot(time, total_area)
    xlabel('Years')
    ylabel(r'Total Sea Ice Area (million km$^2$)')
    grid(True)
    savefig('seaice_area.png')

    print 'Plotting total sea ice volume'
    clf()
    plot(time, total_volume)
    xlabel('Years')
    ylabel(r'Total Sea Ice Volume (million km$^3$)')
    grid(True)
    savefig('seaice_volume.png')

    print 'Saving results to log file'
    f = open(log_path, 'w')
    f.write('Time (years):\n')
    for elm in time:
        f.write(str(elm) + '\n')
    f.write('Total Sea Ice Area (million km^2):\n')
    for elm in total_area:
        f.write(str(elm) + '\n')
    f.write('Total Sea Ice Volume (million km^3):\n')
    for elm in total_volume:
        f.write(str(elm) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to CICE history file: ")
    log_path = raw_input("Path to logfile to save values and/or read previously calculated values: ")
    timeseries_seaice(file_path, log_path)

    

    
