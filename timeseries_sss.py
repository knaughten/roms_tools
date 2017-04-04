from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from cartesian_grid_2d import *

# Calculate and plot timeseries of area-averaged sea surface salinity and
# surface salt flux during a ROMS simulation.
# Input:
# file_path = path to ROMS averages file
# log_path = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
def timeseries_sss (file_path, log_path):

    time = []
    avg_sss = []
    avg_ssflux = []
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
                avg_sss.append(float(line))
            except(ValueError):
                break
        for line in f:
            avg_ssflux.append(float(line))
        f.close()

    print 'Analysing grid'
    id = Dataset(file_path, 'r')
    lon = id.variables['lon_rho'][:-15,1:-1]
    lat = id.variables['lat_rho'][:-15,1:-1]
    zice = id.variables['zice'][:-15,1:-1]
    # Calculate area on the tracer grid
    dx, dy = cartesian_grid_2d(lon, lat)
    dA = ma.masked_where(zice!=0, dx*dy)
    # Read time values and convert from seconds to years
    new_time = id.variables['ocean_time'][:]/(365.25*24*60*60)
    # Concatenate with time values from log file
    for t in range(size(new_time)):
        time.append(new_time[t])

    print 'Reading data'
    # Read surface salinity and salt flux
    # Throw away overlapping periodic boundary and northern sponge layer
    sss = id.variables['salt'][:,-1,:-15,1:-1]
    ssflux = id.variables['ssflux'][:,:-15,1:-1]
    id.close()

    # Build timeseries
    for t in range(size(new_time)):
        avg_sss.append(sum(sss[t,:,:]*dA)/sum(dA))
        avg_ssflux.append(sum(ssflux[t,:,:]*dA)/sum(dA))

    print 'Plotting'
    clf()
    plot(time, avg_sss)
    xlabel('Years')
    ylabel('Average sea surface salinity (psu)')
    grid(True)
    savefig('avg_sss.png')

    clf()
    plot(time, avg_ssflux)
    xlabel('Years')
    ylabel(r'Average surface salt flux (kg/m$^2$/s)')
    grid(True)
    savefig('avg_ssflux.png')

    print 'Saving results to log file'
    f = open(log_path, 'w')
    f.write('Time (years):\n')
    for elm in time:
        f.write(str(elm) + '\n')
    f.write('Average sea surface salinity (psu):\n')
    for elm in avg_sss:
        f.write(str(elm) + '\n')
    f.write('Average surface salt flux (kg/m^2/s):\n')
    for elm in avg_ssflux:
        f.write(str(elm) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to ocean history/averages file: ")
    log_path = raw_input("Path to logfile to save values and/or read previously calculated values: ")
    timeseries_sss(file_path, log_path)
    

    
