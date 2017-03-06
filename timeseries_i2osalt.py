from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from cartesian_grid_2d import *

def timeseries_i2osalt (file_path, log_path):

    rho_fw = 1000.0
    rho_sw = 1026.0
    mps_to_cmpday = 8.64e6

    time = []
    avg_ssflux = []
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
            avg_ssflux.append(float(line))
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
    fresh_ai = id.variables['fresh_ai'][:,:-15,:]
    sss = id.variables['sss'][:,:-15,:]
    rain_ai = id.variables['rain_ai'][:,:-15,:]
    fsalt_ai = id.variables['fsalt_ai'][:,:-15,:]
    id.close()

    for t in range(size(new_time)):
        avg_ssflux.append(sum(-1/rho_fw*((fresh_ai[t,:,:]-rain_ai[t,:,:])*sss[t,:,:]*rho_sw/mps_to_cmpday - fsalt_ai[t,:,:]*1e3)*dA)/sum(dA))

    print 'Plotting'
    clf()
    plot(time, avg_ssflux)
    xlabel('Years')
    ylabel(r'Average sea ice to ocean salt flux (kg/m$^2$/s)')
    grid(True)
    savefig('avg_i2osalt.png')

    print 'Saving results to log file'
    f = open(log_path, 'w')
    f.write('Time (years):\n')
    for elm in time:
        f.write(str(elm) + '\n')
    f.write('Average sea ice to ocean salt flux (kg/m^2/s):\n')
    for elm in avg_ssflux:
        f.write(str(elm) + '\n')
    f.close()


if __name__ == "__main__":

    file_path = raw_input("Path to CICE history file: ")
    log_path = raw_input("Path to logfile to save values and/or read previously calculated values: ")
    timeseries_i2osalt(file_path, log_path)
    
