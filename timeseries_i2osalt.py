from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from cartesian_grid_2d import *

# Calculate and plot timeseries of area-averaged sea ice to ocean salt flux
# during a ROMS-CICE simulation.
# Input:
# file_path = path to CICE history file
# log_path = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
def timeseries_i2osalt (file_path, log_path):

    # Density of freshwater
    rho_fw = 1000.0
    # Density of seawater
    rho_sw = 1025.0
    # Conversion from m/s to cm/day
    mps_to_cmpday = 8.64e6

    time = []
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
    # Read freshwater, salt, and rain fluxes (all scaled by aice) and
    # sea surface salinity
    # Throw away northern sponge layer
    fresh_ai = id.variables['fresh_ai'][:,:-15,:]
    sss = id.variables['sss'][:,:-15,:]
    rain_ai = id.variables['rain_ai'][:,:-15,:]
    fsalt_ai = id.variables['fsalt_ai'][:,:-15,:]
    id.close()

    # Build timeseries
    for t in range(size(new_time)):
        # Merge CICE's freshwater and salt fluxes as in set_vbc.F
        # Subtract rain because we don't care about that
        # Convert to kg/m^2/s
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


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to CICE history file: ")
    log_path = raw_input("Path to logfile to save values and/or read previously calculated values: ")
    timeseries_i2osalt(file_path, log_path)
    
