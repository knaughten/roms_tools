from netCDF4 import Dataset
from numpy import *
from interp_era2roms import *

# Read one year of ERA-Interim evaporation data (12-hourly), interpolate to the
# ROMS grid, and add to the existing FC forcing file created using
# romscice_atm_subdaily.nc.
# Input: year = integer containing the year to process
#        count = time record in the given year to start with

# This script only processes 50 12-hour timesteps at once to prevent memory 
# overflow, and is designed to be called by a self-submitting batch script. 
# See era_evap.job for an example.

def convert_file (year, count):

    # Make sure input arguments are integers (sometimes the batch script likes
    # to pass them as strings)
    year = int(year)
    count = int(count)

    # Paths of ROMS grid file, input ERA-Interim files, and output ROMS-CICE
    # files; other users will need to change these
    grid_file = '/short/m68/kaa561/metroms_iceshelf/apps/common/grid/circ30S_quarterdegree.nc'
    input_evap_file = '/short/m68/kaa561/metroms_iceshelf/data/subdaily_originals/ER_' + str(year)+ '_subdaily_orig.nc'
    output_evap_file = '/short/m68/kaa561/metroms_iceshelf/data/ERA_Interim/FC_' + str(year) + '_subdaily.nc'
    logfile = str(year) + '.log'

    if count == 0:
        log = open(logfile, 'w')
    else:
        log = open(logfile, 'a')
    log.write('Reading grids\n')
    log.close()

    # Read ROMS latitude and longitude
    grid_fid = Dataset(grid_file, 'r')
    lon_roms = grid_fid.variables['lon_rho'][:,:]
    lat_roms = grid_fid.variables['lat_rho'][:,:]
    grid_fid.close()
    num_lon = size(lon_roms, 1)
    num_lat = size(lon_roms, 0)

    # Open input ER file and read time values
    i_fid = Dataset(input_evap_file, 'r')
    evap_time = i_fid.variables['time'][:] # hours since 1900-01-01 00:00:0.0
    # Convert time units
    evap_time = evap_time/24.0 # days since 1900-01-01 00:00:0.0
    evap_time = evap_time - 92*365 - 22 # days since 1992-01-01 00:00:0.0; note that there were 22 leap years between 1900 and 1992
    evap_time = evap_time - 0.5 # switch from evaporation over the preceding 12 hours to evaporation over the following 12 hours; this is easier for ROMS
    # Also read ERA-Interim latitude and longitude
    lon_era = i_fid.variables['longitude'][:]
    lat_era = i_fid.variables['latitude'][:]
    i_fid.close()

    # Define the variable in the output NetCDF file on the first timestep
    if count == 0:
        log = open(logfile, 'a')
        log.write('Adding evaporation to ' + output_evap_file + '\n')
        log.close()

        o_fid = Dataset(output_evap_file, 'a')
        o_fid.createVariable('evaporation', 'f8', ('time', 'eta_rho', 'xi_rho'))
        o_fid.variables['evaporation'].long_name = 'evaporation rate'
        o_fid.variables['evaporation'].units = 'm_per_12hr'
        o_fid.close()

    log = open(logfile, 'a')
    log.write('Processing 12-hourly data\n')
    log.close()

    for t in range(count/2, (count+100)/2):
        if t >= size(evap_time):
            break
        o_fid = Dataset(output_evap_file, 'a')
        log = open(logfile, 'a')
        log.write('Processing record ' + str(t+1) + ' of ' + str(size(evap_time)) + '\n')
        log.close()
        # Write the current time value to output FC file
        o_fid.variables['time'][t] = evap_time[t]
        # Read data for this timestep
        i_fid = Dataset(input_evap_file, 'r')
        e = -1*transpose(i_fid.variables['e'][t,:,:])
        i_fid.close()
        # Interpolate to ROMS grid and write to output FC file
        evap = interp_era2roms(e, lon_era, lat_era, lon_roms, lat_roms)
        o_fid.variables['evaporation'][t,:,:] = evap
        o_fid.close()

    log = open(logfile, 'a')
    log.write('Finished\n')
    log.close()


        
    
    
        
        
    
    
    
