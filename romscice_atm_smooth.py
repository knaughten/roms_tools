from netCDF4 import Dataset
from numpy import *

# Smooth the given variable in the given NetCDF file (dimensions time x latitude
# x longitude) using a moving average in the time dimension. Towards the
# endpoints of the time dimension, wrap the moving average around, treating the
# input file as a periodic dataset (e.g. 1 year of forcing to be repeated for
# a spinup).
# Input:
# filename = path to NetCDF file to modify
# var = string containing name of variable in filename to smooth; note it
#       must have dimensions time x latitude x longitude
# interval = number of time records to use for moving average
#            If interval is odd, it will be interpreted as a moving average
#            of <interval> timesteps centered on the current timestep, e.g.
#            interval = 11 means average over the current timestep, 5 before,
#            and 5 after.
#            If interval is even, it will be interpreted as a moving average
#            of the current timestep and the <interval> nearest timesteps,
#            e.g. interval = 10 gives the same behaviour as interval = 11
#            described above.

def smooth_var (filename, var, interval):

    # Check if interval is odd or even
    if interval % 2 == 0:
        # Add 1 to make it odd
        interval = interval+1
    # Calculate the number of timesteps on either side of the centre timestep
    margin = (interval-1)/2

    # Open file and save dimensions
    id = Dataset(filename, 'a')
    dim = id.variables[var].shape
    num_time = dim[0]
    num_lat = dim[1]
    num_lon = dim[2]
    # Set up array to store <interval> timesteps for moving average
    data = zeros((interval, num_lat, num_lon))

    # Fill in the data array for the first timestep
    print 'Processing first timestep'
    for n in range(interval):
        # n-margin will be negative for the first margin values
        # (wrapping around to the end of the file, e.g. index -1 means
        # the last index, -2 means second last, etc), then 0, then
        # positive for the next margin values.
        data[n,:,:] = id.variables[var][n-margin,:,:]
    # Overwrite the first timestep to be the time-average of data
    id.variables[var][0,:,:] = mean(data, axis=0)        

    # Loop over remaining timesteps
    for t in range(1, num_time):
        print 'Processing timestep ' + str(t+1) + ' of ' + str(num_time)
        # Shift data to the left by 1
        data[:-1,:,:] = data[1:,:,:]
        # Now overwrite the last index
        if t+margin >= num_time:
            # Nearing the end of the file; loop back around to the beginning
            data[-1,:,:] = id.variables[var][t+margin-num_time,:,:]
        else:
            # Easy, just moving ahead by 1 timestep with no wrapping around
            data[-1,:,:] = id.variables[var][t+margin,:,:]
        # Overwrite this timestep to be the time-average of data
        id.variables[var][t,:,:] = mean(data, axis=0)

    id.close()        


# User-defined values here
if __name__ == "__main__":

    an_file = '/short/m68/kaa561/ROMS-CICE-MCT/data/ERA_Interim/subdaily/30day_smoothed/AN_1995_subdaily.nc'
    fc_file = '/short/m68/kaa561/ROMS-CICE-MCT/data/ERA_Interim/subdaily/30day_smoothed/FC_1995_subdaily.nc'
    # Names of variables to process in each file
    an_var = ['Pair', 'Tair', 'Qair', 'cloud', 'Uwind', 'Vwind']
    fc_var = ['rain', 'snow']
    sdays = 30  # 30-day moving average

    # Loop over each listed variable in each file
    for var in an_var:
        print 'Computing ' + str(sdays) + '-day moving average of ' + var
        # This file has 6-hourly values
        smooth_var(an_file, var, sdays*4)
    for var in fc_var:
        print 'Computing ' + str(sdays) + '-day moving average of ' + var
        # This value has 12-hourly values
        smooth_var(fc_file, var, sdays*2)
