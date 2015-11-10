from netCDF4 import Dataset
from numpy import *

# Fill missing timesteps in the ERA-Interim cloud data with the average of
# 1) the nearest previous existing data from the same time of day, and
# 2) the nearest future existing data from the same time of day.
# If 1) doesn't exist, just fill in the timestep with values from 2), and vice
# versa.
# Input: filename = path to ERA-Interim file, regridded and converted to
#                   ROMS-CICE format as in romscice_atm_subdaily.py

def fix_file (filename):

    # Open the file and read the number of timesteps
    id = Dataset(filename, 'a')
    num_steps = id.variables['time'].shape[0]

    # Loop through timesteps
    for t in range(num_steps):
        cloud = id.variables['cloud'][t,:,:]
        if isnan(amin(cloud)):
            # This timestep is missing
            print 'Fixing timestep ' + str(t+1)
            if t < 4:
                # We can only use data from future timesteps.
                # Find the nearest future timestep which is at the same time
                # of day (i.e. a multiple of 4 6-hour timesteps in the future)
                # and where data exists.
                exists = False
                for tt in range(t,num_steps,4):
                    if not isnan(amin(id.variables['cloud'][tt,:,:])):
                        exists = True
                        break
                if exists:
                    cloud = id.variables['cloud'][tt,:,:]
                else:
                    # No such timestep exists
                    print 'Problem with timestep ' + str(t+1)
                    exit
            elif t >= num_steps-4:
                # We can only use data from past timesteps.
                # Similar to above.
                exists = False
                for tt in range(t,-1,-4):
                    if not isnan(amin(id.variables['cloud'][tt,:,:])):
                        exists = True
                        break
                if exists:
                    cloud = id.variables['cloud'][tt,:,:]
                else:
                    print 'Problem with timestep ' + str(t+1)
                    exit
            else:
                # Find one data point on each side.
                for tt1 in range(t,num_steps,4):
                    if not isnan(amin(id.variables['cloud'][tt1,:,:])):
                        exists1 = True
                        break
                for tt2 in range(t,-1,-4):
                    if not isnan(amin(id.variables['cloud'][tt2,:,:])):
                        exists2 = True
                        break
                if exists1 and exists2:
                    # Success on both sides; save the average
                    cloud = 0.5*(id.variables['cloud'][tt1,:,:] + id.variables['cloud'][tt2,:,:])
                elif exists1:
                    # Only the previous data point exists
                    cloud = id.variables['cloud'][tt1,:,:]
                elif exists2:
                    # Only the future data point exists
                    cloud = id.variables['cloud'][tt2,:,:]
                else:
                    print 'Problem with index ' + str(t+1)
                    exit
            # Overwrite the missing timestep with the new data.
            id.variables['cloud'][t,:,:] = cloud

    id.close()
