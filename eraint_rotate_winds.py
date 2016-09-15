from netCDF4 import Dataset
from numpy import *

# Rotate the winds in a ROMS-CICE forcing file, interpolated from ERA-Interim as
# in romscice_atm_subdaily.py or romscice_atm_monthly.py, so that they are in
# local x-y space (i.e. positive Uwind points directly to the cell to the right;
# positive Vwind points directly to the cell above) rather than lon-lat space,
# which is slightly rotated.
# Input:
# grid_path = path to ROMS grid file
# eraint_path = path to ERA-Interim forcing file, interpolated to ROMS grid
#               as in romscice_atm_subdaily.py or romscice_atm_monthly.py
def eraint_rotate_winds (grid_path, eraint_path):

    id = Dataset(grid_path, 'r')
    angle = id.variables['angle'][:,:]
    id.close()

    id = Dataset(eraint_path, 'a')
    num_time = id.variables['time'].shape[0]

    # Loop over timesteps
    for t in range(num_time):
        print 'Processing timestep ' + str(t+1) + ' of ' + str(num_time)
        u_lonlat = id.variables['Uwind'][t,:,:]
        v_lonlat = id.variables['Vwind'][t,:,:]
        # Simple rotation matrix
        u_xy = u_lonlat*cos(angle) + v_lonlat*sin(angle)
        v_xy = v_lonlat*cos(angle) - u_lonlat*sin(angle)
        id.variables['Uwind'][t,:,:] = u_xy
        id.variables['Vwind'][t,:,:] = v_xy

    id.close()


# Command-line interface
if __name__ == "__main__":

    grid_path = raw_input("Path to ROMS grid file: ")
    eraint_path = raw_input("Path to ERA-Interim forcing file: ")
    eraint_rotate_winds(grid_path, eraint_path)
