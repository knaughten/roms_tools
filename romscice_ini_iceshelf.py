from netCDF4 import Dataset
from numpy import *

# Given a ROMS initial condition file from ECCO2 data, overwrite the salinity
# in ice shelf cavities (which aren't in ECCO2) with salinity from Ben's
# original initial conditions file (not really sure where this data is from,
# but it's more realistic than filling the cavities with constant salinity, and
# certainly more stable).

# Paths (other users need to edit these)
ini_file = '../data/caisom001_ini_1995.nc' # Initial conditions file to edit
is_ini_file = '../data/caisom001_ini_iceshelf.nc' # Ben's file
grid_file = '../apps/common/grid/caisom001_OneQuartergrd.nc' # ROMS grid file

# Read zice and make a boolean mask
grid_id = Dataset(grid_file, 'r')
mask_zice = grid_id.variables['mask_zice'][:,:]
zice_index = mask_zice > 0.9 # better for floating-point than ==1
grid_id.close()

ini_id = Dataset(ini_file, 'a')
is_id = Dataset(is_ini_file, 'r')
num_depth = ini_id.variables['salt'].shape[1]

for k in range(num_depth):
    # For each depth level, replace the salinity values in ice shelf cavities
    print 'Processing depth level ' + str(k+1) + ' of ' + str(num_depth)
    ini_salt = ini_id.variables['salt'][0,k,:,:]
    is_salt = is_id.variables['salt'][0,k,:,:]
    ini_salt[zice_index] = is_salt[zice_index]
    ini_id.variables['salt'][0,k,:,:] = ini_salt
ini_id.close()
is_id.close()



