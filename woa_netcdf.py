from netCDF4 import Dataset
from numpy import *

# Read the World Ocean Atlas temperature and salinity climatology from a text
# file (FESOM input format), and save in a NetCDF file for easier use.

# File paths
in_file = '/short/y99/kaa561/FESOM/annual_woa01_ts.out'
out_file = '/short/y99/kaa561/FEsOM/woa01_ts.nc'

print 'Reading text file'
file = open(in_file, 'r')

# First line contains the number of longitude, latitude, and depth values
sizes = (file.readline()).split()
num_lon = int(sizes[0])
num_lat = int(sizes[1])
num_depth = int(sizes[2])

# Read longitude values
lon = []
while True:
    # Each row has multiple values, so split based on white space
    lon_vals = (file.readline()).split()
    for elm in lon_vals:
        lon.append(float(elm))
    if len(lon) == num_lon:
        # Finished with longitude
        break
# Repeat for latitude
lat = []
while True:
    lat_vals = (file.readline()).split()
    for elm in lat_vals:
        lat.append(float(elm))
    if len(lat) == num_lat:
        break
# Repeat for depth
depth = []
while True:
    depth_vals = (file.readline()).split()
    for elm in depth_vals:
        # Convert from negative to positive depth
        depth.append(-1*float(elm))
    if len(depth) == num_depth:
        break
# Read temperature values (save in one long 1D array for now)
temp = []
while True:
    temp_vals = (file.readline()).split()
    for elm in temp_vals:
        temp.append(float(elm))
    if len(temp) == num_lon*num_lat*num_depth:
        break
# Repeat for salinity
salt = []
while True:
    salt_vals = (file.readline()).split()
    for elm in salt_vals:
        salt.append(float(elm))
    if len(salt) == num_lon*num_lat*num_depth:
        break
file.close()

# Copy contents of the long 1D temp and salt arrays into 3D arrays
# (depth x latitude x longitude)
print 'Reshaping temperature and salinity arrays'
temp_3d = zeros((num_depth, num_lat, num_lon))
salt_3d = zeros((num_depth, num_lat, num_lon))
posn = 0
for i in range(num_lon):
    for j in range(num_lat):
        for k in range(num_depth):
            temp_3d[k,j,i] = temp[posn]
            salt_3d[k,j,i] = salt[posn]
            posn = posn+1

# Output to NetCDF file
print 'Writing NetCDF file'
id = Dataset(out_file, 'w')
id.createDimension('longitude', num_lon)
id.createDimension('latitude', num_lat)
id.createDimension('depth', num_depth)
id.createVariable('longitude', 'f8', ('longitude'))
id.variables['longitude'].units = 'degrees'
id.variables['longitude'][:] = lon
id.createVariable('latitude', 'f8', ('latitude'))
id.variables['latitude'].units = 'degrees'
id.variables['latitude'][:] = lat
id.createVariable('depth', 'f8', ('depth'))
id.variables['depth'].units = 'metres'
id.variables['depth'][:] = depth
id.createVariable('temp', 'f8', ('depth', 'latitude', 'longitude'))
id.variables['temp'].units = 'C'
id.variables['temp'][:,:,:] = temp_3d
id.createVariable('salt', 'f8', ('depth', 'latitude', 'longitude'))
id.variables['salt'].units = 'psu'
id.variables['salt'][:,:,:] = salt_3d
id.close()

