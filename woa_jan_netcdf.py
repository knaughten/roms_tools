from netCDF4 import Dataset
from numpy import *

def woa_jan_netcdf ():

    temp_file = '/short/m68/kaa561/woa_jan/woa13_decav_t01mn01v2.csv'
    salt_file = '/short/m68/kaa561/woa_jan/woa13_decav_t01mn01v2.csv'
    out_file = '/short/m68/kaa561/woa_jan/woa_jan.nc'

    lat = arange(-77.5, 87.5+1, 1)
    lon = arange(-179.5, 179.5+1, 1)
    num_lat = size(lat)
    num_lon = size(lon)

    file = open(temp_file, 'r')
    file.readline()
    depth_data = file.readline()
    index = depth_data.index('0')
    depth_data = depth_data[index:]
    depth_vals = depth_data.split(',')
    depth = []
    for elm in depth_vals:
        depth.append(float(elm))
    depth = array(depth)
    num_depth = size(depth)

    temp_data = zeros([num_depth, num_lat, num_lon])
    temp_data[:,:,:] = -999
    for line in file:
        data = line.split(',')
        lat0 = float(data[0])
        lon0 = float(data[1])
        j = where(lat==lat0)[0][0]
        i = where(lon==lon0)[0][0]
        data = data[2:]
        k = 0
        for elm in data:
            if elm not in ['', '\n']:
                temp_data[k,j,i] = float(elm)
            k += 1
    file.close()
    temp_data = ma.masked_where(temp_data==-999, temp_data)

    salt_data = zeros([num_depth, num_lat, num_lon])
    salt_data[:,:,:] = -999
    file = open(salt_file, 'r')
    file.readline()
    file.readline()
    for line in file:
        data = line.split(',')
        lat0 = float(data[0])
        lon0 = float(data[1])
        j = where(lat==lat0)[0][0]
        i = where(lon==lon0)[0][0]
        data = data[2:]
        k = 0
        for elm in data:
            if elm not in ['', '\n']:
                salt_data[k,j,i] = float(elm)
            k += 1
    file.close()
    salt_data = ma.masked_where(salt_data==-999, salt_data)

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
    id.variables['temp'][:,:,:] = temp_data
    id.createVariable('salt', 'f8', ('depth', 'latitude', 'longitude'))
    id.variables['salt'].units = 'psu'
    id.variables['salt'][:,:,:] = salt_data
    id.close()


if __name__ == "__main__":

    woa_jan_netcdf()
        

    

    

    
