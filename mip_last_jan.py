from netCDF4 import Dataset
from numpy import *
from monthly_avg_roms import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from monthly_avg import *

def mip_last_jan (roms_in_file, roms_out_file, fesom_in_file, fesom_out_file):

    id = Dataset(roms_in_file, 'r')
    lon = id.variables['lon_rho'][:,:]
    lat = id.variables['lat_rho'][:,:]
    num_depth = id.variables['temp'].shape[1]
    id.close()
    num_lat = size(lat,0)
    num_lon = size(lon,1)

    print 'Calculating monthly averages for ROMS'
    roms_temp = monthly_avg_roms(roms_in_file, 'temp', [num_depth, num_lat, num_lon], 0)
    roms_salt = monthly_avg_roms(roms_in_file, 'salt', [num_depth, num_lat, num_lon], 0)

    print 'Writing ' + roms_out_file
    id = Dataset(roms_out_file, 'w')
    id.createDimension('xi_rho', size(lon,1))
    id.createDimension('eta_rho', size(lon,0))
    id.createDimension('s_rho', num_depth)    
    id.createDimension('time', None)
    id.createVariable('lon_rho', 'f8', ('eta_rho', 'xi_rho'))
    id.variables['lon_rho'].long_name = 'longitude of rho-points'
    id.variables['lon_rho'].units = 'degree_east'
    id.variables['lon_rho'][:,:] = lon
    id.createVariable('lat_rho', 'f8', ('eta_rho', 'xi_rho'))
    id.variables['lat_rho'].long_name = 'latitude of rho-points'
    id.variables['lat_rho'].units = 'degree_north'
    id.variables['lat_rho'][:,:] = lat
    id.createVariable('time', 'f8', ('time'))
    id.variables['time'].units = 'month'
    id.variables['time'].description = 'DJF, MAM, JJA, SON'
    id.variables['time'][0] = 1
    id.createVariable('temp', 'f8', ('time', 's_rho', 'eta_rho', 'xi_rho'))
    id.variables['temp'].units = 'degC'
    id.variables['temp'][0,:,:,:] = roms_temp
    id.createVariable('salt', 'f8', ('time', 's_rho', 'eta_rho', 'xi_rho'))
    id.variables['salt'].units = 'psu'
    id.variables['salt'][0,:,:,:] = roms_salt
    id.close()

    print 'Calculating monthly averages for FESOM'
    fesom_temp = monthly_avg(fesom_in_file, 'temp', 0)
    fesom_salt = monthly_avg(fesom_in_file, 'salt', 0)

    print 'Writing ' + fesom_out_file
    id = Dataset(fesom_out_file, 'w')
    id.createDimension('nodes_3d', size(fesom_temp))
    id.createDimension('T', None)
    id.createVariable('temp', 'f8', ('T', 'nodes_3d'))
    id.variables['temp'].description = 'mean potential temperature'
    id.variables['temp'].units = 'degC'
    id.variables['temp'][0,:] = fesom_temp
    id.createVariable('salt', 'f8', ('T', 'nodes_3d'))
    id.variables['salt'].description = 'mean salinity'
    id.variables['salt'].units = 'psu'
    id.variables['salt'][0,:] = fesom_salt
    id.close()
    

# Command-line interface
if __name__ == "__main__":

    roms_in_file = raw_input("Path to ROMS ocean averages file containing January 2016: ")
    roms_out_file = raw_input("Path to desired ROMS file to save monthly average: ")
    fesom_in_file = raw_input("Path to FESOM oce.mean file for 2016: ")
    fesom_out_file = raw_input("Path to desired FESOM file to save monthly average: ")
    mip_last_jan(roms_in_file, roms_out_file, fesom_in_file, fesom_out_file)
