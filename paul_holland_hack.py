from netCDF4 import Dataset
from numpy import *
from cartesian_grid_2d import *

def paul_holland_hack (grid_file):

    total_fw = 1500 # Gt/y
    nbdry = -60 # Apply the freshwater evenly south of here
    sec_per_year = 365.25*24*60*60

    id = Dataset(grid_file, 'r')
    lat = id.variables['lat_rho'][:,:]
    lon = id.variables['lon_rho'][:,:]
    mask_zice = id.variables['mask_zice'][:,:]
    mask_rho = id.variables['mask_rho'][:,:]
    id.close()
    mask = mask_rho - mask_zice

    dx, dy = cartesian_grid_2d(lon, lat)
    ocn_flag = mask == 1
    loc_flag = lat < nbdry
    total_area = sum(dx*dy*ocn_flag)
    print 'Total area = ' + str(total_area) + ' m^2'
    target_area = sum(dx*dy*ocn_flag*loc_flag)
    print 'Area south of 60S = ' + str(target_area) + ' m^2'
    # Multiply by 1e12 to convert from Gt/y to kg/y
    # Divide by sec_per_year to convert from kg/y to kg/s
    # Divide by target area to get kg/m^2/s
    fw_flux = total_fw*1e12/target_area/sec_per_year
    print 'Freshwater flux to add = ' + str(fw_flux) + 'kg/m^2/s'


if __name__ == "__main__":

    grid_file = raw_input("Path to ROMS grid file: ")
    paul_holland_hack(grid_file)
