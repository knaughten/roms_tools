from netCDF4 import Dataset
from numpy import *

def cice_grid (roms_grid_name, cice_grid_name, cice_kmt_name):

    roms = Dataset(roms_grid_name, 'r')
    cice_grid = Dataset(cice_grid_name, 'w')
    cice_kmt = Dataset(cice_kmt_name, 'w')

    ulon = roms.variables['lon_rho'][:,:]
    ulat = roms.variables['lat_rho'][:,:]
    angle = roms.variables['angle'][:]
    kmt = roms.variables['mask_rho'][:,:] - roms.variables['mask_zice'][:,:]

    i = arange(angle.shape[1])
    j = arange(angle.shape[0])

    cice_grid.createDimension('i', size(i))
    cice_grid.createDimension('j', size(j))

    cice_grid.createVariable('ulon', 'f8', ('j', 'i'))
    cice_grid.variables['ulon'].units = 'degree_east'
    cice_grid.variables['ulon'][:,:] = ulon

    cice_grid.createVariable('ulat', 'f8', ('j', 'i'))
    cice_grid.variables['ulat'].units = 'degree_north'
    cice_grid.variables['ulat'][:,:] = ulat

    cice_grid.createVariable('angle', 'f8', ('j', 'i'))
    cice_grid.variables['angle'].units = 'radians'
    cice_grid.variables['angle'][:,:] = angle

    cice_kmt.createDimension('i', size(i))
    cice_kmt.createDimension('j', size(j))

    cice_kmt.createVariable('kmt', 'f8', ('j', 'i'))
    cice_kmt.variables['kmt'].units = '1'
    cice_kmt.variables['kmt'][:,:] = kmt

    roms.close()
    cice_grid.close()
    cice_kmt.close()

    
if __name__ == "__main__":

    roms_grid_name = raw_input("Path to existing ROMS grid file: ")
    cice_grid_name = raw_input("Path to desired CICE grid file: ")
    cice_kmt_name = raw_input("Path to desired CICE kmt file: ")
    cice_grid(roms_grid_name, cice_grid_name, cice_kmt_name)
