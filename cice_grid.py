from netCDF4 import Dataset
from numpy import *

# Convert a ROMS grid file to CICE grid and kmt files. Note the CICE grid
# will have two less points on either dimension because it only wants
# computational cells.
# Input:
# roms_grid_name = path to existing ROMS grid file
# cice_grid_name = path to desired CICE grid file
# cice_kmt_name = path to desired CICE kmt file

def cice_grid (roms_grid_name, cice_grid_name, cice_kmt_name):

    # Open files
    roms = Dataset(roms_grid_name, 'r')
    cice_grid = Dataset(cice_grid_name, 'w')
    cice_kmt = Dataset(cice_kmt_name, 'w')

    # Read variables
    # CICE u-grid corresponds to ROMS psi-grid
    ulon = roms.variables['lon_psi'][1:,1:]
    ulat = roms.variables['lat_psi'][1:,1:]
    # Convert angle (on shared tracer grid) to degrees
    angle = roms.variables['angle'][1:-1,1:-1]*180/pi
    # Mask out ice shelf cavities for sea ice
    kmt = roms.variables['mask_rho'][1:-1,1:-1] - roms.variables['mask_zice'][1:-1,1:-1]

    num_lon = size(ulon,1)
    num_lat = size(ulon,0)

    # Write variables
    cice_grid.createDimension('i', num_lon)
    cice_grid.createDimension('j', num_lat)

    cice_grid.createVariable('ulon', 'f8', ('j', 'i'))
    cice_grid.variables['ulon'].units = 'degree_east'
    cice_grid.variables['ulon'][:,:] = ulon

    cice_grid.createVariable('ulat', 'f8', ('j', 'i'))
    cice_grid.variables['ulat'].units = 'degree_north'
    cice_grid.variables['ulat'][:,:] = ulat

    cice_grid.createVariable('angle', 'f8', ('j', 'i'))
    cice_grid.variables['angle'].units = 'radians'
    cice_grid.variables['angle'][:,:] = angle

    cice_kmt.createDimension('i', num_lon)
    cice_kmt.createDimension('j', num_lat)

    cice_kmt.createVariable('kmt', 'f8', ('j', 'i'))
    cice_kmt.variables['kmt'].units = '1'
    cice_kmt.variables['kmt'][:,:] = kmt

    roms.close()
    cice_grid.close()
    cice_kmt.close()


# Command-line interface    
if __name__ == "__main__":

    roms_grid_name = raw_input("Path to existing ROMS grid file: ")
    cice_grid_name = raw_input("Path to desired CICE grid file: ")
    cice_kmt_name = raw_input("Path to desired CICE kmt file: ")
    cice_grid(roms_grid_name, cice_grid_name, cice_kmt_name)
