from netCDF4 import Dataset
from numpy import *

# Convert a ROMS grid file to CICE grid and kmt files.
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
    ulon_tmp = roms.variables['lon_psi'][:,:]
    ulat_tmp = roms.variables['lat_psi'][:,:]
    # Convert angle to degrees
    angle = roms.variables['angle'][:]*180/pi
    # Mask out ice shelf cavities for sea ice
    kmt = roms.variables['mask_rho'][:,:] - roms.variables['mask_zice'][:,:]

    i = arange(angle.shape[1])
    j = arange(angle.shape[0])

    # Get one more index in each dimension for lat and lon; linearly extend
    ulon = zeros([size(j), size(i)])
    ulon[:-1,:-1] = ulon_tmp
    # Longitude is a bit tricky because of the mod 360
    for jj in range(size(j)):
        ulon_1back = ulon[jj,-2]
        ulon_2back = ulon[jj,-3]
        if ulon_2back - ulon_1back > 300:
            ulon_2back -= 360
        ulon[jj,-1] = 2*ulon_1back - ulon_2back
    for ii in range(size(i)):
        ulon_1back = ulon[-2,ii]
        ulon_2back = ulon[-3,ii]
        if ulon_2back - ulon_1back < -300:
            ulon_2back += 360
        ulon[-1,ii] = 2*ulon_1back - ulon_2back        
    ulat = zeros([size(j), size(i)])
    ulat[:-1,:-1] = ulat_tmp
    # Latitude is easy
    ulat[-1,:] = 2*ulat[-2,:] - ulat[-3,:]
    ulat[:,-1] = 2*ulat[:,-2] - ulat[:,-3]

    # Write variables
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


# Command-line interface    
if __name__ == "__main__":

    roms_grid_name = raw_input("Path to existing ROMS grid file: ")
    cice_grid_name = raw_input("Path to desired CICE grid file: ")
    cice_kmt_name = raw_input("Path to desired CICE kmt file: ")
    cice_grid(roms_grid_name, cice_grid_name, cice_kmt_name)
