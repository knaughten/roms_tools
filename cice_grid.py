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

    r = 6.371e6
    deg2rad = pi/180.0

    # Open files
    roms = Dataset(roms_grid_name, 'r')
    cice_grid = Dataset(cice_grid_name, 'w')
    cice_kmt = Dataset(cice_kmt_name, 'w')

    # Read variables
    # CICE u-grid corresponds to ROMS psi-grid
    lon_psi = roms.variables['lon_psi'][:,:]
    lat_psi = roms.variables['lat_psi'][:,:]
    # Chop off halo
    ulon = lon_psi[1:,1:]
    ulat = lat_psi[1:,1:]

    # Calculate resolution of each cell in metres
    # Add or subtract 360 from longitude values which wrap around so that
    # longitude increases monotonically from west to east
    num_lat = size(lon_psi, 0)
    num_lon = size(lon_psi, 1)
    i = tile(arange(1,num_lon+1), (num_lat,1))
    index1 = nonzero((i > 1200)*(lon_psi < 100))
    lon_psi[index1] = lon_psi[index1] + 360
    index2 = nonzero((i < 200)*(lon_psi > 300))
    lon_psi[index2] = lon_psi[index2] - 360
    # Now get difference in longitude between psi (u-grid) points
    # Also chop off halo in latitude direction
    dlon = lon_psi[1:,1:] - lon_psi[1:,:-1]
    # Similarly for latitude
    dlat = lat_psi[1:,1:] - lat_psi[:-1,1:]
    # dx on northern edge of tracer cell
    htn = r*cos(ulat*deg2rad)*dlon*deg2rad
    # dy on eastern edge of tracer cell
    hte = r*dlat*deg2rad

    # Calculate dx on u-grid from htn
    dxu = zeros(shape(htn))
    dxu[:,:-1] = 0.5*(htn[:,:-1] + htn[:,1:])
    dxu[:,-1] = 0.5*(htn[:,-1] + htn[:,0])  # Periodic boundary
    # Calculate dx on t-grid from htn
    dxt = zeros(shape(htn))
    dxt[1:,:] = 0.5*(htn[:-1,:] + htn[1:,:])
    dxt[0,:] = 2*dxt[1,:] - dxt[2,:]  # Extrapolate southern boundary
    # Calculate dy on u-grid from hte
    dyu = zeros(shape(hte))
    dyu[:-1,:] = 0.5*(hte[:-1,:] + hte[1:,:])
    dyu[-1,:] = 2*dyu[-2,:] - dyu[-3,:]  # Extrapolate northern boundary
    # Calculate dy on t-grid from hte
    dyt = zeros(shape(hte))
    dyt[:,1:] = 0.5*(hte[:,1:] + hte[:,:-1])
    dyt[:,0] = 0.5*(hte[:,0] + hte[:,-1])  # Periodic boundary

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

    cice_grid.createVariable('HTN', 'f8', ('j', 'i'))
    cice_grid.variables['HTN'].units = 'm'
    cice_grid.variables['HTN'][:,:] = htn

    cice_grid.createVariable('HTE', 'f8', ('j', 'i'))
    cice_grid.variables['HTE'].units = 'm'
    cice_grid.variables['HTE'][:,:] = hte

    cice_grid.createVariable('dxt', 'f8', ('j', 'i'))
    cice_grid.variables['dxt'].units = 'm'
    cice_grid.variables['dxt'][:,:] = dxt

    cice_grid.createVariable('dyt', 'f8', ('j', 'i'))
    cice_grid.variables['dyt'].units = 'm'
    cice_grid.variables['dyt'][:,:] = dyt

    cice_grid.createVariable('dxu', 'f8', ('j', 'i'))
    cice_grid.variables['dxu'].units = 'm'
    cice_grid.variables['dxu'][:,:] = dxu

    cice_grid.createVariable('dyu', 'f8', ('j', 'i'))
    cice_grid.variables['dyu'].units = 'm'
    cice_grid.variables['dyu'][:,:] = dyu

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
