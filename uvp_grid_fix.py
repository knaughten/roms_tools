from netCDF4 import Dataset
from numpy import *

# The Matlab scripts I use to generate the ROMS grid assume that the u-grid,
# v-grid, and psi-grid points have the same locations as the rho-grid, and they
# just chop off the last point in the required dimension(s). This is not
# correct. This script fixes that.
# Input: grid_file = path to ROMS grid file to edit
def uvp_grid_fix (grid_file):

    # Read x and y coordinates for rho grid
    id = Dataset(grid_file, 'a')
    x_rho = id.variables['x_rho'][:,:]
    y_rho = id.variables['y_rho'][:,:]    

    # Find x and y coordinates for u-grid: halfway between the rho-points
    # in the x direction
    x_u = 0.5*(x_rho[:,:-1] + x_rho[:,1:])
    y_u = 0.5*(y_rho[:,:-1] + y_rho[:,1:])
    # Find longitude and latitude
    lon_u, lat_u = polarstereo_inv(x_u, y_u)

    # Find x and y coordinates for v-grid: halfway between the rho-points
    # in the y direction
    x_v = 0.5*(x_rho[:-1,:] + x_rho[1:,:])
    y_v = 0.5*(y_rho[:-1,:] + y_rho[1:,:])
    # Find longitude and latitude
    lon_v, lat_v = polarstereo_inv(x_v, y_v)

    # Find x and y coordinates for psi-grid: on the corners, i.e. halfway
    # between the rho-points in both directions
    x_psi = 0.5*(x_v[:,:-1] + x_v[:,1:])
    y_psi = 0.5*(y_u[:-1,:] + y_u[1:,:])
    # Find longitude and latitude
    lon_psi, lat_psi = polarstereo_inv(x_psi, y_psi)

    # Save updated variables
    id.variables['x_u'][:,:] = x_u
    id.variables['y_u'][:,:] = y_u
    id.variables['x_v'][:,:] = x_v
    id.variables['y_v'][:,:] = y_v
    id.variables['x_psi'][:,:] = x_psi
    id.variables['y_psi'][:,:] = y_psi
    id.variables['lon_u'][:,:] = lon_u
    id.variables['lat_u'][:,:] = lat_u
    id.variables['lon_v'][:,:] = lon_v
    id.variables['lat_v'][:,:] = lat_v
    id.variables['lon_psi'][:,:] = lon_psi
    id.variables['lat_psi'][:,:] = lat_psi
    id.close()

    
# Given x and y coordinates on a polar stereographic projection, find the
# longitude and latitude by numerically inverting the polar stereographic
# transformation. Adapted from https://au.mathworks.com/matlabcentral/fileexchange/32907-polar-stereographic-coordinate-transformation--map-to-lat-lon-
# Input:
# x, y = coordinates in x-y space
# Output:
# lon, lat = corresponding longitude and latitude in degrees
def polarstereo_inv (x, y):

    a = 6378.137  # Radius of Earth
    e = sqrt(6.694379852e-3)  # Eccentricity
    phi_c = 71  # Projection is tangent to this latitude
    deg2rad = pi/180

    x = -x
    y = -y
    phi_c = phi_c*deg2rad

    t_c = tan(pi/4 - phi_c/2)/((1-e*sin(phi_c))/(1+e*sin(phi_c)))**(e/2)
    m_c = cos(phi_c)/sqrt(1-e**2*(sin(phi_c))**2)
    rho = sqrt(x**2 + y**2)
    t = rho*t_c/(a*m_c)

    chi = pi/2 - 2*arctan(t)
    lat = chi + (e**2/2 + 5*e**4/24 + e**6/12 + 13*e**8/360)*sin(2*chi) + (7*e**4/48 + 29*e**6/240 + 811*e**8/11520)*sin(4*chi) + (7*e**6/120+81*e**8/1120)*sin(6*chi) + (4279*e**8/161280)*sin(8*chi)
    lon = arctan2(x, -y)

    lat = -lat/deg2rad
    lon = -lon/deg2rad
    index = lon < 0
    lon[index] += 360

    return lon, lat


# Command-line interface
if __name__ == '__main__':

    grid_file = raw_input("Path to ROMS grid file: ")
    uvp_grid_fix(grid_file)
