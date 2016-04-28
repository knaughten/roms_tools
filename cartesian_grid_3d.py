from numpy import *
from cartesian_grid_2d import *
from calc_z import *

# Given ROMS grid variables, calculate Cartesian integrands (dx, dy, dz) as
# well as z-coordinates.
# Input:
# lon, lat, h, zice = 2D arrays containing values for latitude, longitude,
#                     bathymetry, and ice shelf draft. All have dimension
#                     latitude x longitude.
# theta_s, theta_b, hc, N = scalar parameters (check your grid file and roms.in)
# zeta = optional 2D array containing values for sea surface height at the
#        desired timestep
# Output:
# dx, dy, dz, z = 3D arrays (dimension depth x latitude x longitude) containing
#                 Cartesian integrands (dx, dy, dz) and z-coordinates (z).
def cartesian_grid_3d (lon, lat, h, zice, theta_s, theta_b, hc, N, zeta=None):

    # Calculate 2D dx and dy in another script
    dx, dy = cartesian_grid_2d(lon, lat)
    # Copy into 3D arrays, same at each depth level
    dx = tile(dx, (N,1,1))
    dy = tile(dy, (N,1,1))
    # Save horizontal dimensions
    num_lat = size(lon, 0)
    num_lon = size(lon, 1)

    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused
    z, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N, zeta)
    # We have z at the midpoint of each cell, now find it on the top and
    # bottom edges of each cell
    z_edges = zeros((N+1, num_lat, num_lon))
    z_edges[1:-1,:,:] = 0.5*(z[0:-1,:,:] + z[1:,:,:])
    # At surface, z=zice; at bottom, extrapolate
    z_edges[-1,:,:] = zice[:,:]
    z_edges[0,:,:] = 2*z[0,:,:] - z_edges[1,:,:]
    # Now find dz
    dz = z_edges[1:,:,:] - z_edges[0:-1,:,:]    

    return dx, dy, dz, z
