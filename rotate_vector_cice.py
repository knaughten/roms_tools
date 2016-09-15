from netCDF4 import Dataset
from numpy import *

# Given a 2D vector in x-y space on the CICE grid, rotate it to lon-lat space.
# Input:
# u, v = x and y components of the vector on the CICE grid
# angle = angle between the CICE x-axis and east, at each point, in radians
# Output:
# u_lonlat, v_lonlat = components of the vector with respect to lon-lat space
def rotate_vector_cice (u, v, angle):

    u_lonlat = u*cos(-angle) + v*sin(-angle)
    v_lonlat = v*cos(-angle) - u*sin(-angle)

    return u_lonlat, v_lonlat
