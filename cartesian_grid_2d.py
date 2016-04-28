from numpy import *

# Given ROMS grid variables, calculate Cartesian integrands dx and dy.
# Input:
# lon, lat = 2D arrays containing values for latitude and longitude, with
#            dimension latitude x longitude
# Output:
# dx, dy = 2D arrays (dimension latitude x longitude) containing Cartesian
#          integrands
def cartesian_grid_2d (lon, lat):

    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0

    # Save horizontal dimensions
    num_lat = size(lon, 0)
    num_lon = size(lon, 1)

    # Add or subtract 360 from longitude values which wrap around so that
    # longitude increases monotonically from west to east
    i = tile(arange(1,num_lon+1), (num_lat,1))
    index1 = nonzero((i > 1200)*(lon < 100))
    lon[index1] = lon[index1] + 360
    index2 = nonzero((i < 200)*(lon > 300))
    lon[index2] = lon[index2] - 360

    # Interpolate to get longitude at the edges of each cell
    w_bdry = 0.5*(lon[:,0] + lon[:,-1] - 360)
    middle_lon = 0.5*(lon[:,0:-1] + lon[:,1:])
    e_bdry = 0.5*(lon[:,0] + 360 + lon[:,-1])
    lon_edges = concatenate((w_bdry[:,None], middle_lon, e_bdry[:,None]), axis=1)
    # Subtract to get the change in longitude over each cell
    dlon = abs(lon_edges[:,1:] - lon_edges[:,0:-1])

    # Similarly for latitude; linearly extrapolate for latitude at edges of
    # N/S boundary cells
    middle_lat = 0.5*(lat[0:-1,:] + lat[1:,:])
    s_bdry = 2*lat[0,:] - middle_lat[0,:]
    n_bdry = 2*lat[-1,:] - middle_lat[-1,:]
    lat_edges = concatenate((s_bdry[None,:], middle_lat, n_bdry[None,:]), axis=0)
    dlat = lat_edges[1:,:] - lat_edges[0:-1,:]

    # dx = r*cos(lat)*dlon where lat and dlon are converted to radians
    dx = r*cos(lat*deg2rad)*dlon*deg2rad
    # dy = r*dlat where dlat is converted to radians
    dy = r*dlat*deg2rad

    return dx, dy
