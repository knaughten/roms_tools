from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from calc_z import *

# Calculate the total ocean heat content at each timestep of the given ocean
# averages file.
def plot_ohc (file_path, grid_path):

    # Grid parameters
    theta_s = 0.9
    theta_b = 4.0
    hc = 40
    N = 31
    # Radius of the Earth in m
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Specific heat of polar seawater (J/K/kg)
    cp = 3974.0
    # Celsius to Kelvin conversion constant
    celsius2kelvin = 273.15

    # Read grid variables
    id = Dataset(grid_path, 'r')
    h = id.variables['h'][:,:]
    zice = id.variables['zice'][:,:]
    lon = id.variables['lon_rho'][:,:]
    lat = id.variables['lat_rho'][:,:]
    mask = id.variables['mask_rho'][:,:]
    id.close()

    # Mask lat and lon at land points
    lon = ma.masked_where(mask==0, lon)
    lat = ma.masked_where(mask==0, lat)
    # Save dimensions
    num_lat = size(lon, 0)
    num_lon = size(lon, 1)

    # Add or subtract 360 from longitude values which wrap around
    # so that longitude increases monotonically from west to east
    i = tile(arange(1, num_lon+1), (num_lat, 1))
    index1 = nonzero((i > 1200)*(lon < 100))
    lon[index1] = lon[index1] + 360
    index2 = nonzero((i < 200)*(lon > 300))
    lon[index2] = lon[index2] - 360

    # Interpolate to get longitude at the edges of each cell
    w_bdry = (lon[:,0] + lon[:,num_lon-1] - 360)/2
    middle_lon = (lon[:,0:num_lon-1] + lon[:,1:num_lon])/2
    e_bdry = (lon[:,0] + 360 + lon[:,num_lon-1])/2
    lon_edges = ma.concatenate((w_bdry[:,None], middle_lon, e_bdry[:,None]), axis=1)
    # Subtract to get the change in longitude over each cell
    dlon = abs(lon_edges[:,1:num_lon+1] - lon_edges[:,0:num_lon])

    # Similarly for latitude
    s_bdry = lat[0,:]
    middle_lat = (lat[0:num_lat-1,:] + lat[1:num_lat,:])/2
    n_bdry = lat[num_lat-1,:]*0 - 50
    lat_edges = ma.concatenate((s_bdry[None,:], middle_lat, n_bdry[None,:]))
    dlat = lat_edges[1:num_lat+1,:] - lat_edges[0:num_lat,:]

    # Convert from spherical to Cartesian coordinates
    # dy = r*dlat where dlat is converted to radians
    dy_2d = r*dlat*pi/180.0
    # Copy into a 3D array, same at each depth level
    dy = tile(dy_2d, (N,1,1))
    # dx = r*cos(lat)*dlon where lat and dlon are converted to radians
    dx_2d = r*cos(pi*lat/180.0)*dlon*pi/180.0
    dx = tile(dx_2d, (N,1,1))

    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    z, sc_r, Cs_r = calc_z(h, zice, lon, lat, theta_s, theta_b, hc, N)
    # We have z at the midpoint of each cell, now find it on the top and
    # bottom edges of each cell
    z_edges = zeros((size(z,0)+1, size(z,1), size(z,2)))
    z_edges[1:-1,:,:] = 0.5*(z[0:-1,:,:] + z[1:,:,:])
    # At surface, z = 0; at bottom, set z to be the same as the midpoint of
    # the deepest cell
    z_edges[0,:,:] = z[0,:,:]
    # Now find dz
    dz = z_edges[1:,:,:] - z_edges[0:-1,:,:]
    dV = dx*dy*dz

    # Read time data and convert from seconds to years
    id = Dataset(file_path, 'r')
    time = id.variables['ocean_time'][:]/(365*24*60*60)
    # Read reference density
    rho0 = id.variables['rho0'][:]

    ohc = zeros(size(time))
    for t in range(size(time)):
        # Read temp and rho at this timestep
        temp = id.variables['temp'][t,:,:,:] + celsius2kelvin
        rho = id.variables['rho'][t,:,:,:] + rho0
        # Integrate temp*rho*cp over volume to get OHC
        ohc[t] = sum(temp*rho*cp*dV)
    id.close()

    # Plot results
    clf()
    plot(time, ohc)
    xlabel('Years')
    ylabel('Ocean Heat Content (J)')
    show()


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to ocean averages file: ")
    grid_path = raw_input("Path to grid file: ")
    plot_ohc(file_path, grid_path)

    

    
    
    
