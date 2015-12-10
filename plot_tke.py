from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from calc_z import *

# Calculate and plot the total kinetic energy at each timestep of the given 
# ocean averages file.
def plot_tke (file_path, grid_path):

    # Grid parameters
    theta_s = 0.9
    theta_b = 4.0
    hc = 40
    N = 31
    # Radius of the Earth in m
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0

    # Read grid variables
    id = Dataset(grid_path, 'r')
    h = id.variables['h'][:,:]
    zice = id.variables['zice'][:,:]
    # Interpolate h and zice to u and v grids
    h_u = 0.5*(h[:,0:-1] + h[:,1:])
    h_v = 0.5*(h[0:-1,:] + h[1:,:])
    zice_u = 0.5*(zice[:,0:-1] + zice[:,1:])
    zice_v = 0.5*(zice[0:-1,:] + zice[1:,:])
    lon_u = id.variables['lon_u'][:,:]
    lat_u = id.variables['lat_u'][:,:]
    mask_u = id.variables['mask_u'][:,:]
    lon_v = id.variables['lon_v'][:,:]
    lat_v = id.variables['lat_v'][:,:]
    mask_v = id.variables['mask_v'][:,:]
    id.close()

    # Mask lat and lon at land points
    lon_u = ma.masked_where(mask_u==0, lon_u)
    lat_u = ma.masked_where(mask_u==0, lat_u)
    lon_v = ma.masked_where(mask_v==0, lon_v)
    lat_v = ma.masked_where(mask_v==0, lat_v)
    # Save dimensions
    num_lat_u = size(lon_u, 0)
    num_lon_u = size(lon_u, 1)
    num_lat_v = size(lon_v, 0)
    num_lon_v = size(lon_v, 1)

    # Add or subtract 360 from longitude values which wrap around
    # so that longitude increases monotonically from west to east
    i_u = tile(arange(1, num_lon_u+1), (num_lat_u, 1))
    index1_u = nonzero((i_u > 1200)*(lon_u < 100))
    lon_u[index1_u] = lon_u[index1_u] + 360
    index2_u = nonzero((i_u < 200)*(lon_u > 300))
    lon_u[index2_u] = lon_u[index2_u] - 360
    # Repeat for v grid
    i_v = tile(arange(1, num_lon_v+1), (num_lat_v, 1))
    index1_v = nonzero((i_v > 1200)*(lon_v < 100))
    lon_v[index1_v] = lon_v[index1_v] + 360
    index2_v = nonzero((i_v < 200)*(lon_v > 300))
    lon_v[index2_v] = lon_v[index2_v] - 360

    # Interpolate to get longitude at the edges of each cell
    w_bdry_u = (lon_u[:,0] + lon_u[:,-1] - 360)/2
    middle_lon_u = (lon_u[:,0:-1] + lon_u[:,1:])/2
    e_bdry_u = (lon_u[:,0] + 360 + lon_u[:,-1])/2
    lon_u_edges = ma.concatenate((w_bdry_u[:,None], middle_lon_u, e_bdry_u[:,None]), axis=1)
    # Subtract to get the change in longitude over each cell
    dlon_u = abs(lon_u_edges[:,1:] - lon_u_edges[:,0:-1])
    # Repeat for v grid
    w_bdry_v = (lon_v[:,0] + lon_v[:,-1] - 360)/2
    middle_lon_v = (lon_v[:,0:-1] + lon_v[:,1:])/2
    e_bdry_v = (lon_v[:,0] + 360 + lon_v[:,-1])/2
    lon_v_edges = ma.concatenate((w_bdry_v[:,None], middle_lon_v, e_bdry_v[:,None]), axis=1)
    dlon_v = abs(lon_v_edges[:,1:] - lon_v_edges[:,0:-1])

    # Similarly for latitude
    s_bdry_u = lat_u[0,:]
    middle_lat_u = (lat_u[0:-1,:] + lat_u[1:,:])/2
    n_bdry_u = lat_u[-1,:]*0 - 50
    lat_u_edges = ma.concatenate((s_bdry_u[None,:], middle_lat_u, n_bdry_u[None,:]))
    dlat_u = lat_u_edges[1:,:] - lat_u_edges[0:-1,:]
    # Repeat for v grid
    s_bdry_v = lat_v[0,:]
    middle_lat_v = (lat_v[0:-1,:] + lat_v[1:,:])/2
    n_bdry_v = lat_v[-1,:]*0 - 50
    lat_v_edges = ma.concatenate((s_bdry_v[None,:], middle_lat_v, n_bdry_v[None,:]))
    dlat_v = lat_v_edges[1:,:] - lat_v_edges[0:-1,:]

    # Convert from spherical to Cartesian coordinates
    # dy = r*dlat where dlat is converted to radians
    dy_u_2d = r*dlat_u*pi/180.0
    dy_v_2d = r*dlat_v*pi/180.0
    # Copy into a 3D array, same at each depth level
    dy_u = tile(dy_u_2d, (N,1,1))
    dy_v = tile(dy_v_2d, (N,1,1))
    # dx = r*cos(lat)*dlon where lat and dlon are converted to radians
    dx_u_2d = r*cos(pi*lat_u/180.0)*dlon_u*pi/180.0
    dx_v_2d = r*cos(pi*lat_v/180.0)*dlon_v*pi/180.0
    dx_u = tile(dx_u_2d, (N,1,1))
    dx_v = tile(dx_v_2d, (N,1,1))

    # Get a 3D array of z-coordinates on u and v grids
    # sc_r and Cs_r are unused in this script
    z_u, sc_r, Cs_r = calc_z(h_u, zice_u, lon_u, lat_u, theta_s, theta_b, hc, N)
    z_v, sc_r, Cs_r = calc_z(h_v, zice_v, lon_v, lat_v, theta_s, theta_b, hc, N)
    # We have z at the midpoint of each cell, now find it on the top and
    # bottom edges of each cell
    z_u_edges = zeros((size(z_u,0)+1, size(z_u,1), size(z_u,2)))
    z_u_edges[1:-1,:,:] = 0.5*(z_u[0:-1,:,:] + z_u[1:,:,:])
    # At surface, z = 0; at bottom, set z to be the same as the midpoint of
    # the deepest cell
    z_u_edges[0,:,:] = z_u[0,:,:]
    # Now find dz
    dz_u = z_u_edges[1:,:,:] - z_u_edges[0:-1,:,:]
    # Repeat on the v grid
    z_v_edges = zeros((size(z_v,0)+1, size(z_v,1), size(z_v,2)))
    z_v_edges[1:-1,:,:] = 0.5*(z_v[0:-1,:,:] + z_v[1:,:,:])
    z_v_edges[0,:,:] = z_v[0,:,:]
    dz_v = z_v_edges[1:,:,:] - z_v_edges[0:-1,:,:]
    # Calculate dV for each grid
    dV_u = dx_u*dy_u*dz_u
    dV_v = dx_v*dy_v*dz_v

    # Read time data and convert from seconds to years
    id = Dataset(file_path, 'r')
    time = id.variables['ocean_time'][:]/(365*24*60*60)
    # Read reference density
    rho0 = id.variables['rho0'][:]

    avgke = zeros(size(time))
    for t in range(size(time)):
        # Read u, v, and rho at this timestep
        u = id.variables['u'][t,:,:,:]
        v = id.variables['v'][t,:,:,:]
        rho = id.variables['rho'][t,:,:,:] + rho0
        # Interpolate rho onto u and v grids
        rho_u = 0.5*(rho[:,:,0:-1] + rho[:,:,1:])
        rho_v = 0.5*(rho[:,0:-1,:] + rho[:,1:,:])
        # Integrate 0.5*rho*vel^2 over volume on each grid, add for TKE
        avgke[t] = sum(0.5*rho_u*u**2*dV_u) + sum(0.5*rho_v*v**2*dV_v)
    id.close()

    # Plot results
    clf()
    plot(time, avgke)
    xlabel('Years')
    ylabel('Southern Ocean Total Kinetic Energy (J)')
    show()


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to ocean averages file: ")
    grid_path = raw_input("Path to grid file: ")
    plot_tke(file_path, grid_path)

    

    
    
    
