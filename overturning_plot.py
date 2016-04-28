from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from calc_z import *

# Calculate the meridional overturning streamfunction at the last timestep
# of the given ROMS history/averages file and make a contour plot in z-space.
# Input:
# file_path = path to ROMS history/averages file
# fig_name = filename for figure 
def overturning_plot (file_path, fig_name):

    # Grid parameters
    theta_s = 0.9
    theta_b = 4.0
    hc = 40
    N = 31
    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0

    # Read v and grid variables
    id = Dataset(file_path, 'r')
    v = id.variables['v'][-1,:,:,:]  # -1 means last timestep
    h = id.variables['h'][:,:]
    zice = id.variables['zice'][:,:]
    zeta = id.variables['zeta'][-1,:,:]
    lon = id.variables['lon_v'][:,:]
    lat = id.variables['lat_v'][:,:]
    id.close()

    num_lat = size(lon,0)
    num_lon = size(lon,1)

    # Add or subtract 360 from longitude values which wrap around
    # so that longitude increases monotonically from west to east
    i = tile(arange(1,num_lon+1), (num_lat,1))
    index1 = nonzero((i > 1200)*(lon < 100))
    lon[index1] = lon[index1] + 360
    index2 = nonzero((i < 200)*(lon > 300))
    lon[index2] = lon[index2] - 360

    # Interpolate to get longitude at the edges of each cell
    w_bdry = 0.5*(lon[:,0] + lon[:,-1] - 360)
    middle_lon = 0.5*(lon[:,0:-1] + lon[:,1:])
    e_bdry = 0.5*(lon[:,0] + 360 + lon[:,-1])
    lon_edges = ma.concatenate((w_bdry[:,None], middle_lon, e_bdry[:,None]), axis=1)
    # Subtract to get the change in longitude over each cell
    dlon = abs(lon_edges[:,1:] - lon_edges[:,0:-1])
    # dx = r*cos(lat)*dlon where lat and dlon are converted to radians
    dx_2d = r*cos(lat*deg2rad)*dlon*deg2rad
    # Copy into a 3D array of dimension depth x latitude x longitude
    dx = tile(dx_2d, (N,1,1))

    # Interpolate h, zice, and zeta to the v-grid
    h_v = 0.5*(h[0:-1,:] + h[1:,:])
    zice_v = 0.5*(zice[0:-1,:] + zice[1:,:])
    zeta_v = 0.5*(zeta[0:-1,:] + zeta[1:,:])

    # Get a 3D grid of z-coordinates; sc_r and Cs_r are unused in this script
    z, sc_r, Cs_r = calc_z(h_v, zice_v, theta_s, theta_b, hc, N, zeta_v)
    # We have z at the midpoint of each cell, now find it on the top and
    # bottom edges of each cell
    z_edges = zeros((size(z,0)+1, size(z,1), size(z,2)))
    z_edges[1:-1,:,:] = 0.5*(z[0:-1,:,:] + z[1:,:,:])
    # At the surface, z = zice; at the bottom, extrapolate
    z_edges[-1,:,:] = zice_v[:,:]
    z_edges[0,:,:] = 2*z[0,:,:] - z_edges[1,:,:]
    # Now find dz
    dz = z_edges[1:,:,:] - z_edges[0:-1,:,:]

    # Calculate transport in each cell
    transport = v*dx*dz
    # Definite integral over longitude
    transport = sum(transport, axis=2)
    # Indefinite integral over depth; flip before and after so the integral
    # starts at the surface, not the bottom. Also convert to Sv.
    transport = flipud(cumsum(flipud(transport), axis=0))*1e-6

    # Calculate latitude and z coordinates, averaged over longitude,
    # for plotting
    avg_lat = mean(lat, axis=1)
    avg_lat = tile(avg_lat, (N,1))
    avg_z = mean(z, axis=2)

    # Centre colour scale on 0
    max_val = amax(abs(transport))
    lev = linspace(-max_val, max_val, num=40)

    # Make the plot
    figure(figsize=(16,8))
    contourf(avg_lat, avg_z, transport, lev, cmap='RdBu_r')
    colorbar()
    xlabel('Latitude')
    ylabel('Depth (m)')
    title('Meridional Overturning Streamfunction (Sv)')

    savefig(fig_name)


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to ROMS history/averages file: ")
    fig_name = raw_input("File name for figure: ")
    overturning_plot(file_path, fig_name)
