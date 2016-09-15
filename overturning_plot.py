from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from cartesian_grid_3d import *
from rotate_vector_roms import *

# Calculate the meridional overturning streamfunction at the last timestep
# of the given ROMS history/averages file and make a contour plot in z-space.
# Input:
# grid_path = path to ROMS grid file
# file_path = path to ROMS history/averages file
# fig_name = filename for figure 
def overturning_plot (grid_path, file_path, fig_name):

    # Grid parameters
    theta_s = 0.9
    theta_b = 4.0
    hc = 40
    N = 31

    # Read angle from the grid file
    grid_id = Dataset(grid_path, 'r')
    angle = grid_id.variables['angle'][:-15,:]
    grid_id.close()
    # Read grid variables
    id = Dataset(file_path, 'r')
    v = id.variables['v'][-1,:,:-15,:-3]  # -1 means last timestep
    h = id.variables['h'][:-15,:-3]
    zice = id.variables['zice'][:-15,:-3]
    zeta = id.variables['zeta'][-1,:-15,:-3]
    lon = id.variables['lon_rho'][:-15,:-3]
    lat = id.variables['lat_rho'][:-15,:-3]
    # Read both velocities in x-y space
    u_xy = id.variables['u'][-1,:,:-15,:]
    v_xy = id.variables['v'][-1,:,:-15,:]    
    id.close()

    # Rotate velocities to lat-lon space
    v = ma.empty([N,v_xy.shape[1]+1,v_xy.shape[2]])
    for k in range(N):
        u_lonlat, v_lonlat = rotate_vector_roms(u_xy[k,:,:], v_xy[k,:,:], angle)
        v[k,:,:] = v_lonlat[:,:]
    # Throw away the periodic boundary overlap
    v = v[:,:,:-3]

    # Calculate Cartesian integrands and z-coordinates
    dx, dy, dz, z = cartesian_grid_3d(lon, lat, h, zice, theta_s, theta_b, hc, N, zeta)

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

    #savefig(fig_name)
    show()


# Command-line interface
if __name__ == "__main__":

    grid_path = raw_input("Path to ROMS grid file: ")
    file_path = raw_input("Path to ROMS history/averages file: ")
    fig_name = raw_input("File name for figure: ")
    overturning_plot(grid_path, file_path, fig_name)
