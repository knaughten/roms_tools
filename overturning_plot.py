from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from cartesian_grid_3d import *

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

    # Read v and grid variables
    id = Dataset(file_path, 'r')
    v = id.variables['v'][-1,:,:-15,:-3]  # -1 means last timestep
    h = id.variables['h'][:-15,:-3]
    zice = id.variables['zice'][:-15,:-3]
    zeta = id.variables['zeta'][-1,:-15,:-3]
    lon_v = id.variables['lon_v'][:-15,:-3]
    lat_v = id.variables['lat_v'][:-15,:-3]
    id.close()

    # Interpolate h, zice, and zeta to the v-grid
    h_v = 0.5*(h[0:-1,:] + h[1:,:])
    zice_v = 0.5*(zice[0:-1,:] + zice[1:,:])
    zeta_v = 0.5*(zeta[0:-1,:] + zeta[1:,:])

    # Calculate Cartesian integrands and z-coordinates
    dx, dy, dz, z = cartesian_grid_3d(lon_v, lat_v, h_v, zice_v, theta_s, theta_b, hc, N, zeta_v)

    # Calculate transport in each cell
    transport = v*dx*dz
    # Definite integral over longitude
    transport = sum(transport, axis=2)
    # Indefinite integral over depth; flip before and after so the integral
    # starts at the surface, not the bottom. Also convert to Sv.
    transport = flipud(cumsum(flipud(transport), axis=0))*1e-6

    # Calculate latitude and z coordinates, averaged over longitude,
    # for plotting
    avg_lat = mean(lat_v, axis=1)
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

    file_path = raw_input("Path to ROMS history/averages file: ")
    fig_name = raw_input("File name for figure: ")
    overturning_plot(file_path, fig_name)
