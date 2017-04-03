from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from calc_z import *

def gadv_freezingpt_slice ():

    # Path to ocean history file
    file_path = '/short/m68/kaa561/advection_bitreproducible/u3/ocean_his_7aug.nc'
    # Timestep to plot
    tstep = 1 #220
    # i-index to plot (1-based)
    i_val = 10
    # Deepest depth to plot
    depth_min = -100
    # Bounds on colour scale
    scale_max = 0.08
    scale_tick = 0.04
    # Bounds on latitudes to plot
    lat_min = -70.3
    lat_max = -66.5
    save = True
    fig_name = 'kn_fig1.png'

    # Grid parameters
    theta_s = 4.0
    theta_b = 0.9
    hc = 40
    N = 31

    # Read temperature, salinity, and grid variables
    id = Dataset(file_path, 'r')
    temp = id.variables['temp'][tstep-1,:,:-15,i_val-1]
    salt = id.variables['salt'][tstep-1,:,:-15,i_val-1]
    h = id.variables['h'][:-15,:]
    zice = id.variables['zice'][:-15,:]
    # Sea surface height is time-dependent
    zeta = id.variables['zeta'][tstep-1,:-15,:]
    lon_2d = id.variables['lon_rho'][:-15,:]
    lat_2d = id.variables['lat_rho'][:-15,:]
    id.close()

    # Calculate freezing point as seen by supercooling code
    tfr = salt/(-18.48 + salt*18.48/1000.0)
    # Calculate difference from freezing point
    deltat = temp - tfr

    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    z_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N, zeta)
    # Select depth and latitude at the given i-index
    z = z_3d[:,:,i_val-1]
    lat = tile(lat_2d[:,i_val-1], (N,1))

    # Plot (pcolor not contour to show what each individual cell is doing)
    fig, ax = subplots(figsize=(12,6))
    pcolor(lat, z, deltat, vmin=-scale_max, vmax=scale_max, cmap='RdBu_r')
    cbar = colorbar(extend='both', ticks=arange(-scale_max, scale_max+scale_tick, scale_tick))
    cbar.ax.tick_params(labelsize=14)
    title(r'Difference from freezing point ($^{\circ}$C) near 8$^{\circ}$E: UP3', fontsize=18)
    xlabel('Latitude', fontsize=16)
    ylabel('Depth (m)', fontsize=16)
    xlim([lat_min, lat_max])
    ylim([depth_min, 0])

    lat_ticks = arange(lat_min+0.5, lat_max+0.5, 1)
    xticks(lat_ticks)
    lat_labels = []
    for val in lat_ticks:
        lat_labels.append(str(int(round(-val))) + r'$^{\circ}$S')
    ax.set_xticklabels(lat_labels, fontsize=14)
    depth_ticks = arange(depth_min, 0+25, 25)
    yticks(depth_ticks)
    depth_labels = []
    for val in depth_ticks:
        depth_labels.append(str(int(round(-val))))
    ax.set_yticklabels(depth_labels, fontsize=14)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    gadv_freezingpt_slice()
