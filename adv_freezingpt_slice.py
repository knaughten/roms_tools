from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from calc_z import *

# Plot the difference from the freezing temperature at a specific timestep
# through the Weddell Sea in the worst-performing advection experiment.
# There is no time-averaging, spatial averaging, or interpolation. This shows
# off the spurious supercooling.
def adv_freezingpt_slice ():

    # Path to ocean history file
    file_path = '/short/m68/kaa561/ROMS-CICE-MCT/tmproms/run/advection/c4_lowdif/ocean_his_0001.nc'
    # Timestep to plot
    tstep = 189
    # i-index to plot (1-based)
    i_val = 1250
    # Deepest depth to plot
    depth_min = -100
    # Bounds on colour scale
    colour_bounds = [-0.3, 0.3]
    # Bounds on latitudes to plot
    lat_min = -78
    lat_max = -72
    save = True
    fig_name = 'adv_freezingpt_slice.png'

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
    tfr = -0.054*salt
    # Calculate difference from freezing point
    deltat = temp - tfr

    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    z_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N, zeta)
    # Select depth and latitude at the given i-index
    z = z_3d[:,:,i_val-1]
    lat = tile(lat_2d[:,i_val-1], (N,1))

    # Determine colour bounds
    if colour_bounds is not None:
        # Specified by user
        scale_min = colour_bounds[0]
        scale_max = colour_bounds[1]
        if scale_min == -scale_max:
            # Centered on zero; use a red-yellow-blue colour scale
            colour_map = 'RdYlBu_r'
        else:
            # Use a rainbow colour scale
            colour_map = 'jet'
    else:
        # Determine automatically
        scale_min = amin(deltat)
        scale_max = amax(deltat)
        colour_map = 'jet'

    # Plot (pcolor not contour to show what each individual cell is doing)
    fig = figure(figsize=(12,6))
    pcolor(lat, z, deltat, vmin=scale_min, vmax=scale_max, cmap=colour_map)
    colorbar()
    title(r'Difference from freezing point ($^{\circ}$C) in Weddell Sea, 7 July')
    xlabel('Latitude')
    ylabel('Depth (m)')
    xlim([lat_min, lat_max])
    ylim([depth_min, 0])

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    adv_freezingpt_slice()
