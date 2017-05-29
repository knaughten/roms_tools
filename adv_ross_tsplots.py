from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from calc_z import *
from interp_lon_roms import *

# Create a 2x2 plot showing zonal slices of temperature and salinity through
# 180E (Ross Sea) at the end of the U3_LIM and C4_LD simulations.
def adv_ross_tsplots ():

    # Paths to simulation directories
    paths = ['/short/m68/kaa561/advection/u3_lim/', '/short/m68/kaa561/advection/c4_l/']
    # Titles for plotting
    labels = [r'a) Temperature ($^{\circ}$C), U3_LIM', r'b) Temperature ($^{\circ}$C), C4_LD - U3_LIM', 'c) Salinity (psu), U3_LIM', 'd) Salinity (psu), C4_LD - U3_LIM']
    # File name: daily average for 31 December
    file_tail = 'ocean_avg_0001.nc'
    var_names = ['temp', 'salt']
    # If 31 December doesn't have its own file, put the time index here
    tstep = 366
    # Longitude to interpolate to
    lon0 = 180
    # Deepest depth to plot
    depth_min = -350 
    # Bounds and ticks for colour scales
    abs_min = [-2, 33.8]
    abs_max = [3, 34.8]    
    abs_ticks = [1, 0.2]
    anom_min = [-1.5, -0.15]
    anom_max = [1.5, 0.15]
    anom_ticks = [0.5, 0.05]
    # Bounds on latitude
    lat_min = -80
    lat_max = -60
    # Grid parameters
    theta_s = 4.0
    theta_b = 0.9
    hc = 40
    N = 31
    Vstretching=2

    # Set up figure
    fig = figure(figsize=(18,12))
    # Loop over variables (temp and salt)
    for var in range(2):                     
        # Read 3D variable
        id = Dataset(paths[0]+file_tail, 'r')
        data_3d = id.variables[var_names[var]][tstep-1,:,:,:]
        # Also read sea surface height
        zeta = id.variables['zeta'][tstep-1,:,:]
        if var==0:
            # For the first simulation, read grid variables
            h = id.variables['h'][:,:]
            zice = id.variables['zice'][:,:]
            lon_2d = id.variables['lon_rho'][:,:]
            lat_2d = id.variables['lat_rho'][:,:]
        id.close()
        # Calculate 3D z-coordinates
        z_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N, zeta, Vstretching)
        # Interpolate temp/salt, z-coordinates, and latitude to 180E
        data0, z, lat = interp_lon_roms(data_3d, z_3d, lat_2d, lon_2d, lon0)
        ax = fig.add_subplot(2, 2, 2*var+1)
        # Shade data (pcolor not contourf so we don't misrepresent the
        # model grid)
        img = pcolor(lat, z, data0, vmin=abs_min[var], vmax=abs_max[var], cmap='jet')
        # Add title
        title(labels[2*var], fontsize=24)
        # Label axes
        if var == 1:
            xlabel('Latitude', fontsize=18)
        ylabel('Depth (m)', fontsize=18)
        xlim([lat_min, lat_max])
        ylim([depth_min, 0])
        # Add colorbars for each variable
        if var == 0:
            cbaxes = fig.add_axes([0.02, 0.575, 0.01, 0.3])
        elif var == 1:
            cbaxes = fig.add_axes([0.02, 0.125, 0.01, 0.3])
        cbar = colorbar(img, ticks=arange(abs_min[var], abs_max[var]+abs_ticks[var], abs_ticks[var]), cax=cbaxes, extend='both')
        cbar.ax.tick_params(labelsize=18)
        # Set ticks the way we want them
        lat_ticks = arange(lat_min+5, lat_max+5, 5)
        ax.set_xticks(lat_ticks)
        lat_labels = []
        for val in lat_ticks:
            lat_labels.append(str(int(round(-val))) + r'$^{\circ}$S')
        ax.set_xticklabels(lat_labels, fontsize=18)
        depth_ticks = range(depth_min+50, 0+100, 100)
        ax.set_yticks(depth_ticks)
        depth_labels = []
        for val in depth_ticks:
            depth_labels.append(str(int(round(-val))))
        ax.set_yticklabels(depth_labels, fontsize=18)

        # Now calculate anomalies for C4_LD simulation
        id = Dataset(paths[1]+file_tail, 'r')
        data_3d = id.variables[var_names[var]][tstep-1,:,:,:]
        # Also read sea surface height
        zeta = id.variables['zeta'][tstep-1,:,:]
        id.close()
        # Calculate 3D z-coordinates
        z_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N, zeta, Vstretching)
        # Interpolate temp/salt, z-coordinates, and latitude to 180E
        data1, z, lat = interp_lon_roms(data_3d, z_3d, lat_2d, lon_2d, lon0)
        # Get anomaly
        data = data1-data0
        ax = fig.add_subplot(2, 2, 2*var+2)
        # Shade data (pcolor not contourf so we don't misrepresent the
        # model grid)
        img = pcolor(lat, z, data, vmin=anom_min[var], vmax=anom_max[var], cmap='RdBu_r')
        # Add title
        title(labels[2*var+1], fontsize=24)
        # Label axes
        if var == 1:
            xlabel('Latitude', fontsize=18)
        xlim([lat_min, lat_max])
        ylim([depth_min, 0])
        # Add colorbars for each variable
        if var == 0:
            cbaxes = fig.add_axes([0.93, 0.575, 0.01, 0.3])
        elif var == 1:
            cbaxes = fig.add_axes([0.93, 0.125, 0.01, 0.3])
        cbar = colorbar(img, ticks=arange(anom_min[var], anom_max[var]+anom_ticks[var], anom_ticks[var]), cax=cbaxes, extend='both')
        cbar.ax.tick_params(labelsize=18)
        # Set ticks the way we want them
        lat_ticks = arange(lat_min+5, lat_max+5, 5)
        ax.set_xticks(lat_ticks)
        lat_labels = []
        for val in lat_ticks:
            lat_labels.append(str(int(round(-val))) + r'$^{\circ}$S')
        ax.set_xticklabels(lat_labels, fontsize=18)
        depth_ticks = range(depth_min+50, 0+100, 100)
        ax.set_yticks(depth_ticks)
        depth_labels = []
        for val in depth_ticks:
            depth_labels.append(str(int(round(-val))))
        ax.set_yticklabels(depth_labels, fontsize=18)

    # Main title
    suptitle(r'180$^{\circ}$E, 31 December (Ross Sea)', fontsize=30)
    #fig.show()
    fig.savefig('adv_ross_tsplots.png')


# Command-line interface
if __name__ == "__main__":

    adv_ross_tsplots()
