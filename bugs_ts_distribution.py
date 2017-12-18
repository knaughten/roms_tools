from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from cartesian_grid_3d import *
from unesco import *

def bugs_ts_distribution (grid_file, ini_file, upwind_file, akima_file, split_file):

    # Only consider regions with bathymetry deeper than 1500 m, not in ice
    # shelf cavities, and cells deeper than 200 m
    h0 = 1500
    z0 = 200
    # Bounds on temperature and salinity bins
    min_salt = 33.8
    max_salt = 36.6
    min_temp = -2
    max_temp = 21
    # Bounds to actually plot
    min_salt_plot = 34
    max_salt_plot = 34.7
    min_temp_plot = 2
    max_temp_plot = 8
    # Number of temperature and salinity bins
    num_bins = 1000
    # ROMS vertical grid parameters
    theta_s = 7.0
    theta_b = 2.0
    hc = 250
    N = 31

    print 'Setting up bins'
    # Calculate boundaries of temperature bins
    temp_bins = linspace(min_temp, max_temp, num=num_bins)
    # Calculate centres of temperature bins (for plotting)
    temp_centres = 0.5*(temp_bins[:-1] + temp_bins[1:])
    # Repeat for salinity
    salt_bins = linspace(min_salt, max_salt, num=num_bins)
    salt_centres = 0.5*(salt_bins[:-1] + salt_bins[1:])
    # Set up 2D arrays of temperature bins x salinity bins to increment with
    # volume of water masses
    ts_vals_ini = zeros([size(temp_centres), size(salt_centres)])
    ts_vals_upwind = zeros([size(temp_centres), size(salt_centres)])
    ts_vals_akima = zeros([size(temp_centres), size(salt_centres)])
    ts_vals_split = zeros([size(temp_centres), size(salt_centres)])
    # Get 2D versions of the temperature and salinity bins
    salt_2d, temp_2d = meshgrid(salt_centres, temp_centres)
    # Calculate potential density of each combination of temperature and
    # salinity bins
    density = unesco(temp_2d, salt_2d, zeros(shape(temp_centres)))-1000
    # Density contours to plot
    density_lev = arange(26.5, 28+0.25, 0.25)

    print 'Reading grid'
    id = Dataset(grid_file, 'r')
    lon = id.variables['lon_rho'][:,:]
    lat = id.variables['lat_rho'][:,:]
    h = id.variables['h'][:,:]
    zice = id.variables['zice'][:,:]
    id.close()
    num_lat = size(lat, 0)
    num_lon = size(lon, 1)
    # Get integrands on 3D grid
    dx, dy, dz, z = cartesian_grid_3d(lon, lat, h, zice, theta_s, theta_b, hc, N)
    # Get volume integrand
    dV = dx*dy*dz

    print 'Reading data'
    # Read temp and salt from the first history output step, output before
    # timestepping even starts (i.e. initial conditions for January 1992)
    id = Dataset(ini_file, 'r')
    ini_temp = id.variables['temp'][0,:,:,:]
    ini_salt = id.variables['salt'][0,:,:,:]
    id.close()
    # Then read temp and salt averaged over the last January for each simulation
    id = Dataset(upwind_file, 'r')
    upwind_temp = id.variables['temp'][0,:,:,:]
    upwind_salt = id.variables['salt'][0,:,:,:]
    id.close()
    id = Dataset(akima_file, 'r')
    akima_temp = id.variables['temp'][0,:,:,:]
    akima_salt = id.variables['salt'][0,:,:,:]
    id.close()
    id = Dataset(split_file, 'r')
    split_temp = id.variables['temp'][0,:,:,:]
    split_salt = id.variables['salt'][0,:,:,:]
    id.close()

    print 'Binning temperature and salinity'
    # Loop over 2D grid boxes
    for j in range(num_lat):
        for i in range(num_lon):
            # Check for land mask
            if ini_temp[0,j,i] is ma.masked:
                continue
            # Check for ice shelf cavity
            if zice[j,i] < 0:
                continue
            # Check for too-shallow bathymetry
            if h[j,i] < h0:
                continue
            for k in range(N):
                # Check for too-shallow cells
                if abs(z[k,j,i]) < z0:
                    continue
                # First categorise the initial data
                # Figure out which bins this falls into
                temp_index = nonzero(temp_bins > ini_temp[k,j,i])[0][0] - 1
                salt_index = nonzero(salt_bins > ini_salt[k,j,i])[0][0] - 1
                # Increment bins with volume
                ts_vals_ini[temp_index, salt_index] += dV[k,j,i]
                # Upwind simulation
                temp_index = nonzero(temp_bins > upwind_temp[k,j,i])[0][0] - 1
                salt_index = nonzero(salt_bins > upwind_salt[k,j,i])[0][0] - 1
                ts_vals_upwind[temp_index, salt_index] += dV[k,j,i]
                # Akima simulation
                temp_index = nonzero(temp_bins > akima_temp[k,j,i])[0][0] - 1
                salt_index = nonzero(salt_bins > akima_salt[k,j,i])[0][0] - 1
                ts_vals_akima[temp_index, salt_index] += dV[k,j,i]
                # Split simulation
                temp_index = nonzero(temp_bins > split_temp[k,j,i])[0][0] - 1
                salt_index = nonzero(salt_bins > split_salt[k,j,i])[0][0] - 1
                ts_vals_split[temp_index, salt_index] += dV[k,j,i]
    # Mask bins with zero volume
    ts_vals_ini = ma.masked_where(ts_vals_ini==0, ts_vals_ini)
    ts_vals_upwind = ma.masked_where(ts_vals_upwind==0, ts_vals_upwind)
    ts_vals_akima = ma.masked_where(ts_vals_akima==0, ts_vals_akima)
    ts_vals_split = ma.masked_where(ts_vals_split==0, ts_vals_split)

    # Find the volume bounds for plotting
    min_val = log(amin(array([amin(ts_vals_ini), amin(ts_vals_upwind), amin(ts_vals_akima), amin(ts_vals_split)])))
    max_val = log(amax(array([amax(ts_vals_ini), amax(ts_vals_upwind), amax(ts_vals_akima), amax(ts_vals_split)])))

    print 'Plotting'
    fig = figure(figsize=(14,24))
    gs = GridSpec(2,2)
    gs.update(left=0.1, right=0.9, bottom=0.12, top=0.95, wspace=0.05, hspace=0.12)
    # Initial conditions
    ax = subplot(gs[0,0])
    # Plot with log scale
    pcolor(salt_centres, temp_centres, log(ts_vals_ini), vmin=min_val, vmax=max_val, cmap='jet')
    # Add density contours
    cs = contour(salt_centres, temp_centres, density, density_lev, colors=(0.6,0.6,0.6), linestyles='dotted')
    clabel(cs, inline=1, fontsize=14, color=(0.6,0.6,0.6), fmt='%1.1f')
    xlim([min_salt_plot, max_salt_plot])
    ylim([min_temp_plot, max_temp_plot])
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ylabel(r'Temperature ($^{\circ}$C)', fontsize=22)
    title('Initial conditions', fontsize=26)
    # Upwind
    ax = subplot(gs[0,1])
    pcolor(salt_centres, temp_centres, log(ts_vals_upwind), vmin=min_val, vmax=max_val, cmap='jet')
    cs = contour(salt_centres, temp_centres, density, density_lev, colors=(0.6,0.6,0.6), linestyles='dotted')
    clabel(cs, inline=1, fontsize=14, color=(0.6,0.6,0.6), fmt='%1.1f')
    xlim([min_salt_plot, max_salt_plot])
    ylim([min_temp_plot, max_temp_plot])
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    title('Upwind third-order advection', fontsize=26)
    # Akima
    ax = subplot(gs[1,0])
    pcolor(salt_centres, temp_centres, log(ts_vals_akima), vmin=min_val, vmax=max_val, cmap='jet')
    cs = contour(salt_centres, temp_centres, density, density_lev, colors=(0.6,0.6,0.6), linestyles='dotted')
    clabel(cs, inline=1, fontsize=14, color=(0.6,0.6,0.6), fmt='%1.1f')
    xlim([min_salt_plot, max_salt_plot])
    ylim([min_temp_plot, max_temp_plot])
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    xlabel('Salinity (psu)', fontsize=22)
    ylabel(r'Temperature ($^{\circ}$C)', fontsize=22)
    title('Akima advection', fontsize=26)
    # Split
    ax = subplot(gs[1,1])
    img = pcolor(salt_centres, temp_centres, log(ts_vals_split), vmin=min_val, vmax=max_val, cmap='jet')
    cs = contour(salt_centres, temp_centres, density, density_lev, colors=(0.6,0.6,0.6), linestyles='dotted')
    clabel(cs, inline=1, fontsize=14, color=(0.6,0.6,0.6), fmt='%1.1f')
    xlim([min_salt_plot, max_salt_plot])
    ylim([min_temp_plot, max_temp_plot])
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    xlabel('Salinity (psu)', fontsize=22)
    title('RSUP3 advection', fontsize=26)
    # Colorbar at bottom
    cbaxes = fig.add_axes([0.35, 0.04, 0.3, 0.02])
    cbar = colorbar(img, cax=cbaxes, orientation='horizontal')
    cbar.ax.tick_params(labelsize=14)
    text(0.5, 0.01, 'log of volume', fontsize=20, transform=fig.transFigure, ha='center')
    # Main title
    suptitle('Deep water masses after 25 years: AAIW', fontsize=30)
    fig.show()
    fig.savefig('bugs_ts_distribution.png')


# Command-line interface
if __name__ == '__main__':

    grid_file = raw_input("Path to ROMS grid file: ")
    ini_file = raw_input("Path to first ocean_his file, with initial conditions in first time index: ")
    upwind_file = raw_input("Path to January 2016 average file for upwind third-order advection: ")
    akima_file = raw_input("Path to January 2016 average file for Akima advection: ")
    split_file = raw_input("Path to January 2016 average file for RSUP3 advection: ")
    bugs_ts_distribution(grid_file, ini_file, upwind_file, akima_file, split_file)

    
    
    
    

    
