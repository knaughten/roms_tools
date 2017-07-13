from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.colors import *
from cartesian_grid_3d import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from fesom_grid import *
from unesco import *

# Make a 2x1 plot of T/S distributions south of 65S, colour-coded based on
# depth, in MetROMS (left) and FESOM (right). Include the surface freezing
# point and density contours.
# Input:
# roms_grid = path to ROMS grid file
# roms_file = path to time-averaged ROMS file containing temperature and
#             salinity (I used 2002-2016 average)
# fesom_mesh_path = path to FESOM mesh directory (I used high-res)
# fesom_file = path to time-averaged FESOM file containing temperature and
#              salinity, over the same period as roms_file
def mip_ts_distribution (roms_grid, roms_file, fesom_mesh_path, fesom_file):

    # Northern boundary of water masses to consider
    nbdry = -65
    # Number of temperature and salinity bins
    num_bins = 1000
    # Bounds on temperature and salinity bins (pre-computed, change if needed)
    min_salt = 32.3
    max_salt = 35.1
    min_temp = -3.1
    max_temp = 3.8
    # Bounds to actually plot
    min_salt_plot = 33.25
    max_salt_plot = 35.0
    min_temp_plot = -3
    max_temp_plot = 3.8
    # FESOM grid generation parameters
    circumpolar = False
    cross_180 = False
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
    # Set up 2D arrays of temperature bins x salinity bins to hold average
    # depth of water masses, weighted by volume
    ts_vals_roms = zeros([size(temp_centres), size(salt_centres)])
    ts_vals_fesom = zeros([size(temp_centres), size(salt_centres)])
    # Also arrays to integrate volume
    volume_roms = zeros([size(temp_centres), size(salt_centres)])
    volume_fesom = zeros([size(temp_centres), size(salt_centres)])
    # Calculate surface freezing point as a function of salinity as seen by
    # each sea ice model
    freezing_pt_roms = salt_centres/(-18.48 + 18.48/1e3*salt_centres)
    freezing_pt_fesom = -0.0575*salt_centres + 1.7105e-3*sqrt(salt_centres**3) - 2.155e-4*salt_centres**2
    # Get 2D versions of the temperature and salinity bins
    salt_2d, temp_2d = meshgrid(salt_centres, temp_centres)
    # Calculate potential density of each combination of temperature and
    # salinity bins
    density = unesco(temp_2d, salt_2d, zeros(shape(temp_centres)))-1000
    # Density contours to plot
    density_lev = arange(26.6, 28.4, 0.2)

    print 'Processing ROMS'
    # Read ROMS grid variables we need
    id = Dataset(roms_grid, 'r')
    roms_lon = id.variables['lon_rho'][:,:]
    roms_lat = id.variables['lat_rho'][:,:]
    roms_h = id.variables['h'][:,:]
    roms_zice = id.variables['zice'][:,:]
    id.close()
    num_lat = size(roms_lat, 0)
    num_lon = size(roms_lon, 1)
    # Get integrands on 3D grid
    roms_dx, roms_dy, roms_dz, roms_z = cartesian_grid_3d(roms_lon, roms_lat, roms_h, roms_zice, theta_s, theta_b, hc, N)
    # Get volume integrand
    dV = roms_dx*roms_dy*roms_dz
    # Read ROMS output
    id = Dataset(roms_file, 'r')
    roms_temp = id.variables['temp'][0,:,:,:]
    roms_salt = id.variables['salt'][0,:,:,:]
    id.close()
    # Loop over 2D grid boxes
    for j in range(num_lat):
        for i in range(num_lon):
            # Check for land mask
            if roms_temp[0,j,i] is ma.masked:
                continue
            # Check if we're in the region of interest
            if roms_lat[j,i] < nbdry:
                # Loop downward
                for k in range(N):
                    # Figure out which bins this falls into
                    temp_index = nonzero(temp_bins > roms_temp[k,j,i])[0][0] - 1
                    salt_index = nonzero(salt_bins > roms_salt[k,j,i])[0][0] - 1
                    # Integrate depth*dV in this bin
                    ts_vals_roms[temp_index, salt_index] += -roms_z[k,j,i]*dV[k,j,i]
                    volume_roms[temp_index, salt_index] += dV[k,j,i]
    # Mask bins with zero volume
    ts_vals_roms = ma.masked_where(volume_roms==0, ts_vals_roms)
    volume_roms = ma.masked_where(volume_roms==0, volume_roms)
    # Convert depths from integrals to volume-averages
    ts_vals_roms /= volume_roms

    print 'Processing FESOM'
    # Make FESOM grid elements
    elements = fesom_grid(fesom_mesh_path, circumpolar, cross_180)
    # Read temperature and salinity at each 3D node
    id = Dataset(fesom_file, 'r')
    fesom_temp = id.variables['temp'][0,:]
    fesom_salt = id.variables['salt'][0,:]
    id.close()
    # Loop over elements
    for elm in elements:
        # See if we're in the region of interest
        if all(elm.lat < nbdry):
            # Get area of 2D triangle
            area = elm.area()
            nodes = [elm.nodes[0], elm.nodes[1], elm.nodes[2]]
            # Loop downward
            while True:
                if nodes[0].below is None or nodes[1].below is None or nodes[2].below is None:
                    # We've reached the bottom
                    break
                # Calculate average temperature, salinity, depth, and layer
                # thickness over this 3D triangular prism
                temp_vals = []
                salt_vals = []
                depth_vals = []
                dz = []
                for i in range(3):
                    # Average temperature over 6 nodes
                    temp_vals.append(fesom_temp[nodes[i].id])
                    temp_vals.append(fesom_temp[nodes[i].below.id])
                    # Average salinity over 6 nodes
                    salt_vals.append(fesom_salt[nodes[i].id])
                    salt_vals.append(fesom_salt[nodes[i].below.id])
                    # Average depth over 6 nodes
                    depth_vals.append(nodes[i].depth)
                    depth_vals.append(nodes[i].below.depth)
                    # Average dz over 3 vertical edges
                    dz.append(abs(nodes[i].depth - nodes[i].below.depth))
                    # Get ready for next repetition of loop
                    nodes[i] = nodes[i].below
                temp_elm = mean(array(temp_vals))
                salt_elm = mean(array(salt_vals))
                depth_elm = mean(array(depth_vals))
                # Calculate volume of 3D triangular prism
                volume = area*mean(array(dz))
                # Figure out which bins this falls into
                temp_index = nonzero(temp_bins > temp_elm)[0][0] - 1
                salt_index = nonzero(salt_bins > salt_elm)[0][0] - 1
                # Integrate depth*volume in this bin
                ts_vals_fesom[temp_index, salt_index] += depth_elm*volume
                volume_fesom[temp_index, salt_index] += volume
    # Mask bins with zero volume
    ts_vals_fesom = ma.masked_where(volume_fesom==0, ts_vals_fesom)
    volume_fesom = ma.masked_where(volume_fesom==0, volume_fesom)
    # Convert depths from integrals to volume-averages
    ts_vals_fesom /= volume_fesom

    # Find the maximum depth for plotting
    max_depth = max(amax(ts_vals_roms), amax(ts_vals_fesom))
    # Make a nonlinear scale
    bounds = linspace(0, max_depth**(1.0/2.5), num=100)**2.5
    norm = BoundaryNorm(boundaries=bounds, ncolors=256)
    # Set labels for density contours
    manual_locations = [(33.4, 3.0), (33.65, 3.0), (33.9, 3.0), (34.2, 3.0), (34.45, 3.5), (34.65, 3.25), (34.9, 3.0), (34.8, 0)]

    print "Plotting"
    fig = figure(figsize=(20,12))
    # ROMS
    ax = fig.add_subplot(1, 2, 1)
    pcolor(salt_centres, temp_centres, ts_vals_roms, norm=norm, vmin=0, vmax=max_depth, cmap='jet')
    # Add surface freezing point line
    plot(salt_centres, freezing_pt_roms, color='black', linestyle='dashed')
    # Add density contours
    cs = contour(salt_centres, temp_centres, density, density_lev, colors=(0.6,0.6,0.6), linestyles='dotted')
    clabel(cs, inline=1, fontsize=14, color=(0.6,0.6,0.6), fmt='%1.1f', manual=manual_locations)
    xlim([min_salt_plot, max_salt_plot])
    ylim([min_temp_plot, max_temp_plot])
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    xlabel('Salinity (psu)', fontsize=22)
    ylabel(r'Temperature ($^{\circ}$C)', fontsize=22)
    title('MetROMS', fontsize=28)
    # FESOM
    ax = fig.add_subplot(1, 2, 2)
    img = pcolor(salt_centres, temp_centres, ts_vals_fesom, norm=norm, vmin=0, vmax=max_depth, cmap='jet')
    plot(salt_centres, freezing_pt_fesom, color='black', linestyle='dashed')
    cs = contour(salt_centres, temp_centres, density, density_lev, colors=(0.6,0.6,0.6), linestyles='dotted')
    clabel(cs, inline=1, fontsize=14, color=(0.6,0.6,0.6), fmt='%1.1f', manual=manual_locations)
    xlim([min_salt_plot, max_salt_plot])
    ylim([min_temp_plot, max_temp_plot])
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    xlabel('Salinity (psu)', fontsize=22)
    title('FESOM (high-res)', fontsize=28)
    # Add a colourbar on the right
    cbaxes = fig.add_axes([0.93, 0.2, 0.02, 0.6])
    cbar = colorbar(img, cax=cbaxes, ticks=[0,50,100,200,500,1000,2000,4000])
    cbar.ax.tick_params(labelsize=18)
    # Add the main title
    suptitle('Water masses south of 65$^{\circ}$S: depth (m), 2002-2016 average', fontsize=30)
    subplots_adjust(wspace=0.1)
    #fig.show()
    fig.savefig('ts_distribution_orig.png')


# Command-line interface
if __name__ == "__main__":

    roms_grid = raw_input("Path to ROMS grid file: ")
    roms_file = raw_input("Path to time-averaged ROMS file containing temperature and salinity: ")
    fesom_mesh_path = raw_input("Path to FESOM mesh directory: ")
    fesom_file = raw_input("Path to time-averaged FESOM file containing temperature and salinity: ")
    mip_ts_distribution(roms_grid, roms_file, fesom_mesh_path, fesom_file)
    
            
      
                                         

    
