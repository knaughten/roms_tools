from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from calc_z import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from patches import *
# This will use the FESOM version of unesco.py for both MetROMS and FESOM,
# luckily it's identical
from unesco import *

# Make a 4x2 plot showing summer (DJF) and winter (JJA) mixed layer depth
# (defined as in Sallee et al 2013: depth at which potential density is 0.03
# kg/m^3 higher than at the surface) comparing MetROMS (top) and FESOM (bottom)
# where each plot is shown zoomed out to include the whole ACC, and zoomed in
# to the continental shelf.
# Input:
# roms_grid = path to ROMS grid file
# roms_seasonal_file = path to seasonal climatology of ROMS 3D temperature and
#                      salinity, precomputed using seasonal_climatology_roms.py
# fesom_mesh_path = path to FESOM mesh directory
# fesom_seasonal_file = path to seasonal climatology of FESOM 3D temperature and
#                       salinity, precomputed using seasonal_climatology.py in
#                       the "fesomtools" repository
def mip_mld (roms_grid, roms_seasonal_file, fesom_mesh_path, fesom_seasonal_file):

    # Definition of mixed layer depth: where potential density exceeds
    # surface density by this amount (kg/m^3) as in Sallee et al 2013
    density_anom = 0.03
    # Northern boundary for ACC plot: 30S
    nbdry1 = -30 + 90
    # Northern boundary for continental shelf plot: 63S
    nbdry2 = -64 + 90
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # FESOM parameters
    circumpolar = True
    mask_cavities = False
    # ROMS parameters
    theta_s = 7.0
    theta_b = 2.0
    hc = 250
    N = 31
    # Season names
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    # Maximum for colour scale in each season
    max_bound_summer = 150
    max_bound_winter = 600
    # Longitude labels for first panel
    lon_ticks = array([-120, -60, 60, 120])
    lat_ticks = array([-28, -25, -25, -28])
    lon_labels = [r'120$^{\circ}$W', r'60$^{\circ}$W', r'60$^{\circ}$E', r'120$^{\circ}$E']
    lon_rot = [-60, 60, -60, 60]

    print 'Processing MetROMS:'
    print 'Reading grid'
    id = Dataset(roms_grid, 'r')
    roms_h = id.variables['h'][:,:]
    roms_zice = id.variables['zice'][:,:]
    roms_lon = id.variables['lon_rho'][:,:]
    roms_lat = id.variables['lat_rho'][:,:]
    id.close()
    # Polar coordinates for plotting
    roms_x = -(roms_lat+90)*cos(roms_lon*deg2rad+pi/2)
    roms_y = (roms_lat+90)*sin(roms_lon*deg2rad+pi/2)
    # Longitude labels
    x_ticks = -(lat_ticks+90)*cos(lon_ticks*deg2rad+pi/2)
    y_ticks = (lat_ticks+90)*sin(lon_ticks*deg2rad+pi/2)
    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    roms_z, sc_r, Cs_r = calc_z(roms_h, roms_zice, theta_s, theta_b, hc, N)
    # Make depth positive
    roms_z = -1*roms_z
    print 'Reading data'
    id = Dataset(roms_seasonal_file, 'r')
    roms_temp = id.variables['temp'][:,:,:,:]
    roms_salt = id.variables['salt'][:,:,:,:]
    id.close()
    print 'Calculating density'
    roms_density = unesco(roms_temp, roms_salt, zeros(shape(roms_temp)))
    print 'Calculating mixed layer depth'
    roms_mld = ma.empty([4, size(roms_lon,0), size(roms_lon,1)])
    # Awful triple loop here, can't find a cleaner way
    for season in range(4):
        print '...' + season_names[season]
        for j in range(size(roms_lon,0)):
            for i in range(size(roms_lon,1)):
                # Get surface density
                density_sfc = roms_density[season,-1,j,i]
                # Get surface depth (only nonzero in ice shelf cavities)
                depth_sfc = roms_z[-1,j,i]
                if density_sfc is ma.masked:
                    # Land
                    roms_mld[season,j,i] = ma.masked
                else:
                    # Loop downward
                    k = size(roms_density,1)-2
                    while True:
                        if k < 0:
                            # Reached the bottom
                            roms_mld[season,j,i] = roms_z[0,j,i]-depth_sfc
                            break
                        if roms_density[season,k,j,i] >= density_sfc + density_anom:
                            # Reached the critical density anomaly
                            roms_mld[season,j,i] = roms_z[k,j,i]-depth_sfc
                            break
                        k -= 1

    print 'Processing FESOM:'
    print 'Building mesh'
    elements, patches = make_patches(fesom_mesh_path, circumpolar, mask_cavities)
    print 'Reading data'
    id = Dataset(fesom_seasonal_file, 'r')
    fesom_temp_nodes = id.variables['temp'][:,:]
    fesom_salt_nodes = id.variables['salt'][:,:]
    id.close()
    print 'Calculating density'
    fesom_density_nodes = unesco(fesom_temp_nodes, fesom_salt_nodes, zeros(shape(fesom_temp_nodes)))
    print 'Calculating mixed layer depth'
    # Set up array for mixed layer depth at each element, at each season
    fesom_mld = zeros([4, len(elements)])
    # Loop over seasons and elements to fill these in
    for season in range(4):
        print '...' + season_names[season]
        mld_season = []
        for elm in elements:
            # Get mixed layer depth at each node
            mld_nodes = []
            for i in range(3):
                node = elm.nodes[i]
                density_sfc = fesom_density_nodes[season,node.id]
                # Save surface depth (only nonzero in ice shelf cavities)
                depth_sfc = node.depth
                temp_depth = node.depth
                curr_node = node.below
                while True:
                    if curr_node is None:
                        # Reached the bottom
                        mld_nodes.append(temp_depth-depth_sfc)
                        break
                    if fesom_density_nodes[season,curr_node.id] >= density_sfc + density_anom:
                        # Reached the critical density anomaly
                        mld_nodes.append(curr_node.depth-depth_sfc)
                        break
                    temp_depth = curr_node.depth
                    curr_node = curr_node.below
            # For this element, save the mean mixed layer depth
            mld_season.append(mean(array(mld_nodes)))
        fesom_mld[season,:] = array(mld_season)

    print 'Plotting'
    fig = figure(figsize=(19,9))
    # MetROMS, summer, ACC
    ax = fig.add_subplot(2, 4, 1, aspect='equal')
    pcolor(roms_x, roms_y, roms_mld[0,:,:], vmin=0, vmax=max_bound_summer, cmap='jet')
    text(-62, 0, 'MetROMS', fontsize=24, ha='right')
    title(season_names[0], fontsize=24)
    xlim([-nbdry1, nbdry1])
    ylim([-nbdry1, nbdry1])
    # Add longitude labels
    for i in range(size(x_ticks)):
        text(x_ticks[i], y_ticks[i], lon_labels[i], ha='center', rotation=lon_rot[i], fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])
    # MetROMS, summer, continental shelf
    ax = fig.add_subplot(2, 4, 2, aspect='equal')
    pcolor(roms_x, roms_y, roms_mld[0,:,:], vmin=0, vmax=max_bound_summer, cmap='jet')
    title(season_names[0], fontsize=24)
    xlim([-nbdry2, nbdry2])
    ylim([-nbdry2, nbdry2])
    ax.set_xticks([])
    ax.set_yticks([])
    # MetROMS, winter, ACC
    ax = fig.add_subplot(2, 4, 3, aspect='equal')
    pcolor(roms_x, roms_y, roms_mld[2,:,:], vmin=0, vmax=max_bound_winter, cmap='jet')
    title(season_names[2], fontsize=24)
    xlim([-nbdry1, nbdry1])
    ylim([-nbdry1, nbdry1])
    ax.set_xticks([])
    ax.set_yticks([])
    # MetROMS, winter, continental shelf
    ax = fig.add_subplot(2, 4, 4, aspect='equal')
    pcolor(roms_x, roms_y, roms_mld[2,:,:], vmin=0, vmax=max_bound_winter, cmap='jet')
    title(season_names[2], fontsize=24)
    xlim([-nbdry2, nbdry2])
    ylim([-nbdry2, nbdry2])
    ax.set_xticks([])
    ax.set_yticks([])
    # FESOM, summer, ACC
    ax = fig.add_subplot(2, 4, 5, aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(fesom_mld[0,:])
    img.set_clim(vmin=0, vmax=max_bound_summer)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry1, nbdry1])
    ylim([-nbdry1, nbdry1])
    ax.set_xticks([])
    ax.set_yticks([])
    text(-62, 0, 'FESOM', fontsize=24, ha='right')
    text(-62, -10, '(high-res)', fontsize=24, ha='right')
    # FESOM, summer, continental shelf
    ax = fig.add_subplot(2, 4, 6, aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(fesom_mld[0,:])
    img.set_clim(vmin=0, vmax=max_bound_summer)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry2, nbdry2])
    ylim([-nbdry2, nbdry2])
    ax.set_xticks([])
    ax.set_yticks([])
    # Add a horizontal colorbar for summer
    cbaxes = fig.add_axes([0.2, 0.04, 0.2, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='max', ticks=arange(0, max_bound_summer+50, 50))
    cbar.ax.tick_params(labelsize=20)
    # FESOM, winter, ACC
    ax = fig.add_subplot(2, 4, 7, aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(fesom_mld[2,:])
    img.set_clim(vmin=0, vmax=max_bound_winter)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry1, nbdry1])
    ylim([-nbdry1, nbdry1])
    ax.set_xticks([])
    ax.set_yticks([])
    # FESOM, winter, continental shelf
    ax = fig.add_subplot(2, 4, 8, aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(fesom_mld[2,:])
    img.set_clim(vmin=0, vmax=max_bound_winter)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry2, nbdry2])
    ylim([-nbdry2, nbdry2])
    ax.set_xticks([])
    ax.set_yticks([])
    # Add a horizontal colorbar for winter
    cbaxes = fig.add_axes([0.6, 0.04, 0.2, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='max', ticks=arange(0, max_bound_winter+200, 200))
    cbar.ax.tick_params(labelsize=20)
    # Add the main title
    suptitle('Mixed layer depth (m), 2002-2016 average', fontsize=30)
    # Decrease space between plots
    subplots_adjust(wspace=0.025,hspace=0.025)
    #fig.show()
    fig.savefig('mld.png')
                        

# Command-line interface
if __name__ == "__main__":

    roms_grid = raw_input("Path to ROMS grid file: ")
    roms_seasonal_file = raw_input("Path to ROMS seasonal climatology file containing 3D temp and salt: ")
    fesom_mesh_path = raw_input("Path to FESOM mesh directory: ")
    fesom_seasonal_file = raw_input("Path to FESOM seasonal climatology file containing 3D temp and salt: ")
    mip_mld(roms_grid, roms_seasonal_file, fesom_mesh_path, fesom_seasonal_file)

    
                            
                        
                    
                    
            

    
