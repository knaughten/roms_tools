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

# Input:
# roms_grid = path to ROMS grid file
# roms_seasonal_file = path to seasonal climatology of ROMS 3D temperature and
#                      salinity, precomputed using seasonal_climatology_roms.py
# fesom_mesh_path_lr, fesom_mesh_path_hr = path to FESOM mesh directories for
#                     low-res and high-res meshes
# fesom_seasonal_file_lr, fesom_seasonal_file_hr = paths to seasonal
#                         climatologies of FESOM 3D temperature and salinity
#                         for low-res and high-res respectively, precomputed
#                         using seasonal_climatology.py in the "fesomtools"
#                         repository
def mip_mld (roms_grid, roms_seasonal_file, fesom_mesh_path_lr, fesom_seasonal_file_lr, fesom_mesh_path_hr, fesom_seasonal_file_hr):

    # Path to Sallee's observations
    obs_file = '/short/m68/kaa561/Climatology_MLD003_v2017.nc'
    # Days per month
    days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    # Definition of mixed layer depth: where potential density exceeds
    # surface density by this amount (kg/m^3) as in Sallee et al 2013
    density_anom = 0.03
    # Northern boundary for ACC plot: 30S
    nbdry1 = -30 + 90
    # Northern boundary for continental shelf plot: 64S
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

    print 'Processing low-res FESOM:'
    print 'Building mesh'
    elements_lr, patches_lr = make_patches(fesom_mesh_path_lr, circumpolar, mask_cavities)
    print 'Reading data'
    id = Dataset(fesom_seasonal_file_lr, 'r')
    fesom_temp_nodes_lr = id.variables['temp'][:,:]
    fesom_salt_nodes_lr = id.variables['salt'][:,:]
    id.close()
    print 'Calculating density'
    fesom_density_nodes_lr = unesco(fesom_temp_nodes_lr, fesom_salt_nodes_lr, zeros(shape(fesom_temp_nodes_lr)))
    print 'Calculating mixed layer depth'
    # Set up array for mixed layer depth at each element, at each season
    fesom_mld_lr = zeros([4, len(elements_lr)])
    # Loop over seasons and elements to fill these in
    for season in range(4):
        print '...' + season_names[season]
        mld_season = []
        for elm in elements_lr:
            # Get mixed layer depth at each node
            mld_nodes = []
            for i in range(3):
                node = elm.nodes[i]
                density_sfc = fesom_density_nodes_lr[season,node.id]
                # Save surface depth (only nonzero in ice shelf cavities)
                depth_sfc = node.depth
                temp_depth = node.depth
                curr_node = node.below
                while True:
                    if curr_node is None:
                        # Reached the bottom
                        mld_nodes.append(temp_depth-depth_sfc)
                        break
                    if fesom_density_nodes_lr[season,curr_node.id] >= density_sfc + density_anom:
                        # Reached the critical density anomaly
                        mld_nodes.append(curr_node.depth-depth_sfc)
                        break
                    temp_depth = curr_node.depth
                    curr_node = curr_node.below
            # For this element, save the mean mixed layer depth
            mld_season.append(mean(array(mld_nodes)))
        fesom_mld_lr[season,:] = array(mld_season)

    print 'Processing high-res FESOM:'
    print 'Building mesh'
    elements_hr, patches_hr = make_patches(fesom_mesh_path_hr, circumpolar, mask_cavities)
    print 'Reading data'
    id = Dataset(fesom_seasonal_file_hr, 'r')
    fesom_temp_nodes_hr = id.variables['temp'][:,:]
    fesom_salt_nodes_hr = id.variables['salt'][:,:]
    id.close()
    print 'Calculating density'
    fesom_density_nodes_hr = unesco(fesom_temp_nodes_hr, fesom_salt_nodes_hr, zeros(shape(fesom_temp_nodes_hr)))
    print 'Calculating mixed layer depth'
    # Set up array for mixed layer depth at each element, at each season
    fesom_mld_hr = zeros([4, len(elements_hr)])
    # Loop over seasons and elements to fill these in
    for season in range(4):
        print '...' + season_names[season]
        mld_season = []
        for elm in elements_hr:
            # Get mixed layer depth at each node
            mld_nodes = []
            for i in range(3):
                node = elm.nodes[i]
                density_sfc = fesom_density_nodes_hr[season,node.id]
                # Save surface depth (only nonzero in ice shelf cavities)
                depth_sfc = node.depth
                temp_depth = node.depth
                curr_node = node.below
                while True:
                    if curr_node is None:
                        # Reached the bottom
                        mld_nodes.append(temp_depth-depth_sfc)
                        break
                    if fesom_density_nodes_hr[season,curr_node.id] >= density_sfc + density_anom:
                        # Reached the critical density anomaly
                        mld_nodes.append(curr_node.depth-depth_sfc)
                        break
                    temp_depth = curr_node.depth
                    curr_node = curr_node.below
            # For this element, save the mean mixed layer depth
            mld_season.append(mean(array(mld_nodes)))
        fesom_mld_hr[season,:] = array(mld_season)

    print 'Processing obs'
    # Read grid and monthly climatology
    id = Dataset(obs_file, 'r')
    obs_lon = id.variables['lon'][:]
    obs_lat = id.variables['lat'][:]
    obs_mld_monthly = id.variables['ML_Press'][:,:,:]
    id.close()
    # Polar coordinates for plotting
    obs_lon_2d, obs_lat_2d = meshgrid(obs_lon, obs_lat)
    obs_x = -(obs_lat_2d+90)*cos(obs_lon_2d*deg2rad+pi/2)
    obs_y = (obs_lat_2d+90)*sin(obs_lon_2d*deg2rad+pi/2)
    # Integrate seasonal averages
    obs_mld = zeros([4, size(obs_lat), size(obs_lon)])
    ndays = zeros(4)
    for month in range(12):
        if month+1 in [12, 1, 2]:
            # DJF
            season = 0
        elif month+1 in [3, 4, 5]:
            # MAM
            season = 1
        elif month+1 in [6, 7, 8]:
            # JJA
            season = 2
        elif month+1 in [9, 10, 11]:
            # SON
            season = 3
        obs_mld[season,:,:] += obs_mld_monthly[month,:,:]*days_per_month[month]
        ndays[season] += days_per_month[month]
    # Convert from integrals to averages
    for season in range(4):
        obs_mld[season,:,:] = obs_mld[season,:,:]/ndays[season]
    # Apply land mask
    obs_mld = ma.masked_where(isnan(obs_mld), obs_mld)    

    print 'Plotting'
    # ACC
    fig1 = figure(figsize=(18,9))
    # Summer
    # MetROMS
    ax = fig1.add_subplot(2, 4, 1, aspect='equal')
    pcolor(roms_x, roms_y, roms_mld[0,:,:], vmin=0, vmax=max_bound_summer, cmap='jet')
    text(-67, 0, season_names[0], fontsize=24, ha='right')
    title('MetROMS', fontsize=24)
    xlim([-nbdry1, nbdry1])
    ylim([-nbdry1, nbdry1])
    # Add longitude labels
    for i in range(size(x_ticks)):
        text(x_ticks[i], y_ticks[i], lon_labels[i], ha='center', rotation=lon_rot[i], fontsize=12)
    ax.set_xticks([])
    ax.set_yticks([])
    # FESOM low-res
    ax = fig1.add_subplot(2, 4, 2, aspect='equal')
    img = PatchCollection(patches_lr, cmap='jet')
    img.set_array(fesom_mld_lr[0,:])
    img.set_clim(vmin=0, vmax=max_bound_summer)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry1, nbdry1])
    ylim([-nbdry1, nbdry1])
    ax.set_xticks([])
    ax.set_yticks([])
    title('FESOM (low-res)', fontsize=24)
    # FESOM high-res
    ax = fig1.add_subplot(2, 4, 3, aspect='equal')
    img = PatchCollection(patches_hr, cmap='jet')
    img.set_array(fesom_mld_hr[0,:])
    img.set_clim(vmin=0, vmax=max_bound_summer)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry1, nbdry1])
    ylim([-nbdry1, nbdry1])
    ax.set_xticks([])
    ax.set_yticks([])
    title('FESOM (high-res)', fontsize=24)
    # Obs
    ax = fig1.add_subplot(2, 4, 4, aspect='equal')
    img = pcolor(obs_x, obs_y, obs_mld[0,:,:], vmin=0, vmax=max_bound_summer, cmap='jet')    
    xlim([-nbdry1, nbdry1])
    ylim([-nbdry1, nbdry1])
    ax.set_xticks([])
    ax.set_yticks([])
    title('Observations', fontsize=24)
    # Add a colorbar for summer
    cbaxes = fig1.add_axes([0.93, 0.55, 0.02, 0.3])
    cbar = colorbar(img, cax=cbaxes, extend='max', ticks=arange(0, max_bound_summer+50, 50))
    cbar.ax.tick_params(labelsize=20)
    # Winter
    # MetROMS
    ax = fig1.add_subplot(2, 4, 5, aspect='equal')
    pcolor(roms_x, roms_y, roms_mld[2,:,:], vmin=0, vmax=max_bound_winter, cmap='jet')
    text(-67, 0, season_names[2], fontsize=24, ha='right')
    xlim([-nbdry1, nbdry1])
    ylim([-nbdry1, nbdry1])
    ax.set_xticks([])
    ax.set_yticks([])
    # FESOM low-res
    ax = fig1.add_subplot(2, 4, 6, aspect='equal')
    img = PatchCollection(patches_lr, cmap='jet')
    img.set_array(fesom_mld_lr[2,:])
    img.set_clim(vmin=0, vmax=max_bound_winter)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry1, nbdry1])
    ylim([-nbdry1, nbdry1])
    ax.set_xticks([])
    ax.set_yticks([])
    # FESOM high-res
    ax = fig1.add_subplot(2, 4, 7, aspect='equal')
    img = PatchCollection(patches_hr, cmap='jet')
    img.set_array(fesom_mld_hr[2,:])
    img.set_clim(vmin=0, vmax=max_bound_winter)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry1, nbdry1])
    ylim([-nbdry1, nbdry1])
    ax.set_xticks([])
    ax.set_yticks([])
    # Obs
    ax = fig1.add_subplot(2, 4, 8, aspect='equal')
    img = pcolor(obs_x, obs_y, obs_mld[2,:,:], vmin=0, vmax=max_bound_winter, cmap='jet')
    xlim([-nbdry1, nbdry1])
    ylim([-nbdry1, nbdry1])
    ax.set_xticks([])
    ax.set_yticks([])
    # Add a colorbar for winter
    cbaxes = fig1.add_axes([0.93, 0.15, 0.02, 0.3])
    cbar = colorbar(img, cax=cbaxes, extend='max', ticks=arange(0, max_bound_winter+200, 200))
    cbar.ax.tick_params(labelsize=20)
    # Add the main title
    suptitle('Mixed layer depth (m), 2002-2016 average', fontsize=30)
    # Decrease space between plots
    subplots_adjust(wspace=0.025,hspace=0.025)
    fig1.show()
    fig1.savefig('mld_acc.png')

    # Continental shelf
    fig2 = figure(figsize=(13,9))
    # Summer
    # MetROMS
    ax = fig2.add_subplot(2, 3, 1, aspect='equal')
    pcolor(roms_x, roms_y, roms_mld[0,:,:], vmin=0, vmax=max_bound_summer, cmap='jet')
    text(-28, 0, season_names[0], fontsize=24, ha='right')
    title('MetROMS', fontsize=24)
    xlim([-nbdry2, nbdry2])
    ylim([-nbdry2, nbdry2])
    ax.set_xticks([])
    ax.set_yticks([])
    # FESOM low-res
    ax = fig2.add_subplot(2, 3, 2, aspect='equal')
    img = PatchCollection(patches_lr, cmap='jet')
    img.set_array(fesom_mld_lr[0,:])
    img.set_clim(vmin=0, vmax=max_bound_summer)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry2, nbdry2])
    ylim([-nbdry2, nbdry2])
    ax.set_xticks([])
    ax.set_yticks([])
    title('FESOM (low-res)', fontsize=24)
    # FESOM high-res
    ax = fig2.add_subplot(2, 3, 3, aspect='equal')
    img = PatchCollection(patches_hr, cmap='jet')
    img.set_array(fesom_mld_hr[0,:])
    img.set_clim(vmin=0, vmax=max_bound_summer)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry2, nbdry2])
    ylim([-nbdry2, nbdry2])
    ax.set_xticks([])
    ax.set_yticks([])
    title('FESOM (high-res)', fontsize=24)
    # Add a colorbar for summer
    cbaxes = fig2.add_axes([0.93, 0.55, 0.02, 0.3])
    cbar = colorbar(img, cax=cbaxes, extend='max', ticks=arange(0, max_bound_summer+50, 50))
    cbar.ax.tick_params(labelsize=20)
    # Winter
    # MetROMS
    ax = fig2.add_subplot(2, 3, 4, aspect='equal')
    pcolor(roms_x, roms_y, roms_mld[2,:,:], vmin=0, vmax=max_bound_winter, cmap='jet')
    text(-28, 0, season_names[2], fontsize=24, ha='right')
    xlim([-nbdry2, nbdry2])
    ylim([-nbdry2, nbdry2])
    ax.set_xticks([])
    ax.set_yticks([])
    # FESOM low-res
    ax = fig2.add_subplot(2, 3, 5, aspect='equal')
    img = PatchCollection(patches_lr, cmap='jet')
    img.set_array(fesom_mld_lr[2,:])
    img.set_clim(vmin=0, vmax=max_bound_winter)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry2, nbdry2])
    ylim([-nbdry2, nbdry2])
    ax.set_xticks([])
    ax.set_yticks([])
    # FESOM high-res
    ax = fig2.add_subplot(2, 3, 6, aspect='equal')
    img = PatchCollection(patches_hr, cmap='jet')
    img.set_array(fesom_mld_hr[2,:])
    img.set_clim(vmin=0, vmax=max_bound_winter)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry2, nbdry2])
    ylim([-nbdry2, nbdry2])
    ax.set_xticks([])
    ax.set_yticks([])
    # Add a colorbar for winter
    cbaxes = fig2.add_axes([0.93, 0.15, 0.02, 0.3])
    cbar = colorbar(img, cax=cbaxes, extend='max', ticks=arange(0, max_bound_winter+200, 200))
    cbar.ax.tick_params(labelsize=20)
    # Add the main title
    suptitle('Mixed layer depth (m), 2002-2016 average', fontsize=30)
    # Decrease space between plots
    subplots_adjust(wspace=0.025,hspace=0.025)
    fig2.show()
    fig2.savefig('mld_shelf.png')
                        

# Command-line interface
if __name__ == "__main__":

    roms_grid = raw_input("Path to ROMS grid file: ")
    roms_seasonal_file = raw_input("Path to ROMS seasonal climatology file containing 3D temp and salt: ")
    fesom_mesh_path_lr = raw_input("Path to FESOM low-res mesh directory: ")
    fesom_seasonal_file_lr = raw_input("Path to FESOM low-res seasonal climatology file containing 3D temp and salt: ")
    fesom_mesh_path_hr = raw_input("Path to FESOM high-res mesh directory: ")
    fesom_seasonal_file_hr = raw_input("Path to FESOM high-res seasonal climatology file containing 3D temp and salt: ")
    mip_mld(roms_grid, roms_seasonal_file, fesom_mesh_path_lr, fesom_seasonal_file_lr, fesom_mesh_path_hr, fesom_seasonal_file_hr)

    
                            
                        
                    
                    
            

    
