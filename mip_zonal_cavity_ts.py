from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from calc_z import *
from interp_lon_roms import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from fesom_grid import *
from fesom_sidegrid import *

# Create one 3x2 plot for each major ice shelf: zonal slices of temperature
# (top) and salinity (bottom), comparing MetROMS, low-res FESOM, and high-res
# FESOM. Longitudes to slice through, and latitude bounds, are pre-determined.
# Input:
# roms_grid = path to ROMS grid file
# roms_file = path to ROMS time-averaged file containing 3D temp and salt
# fesom_mesh_path_lr, fesom_mesh_path_hr = paths to FESOM mesh directories for
#                     the low-res and high-res simulations
# fesom_file_lr, fesom_file_hr = paths to FESOM time-averaged files containing
#                     3D temp and salt for the low-res and high-res simulations
def mip_zonal_cavity_ts (roms_grid, roms_file, fesom_mesh_path_lr, fesom_file_lr, fesom_mesh_path_hr, fesom_file_hr):

    # Name of each ice shelf
    shelf_names = ['Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    # Beginnings of filenames for figures
    fig_heads = ['larsen_d', 'larsen_c', 'wilkins_georgevi_stange', 'ronne_filchner', 'abbot', 'pig', 'thwaites', 'dotson', 'getz', 'nickerson', 'sulzberger', 'mertz', 'totten_moscowuni', 'shackleton', 'west', 'amery', 'prince_harald', 'baudouin_borchgrevink', 'lazarev', 'nivl', 'fimbul_jelbart_ekstrom', 'brunt_riiser_larsen', 'ross']
    # Longitudes intersecting each ice shelf
    lon0 = [-60, -62, -68, -55, -93, -101, -106, -113, -120, -145, -150, 145, 116, 96, 85, 71, 36, 25, 15, 11, -1, -20, 180]
    # Latitude bounds for each ice shelf
    lat_min = [-73.1, -69.35, -73.1, -82.6, -73.28, -75.4, -75.5, -75, -74.9, -75.9, -77.8, -67.7, -67.17, -66.67, -67.25, -72, -69.7, -71, -70.4, -70.75, -71.83, -75.6, -84.6]
    lat_max = [-72, -66.13, -70, -75.5, -72.3, -74.4, -74.67, -74, -73.5, -75.3, -76.41, -67, -66.5, -64.83, -66.25, -68.5, -68.7, -69.9, -69.33, -69.83, -69.33, -72.9, -77]
    num_shelves = len(shelf_names)
    # ROMS grid parameters
    theta_s = 7.0
    theta_b = 2.0
    hc = 250
    N = 31

    print 'Setting up ROMS'
    # Start with grid
    id = Dataset(roms_grid, 'r')
    h = id.variables['h'][:,:]
    zice = id.variables['zice'][:,:]
    lon_2d = id.variables['lon_rho'][:,:]
    lat_2d = id.variables['lat_rho'][:,:]
    id.close()
    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    z_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N)
    # Read temperature and salinity
    id = Dataset(roms_file, 'r')
    roms_temp_3d = id.variables['temp'][0,:,:,:]
    roms_salt_3d = id.variables['salt'][0,:,:,:]
    id.close()

    print 'Setting up low-res FESOM'
    # Build the regular FESOM grid
    elm2D_lr = fesom_grid(fesom_mesh_path_lr)
    # Read temperature and salinity at every node
    id = Dataset(fesom_file_lr, 'r')
    fesom_temp_nodes_lr = id.variables['temp'][0,:]
    fesom_salt_nodes_lr = id.variables['salt'][0,:]
    id.close()

    print 'Setting up high-res FESOM'
    elm2D_hr = fesom_grid(fesom_mesh_path_hr)
    id = Dataset(fesom_file_hr, 'r')
    fesom_temp_nodes_hr = id.variables['temp'][0,:]
    fesom_salt_nodes_hr = id.variables['salt'][0,:]
    id.close()

    # Loop over ice shelves
    for index in range(num_shelves):
        print 'Processing ' + shelf_names[index]
        # Figure out what to write on the title about longitude
        if lon0[index] < 0:
            lon_string = ' ('+str(-lon0[index])+r'$^{\circ}$W)'
        else:
            lon_string = ' ('+str(lon0[index])+r'$^{\circ}$E)'

        # MetROMS
        # Make sure longitude is between 0 and 360
        roms_lon0 = lon0[index]
        if roms_lon0 < 0:
            roms_lon0 += 360
        # Interpolate to given longitude
        roms_temp, roms_z, roms_lat = interp_lon_roms(roms_temp_3d, z_3d, lat_2d, lon_2d, roms_lon0)
        roms_salt, roms_z, roms_lat = interp_lon_roms(roms_salt_3d, z_3d, lat_2d, lon_2d, roms_lon0)
        # Figure out deepest depth
        flag = (roms_lat >= lat_min[index])*(roms_lat <= lat_max[index])
        depth_min_tmp = amin(roms_z[flag])
        # Round down to nearest 50 metres
        depth_min = floor(depth_min_tmp/50)*50

        # FESOM low-res
        # Build arrays of SideElements making up zonal slices
        selements_temp_lr = fesom_sidegrid(elm2D_lr, fesom_temp_nodes_lr, lon0[index], lat_max[index])
        selements_salt_lr = fesom_sidegrid(elm2D_lr, fesom_salt_nodes_lr, lon0[index], lat_max[index])
        # Build array of quadrilateral patches for the plots, and data values
        # corresponding to each SideElement
        patches_lr = []
        fesom_temp_lr = []
        for selm in selements_temp_lr:
            # Make patch
            coord = transpose(vstack((selm.y, selm.z)))
            patches_lr.append(Polygon(coord, True, linewidth=0.))
            # Save data value
            fesom_temp_lr.append(selm.var)
        fesom_temp_lr = array(fesom_temp_lr)
        # Salinity has same patches but different values
        fesom_salt_lr = []
        for selm in selements_salt_lr:
            fesom_salt_lr.append(selm.var)
        fesom_salt_lr = array(fesom_salt_lr)

        # FESOM high-res
        selements_temp_hr = fesom_sidegrid(elm2D_hr, fesom_temp_nodes_hr, lon0[index], lat_max[index])
        selements_salt_hr = fesom_sidegrid(elm2D_hr, fesom_salt_nodes_hr, lon0[index], lat_max[index])
        patches_hr = []
        fesom_temp_hr = []
        for selm in selements_temp_hr:
            coord = transpose(vstack((selm.y, selm.z)))
            patches_hr.append(Polygon(coord, True, linewidth=0.))
            fesom_temp_hr.append(selm.var)
        fesom_temp_hr = array(fesom_temp_hr)
        fesom_salt_hr = []
        for selm in selements_salt_hr:
            fesom_salt_hr.append(selm.var)
        fesom_salt_hr = array(fesom_salt_hr)

        # Find bounds on each variable
        temp_min = amin(array([amin(roms_temp[flag]), amin(fesom_temp_lr), amin(fesom_temp_hr)]))
        temp_max = amax(array([amax(roms_temp[flag]), amax(fesom_temp_lr), amax(fesom_temp_hr)]))
        salt_min = amin(array([amin(roms_salt[flag]), amin(fesom_salt_lr), amin(fesom_salt_hr)]))
        salt_max = amax(array([amax(roms_salt[flag]), amax(fesom_salt_lr), amax(fesom_salt_hr)]))
        # Plot
        fig = figure(figsize=(24,12))
        # MetROMS temperature
        ax = fig.add_subplot(2, 3, 1)
        pcolor(roms_lat, roms_z, roms_temp, vmin=temp_min, vmax=temp_max, cmap='jet')
        title(r'MetROMS temperature ($^{\circ}$C)', fontsize=20)
        ylabel('Depth (m)', fontsize=16)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        # FESOM low-res temperature
        ax = fig.add_subplot(2, 3, 2)
        img = PatchCollection(patches_lr, cmap='jet')
        img.set_array(fesom_temp_lr)
        img.set_edgecolor('face')
        img.set_clim(vmin=temp_min, vmax=temp_max)
        ax.add_collection(img)
        title(r'FESOM (low-res) temperature ($^{\circ}$C)', fontsize=20)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        # FESOM high-res temperature
        ax = fig.add_subplot(2, 3, 3)
        img = PatchCollection(patches_hr, cmap='jet')
        img.set_array(fesom_temp_hr)
        img.set_edgecolor('face')
        img.set_clim(vmin=temp_min, vmax=temp_max)
        ax.add_collection(img)
        title(r'FESOM (high-res) temperature ($^{\circ}$C)', fontsize=20)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        # Add colorbar for temperature
        cbaxes = fig.add_axes([0.92, 0.575, 0.01, 0.3])
        cbar = colorbar(img, cax=cbaxes)
        cbar.ax.tick_params(labelsize=16)
        # MetROMS salinity
        ax = fig.add_subplot(2, 3, 4)
        pcolor(roms_lat, roms_z, roms_salt, vmin=salt_min, vmax=salt_max, cmap='jet')
        title('MetROMS salinity (psu)', fontsize=20)    
        xlabel('Latitude', fontsize=16)
        ylabel('Depth (m)', fontsize=16)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        # FESOM low-res salinity
        ax = fig.add_subplot(2, 3, 5)
        img = PatchCollection(patches_lr, cmap='jet')
        img.set_array(fesom_salt_lr)
        img.set_edgecolor('face') 
        img.set_clim(vmin=salt_min, vmax=salt_max)
        ax.add_collection(img)
        title(r'FESOM (low-res) salinity (psu)', fontsize=20)
        xlabel('Latitude', fontsize=16)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        # FESOM high-res salinity
        ax = fig.add_subplot(2, 3, 6)
        img = PatchCollection(patches_hr, cmap='jet')
        img.set_array(fesom_salt_hr)
        img.set_edgecolor('face') 
        img.set_clim(vmin=salt_min, vmax=salt_max)
        ax.add_collection(img)
        title(r'FESOM (high-res) salinity (psu)', fontsize=20)
        xlabel('Latitude', fontsize=16)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        # Add colorbar for salinity
        cbaxes = fig.add_axes([0.92, 0.125, 0.01, 0.3])
        cbar = colorbar(img, cax=cbaxes)
        cbar.ax.tick_params(labelsize=16)
        # Main title
        suptitle(shelf_names[index] + lon_string, fontsize=28)
        #fig.show()
        fig.savefig(fig_heads[index] + '_zonal_ts.png')
        

# Command-line interface
if __name__ == "__main__":

    roms_grid = raw_input("Path to ROMS grid file: ")
    roms_file = raw_input("Path to ROMS time-averaged file containing 3D temp and salt: ")
    fesom_mesh_path_lr = raw_input("Path to FESOM low-res mesh directory: ")
    fesom_file_lr = raw_input("Path to FESOM low-res time-averaged file containing 3D temp and salt: ")
    fesom_mesh_path_hr = raw_input("Path to FESOM high-res mesh directory: ")
    fesom_file_hr = raw_input("Path to FESOM high-res time-averaged file containing 3D temp and salt: ")
    mip_zonal_cavity_ts(roms_grid, roms_file, fesom_mesh_path_lr, fesom_file_lr, fesom_mesh_path_hr, fesom_file_hr)
