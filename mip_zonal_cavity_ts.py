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

# Create one 2x2 plot for each major ice shelf: zonal slices of temperature
# (top) and salinity (bottom), comparing MetROMS (left) and FESOM (right).
# Longitudes to slice through, and latitude bounds, are pre-determined.
# Input:
# roms_grid = path to ROMS grid file
# roms_file = path to ROMS time-averaged file containing 3D temp and salt
# fesom_mesh_path = path to FESOM mesh directory
# fesom_file = path to FESOM time-averaged file containing 3D temp and salt
def mip_zonal_cavity_ts (roms_grid, roms_file, fesom_mesh_path, fesom_file):

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

    print 'Setting up FESOM'
    # Build the regular FESOM grid
    elm2D = fesom_grid(fesom_mesh_path)
    # Read temperature and salinity at every node
    id = Dataset(fesom_file, 'r')
    fesom_temp_nodes = id.variables['temp'][0,:]
    fesom_salt_nodes = id.variables['salt'][0,:]
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

        # FESOM
        # Build arrays of SideElements making up zonal slices
        selements_temp = fesom_sidegrid(elm2D, fesom_temp_nodes, lon0[index], lat_max[index])
        selements_salt = fesom_sidegrid(elm2D, fesom_salt_nodes, lon0[index], lat_max[index])
        # Build array of quadrilateral patches for the plots, and data values
        # corresponding to each SideElement
        patches = []
        fesom_temp = []
        for selm in selements_temp:
            # Make patch
            coord = transpose(vstack((selm.y, selm.z)))
            patches.append(Polygon(coord, True, linewidth=0.))
            # Save data value
            fesom_temp.append(selm.var)
        fesom_temp = array(fesom_temp)
        # Salinity has same patches but different values
        fesom_salt = []
        for selm in selements_salt:
            fesom_salt.append(selm.var)
        fesom_salt = array(fesom_salt)

        # Find bounds on each variable
        temp_min = min(amin(roms_temp[flag]), amin(fesom_temp))
        temp_max = max(amax(roms_temp[flag]), amax(fesom_temp))
        salt_min = min(amin(roms_salt[flag]), amin(fesom_salt))
        salt_max = max(amax(roms_salt[flag]), amax(fesom_salt))
        # Plot
        fig = figure(figsize=(18,12))
        # MetROMS temperature
        ax = fig.add_subplot(2, 2, 1)
        pcolor(roms_lat, roms_z, roms_temp, vmin=temp_min, vmax=temp_max, cmap='jet')
        title(r'MetROMS temperature ($^{\circ}$C)', fontsize=20)
        ylabel('Depth (m)', fontsize=16)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        # FESOM temperature
        ax = fig.add_subplot(2, 2, 2)
        img1 = PatchCollection(patches, cmap='jet')
        img1.set_array(fesom_temp)
        img1.set_edgecolor('face')
        ax.add_collection(img1)
        title(r'FESOM temperature ($^{\circ}$C)', fontsize=20)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        # Add colorbar for temperature
        cbaxes1 = fig.add_axes([0.92, 0.575, 0.01, 0.3])
        cbar1 = colorbar(img1, cax=cbaxes1)
        cbar1.ax.tick_params(labelsize=16)
        # MetROMS salinity
        ax = fig.add_subplot(2, 2, 3)
        pcolor(roms_lat, roms_z, roms_salt, vmin=salt_min, vmax=salt_max, cmap='jet')
        title('MetROMS salinity (psu)', fontsize=20)    
        xlabel('Latitude', fontsize=16)
        ylabel('Depth (m)', fontsize=16)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        # FESOM salinity
        ax = fig.add_subplot(2, 2, 4)
        img2 = PatchCollection(patches, cmap='jet')
        img2.set_array(fesom_salt)
        img2.set_edgecolor('face')
        ax.add_collection(img2)
        title(r'FESOM salinity (psu)', fontsize=20)
        xlabel('Latitude', fontsize=16)
        xlim([lat_min[index], lat_max[index]])
        ylim([depth_min, 0])
        # Add colorbar for salinity
        cbaxes2 = fig.add_axes([0.92, 0.125, 0.01, 0.3])
        cbar2 = colorbar(img2, cax=cbaxes2)
        cbar2.ax.tick_params(labelsize=16)
        # Main title
        suptitle(shelf_names[index] + lon_string, fontsize=28)
        #fig.show()
        fig.savefig(fig_heads[index] + '_zonal_ts.png')
        

# Command-line interface
if __name__ == "__main__":

    roms_grid = raw_input("Path to ROMS grid file: ")
    roms_file = raw_input("Path to ROMS time-averaged file containing 3D temp and salt: ")
    fesom_mesh_path = raw_input("Path to FESOM mesh directory: ")
    fesom_file = raw_input("Path to FESOM time-averaged file containing 3D temp and salt: ")
    mip_zonal_cavity_ts(roms_grid, roms_file, fesom_mesh_path, fesom_file)
