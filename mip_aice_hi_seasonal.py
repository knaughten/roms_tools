from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from patches import *

# Make two 4x2 plots, showing seasonal averages of either sea ice concentration
# or effective thickness (concentration*height), comparing MetROMS (top) and
# FESOM (bottom).
# Input:
# cice_seasonal_file = path to seasonal climatology of CICE variables 'aice' and
#                      'hi', pre-computed using seasonal_climatology_cice.py
# fesom_mesh_path = path to FESOM mesh directory
# fesom_seasonal_file = path to seasonal climatology of FESOM variables 'area'
#                       and 'hice', pre-computed using seasonal_climatology.py
#                       in the "fesomtools" repository
def mip_aice_hi_seasonal (cice_seasonal_file, fesom_mesh_path, fesom_seasonal_file):

    # Northern boundary of plot 50S
    nbdry = -50 + 90
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # FESOM parameters
    circumpolar = True
    mask_cavities = True
    # Season names for plot titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    # Colour bounds
    aice_bounds = [0, 1]
    hi_bounds = [0, 1.5]
    
    print 'Processing MetROMS'
    # Read CICE grid data
    id = Dataset(cice_seasonal_file, 'r')
    cice_lon_tmp = id.variables['TLON'][:-15,:]
    cice_lat_tmp = id.variables['TLAT'][:-15,:]
    cice_aice_tmp = id.variables['aice'][:,:-15,:]
    cice_hi_tmp = id.variables['hi'][:,:-15,:]
    id.close()
    # Wrap the periodic boundary by 1 cell
    cice_lon = ma.empty([size(cice_lon_tmp,0), size(cice_lon_tmp,1)+1])
    cice_lat = ma.empty([size(cice_lat_tmp,0), size(cice_lat_tmp,1)+1])
    cice_lon[:,:-1] = cice_lon_tmp
    cice_lon[:,-1] = cice_lon_tmp[:,0]
    cice_lat[:,:-1] = cice_lat_tmp
    cice_lat[:,-1] = cice_lat_tmp[:,0]
    cice_aice = ma.empty([size(cice_aice_tmp,0), size(cice_aice_tmp,1), size(cice_aice_tmp,2)+1])
    cice_aice[:,:,:-1] = cice_aice_tmp
    cice_aice[:,:,-1] = cice_aice_tmp[:,:,0]
    cice_hi = ma.empty([size(cice_hi_tmp,0), size(cice_hi_tmp,1), size(cice_hi_tmp,2)+1])
    cice_hi[:,:,:-1] = cice_hi_tmp
    cice_hi[:,:,-1] = cice_hi_tmp[:,:,0]
    # Multiply sea ice thickness by concentration to get effective thickness
    # (FESOM is already scaled by concentration in the output)
    cice_hi = cice_aice*cice_hi
    # Polar coordinates for plotting
    cice_x = -(cice_lat+90)*cos(cice_lon*deg2rad+pi/2)
    cice_y = (cice_lat+90)*sin(cice_lon*deg2rad+pi/2)

    print 'Processing FESOM'
    # Build FESOM mesh
    elements, patches = make_patches(fesom_mesh_path, circumpolar, mask_cavities)
    # Read data
    id = Dataset(fesom_seasonal_file, 'r')
    fesom_aice_nodes = id.variables['area'][:,:]
    fesom_hi_nodes = id.variables['hice'][:,:]
    id.close()
    # Count the number of elements not in ice shelf cavities
    num_elm = 0
    for elm in elements:
        if not elm.cavity:
            num_elm += 1
    # Set up arrays for element-averages for each season
    fesom_aice = zeros([4, num_elm])
    fesom_hi = zeros([4, num_elm])
    # Loop over elements to fill these in
    i = 0
    for elm in elements:
        if not elm.cavity:
            # Average over 3 component nodes
            fesom_aice[:,i] = (fesom_aice_nodes[:,elm.nodes[0].id] + fesom_aice_nodes[:,elm.nodes[1].id] + fesom_aice_nodes[:,elm.nodes[2].id])/3
            fesom_hi[:,i] = (fesom_hi_nodes[:,elm.nodes[0].id] + fesom_hi_nodes[:,elm.nodes[1].id] + fesom_hi_nodes[:,elm.nodes[2].id])/3
            i += 1

    print 'Plotting sea ice concentration'
    fig = figure(figsize=(20,9))
    # Loop over seasons
    for season in range(4):
        # MetROMS
        ax = fig.add_subplot(2, 4, season+1, aspect='equal')
        pcolor(cice_x, cice_y, cice_aice[season,:,:], vmin=aice_bounds[0], vmax=aice_bounds[1], cmap='jet')
        if season == 0:
            text(-43, 0, 'MetROMS', fontsize=24, ha='right')
        title(season_names[season], fontsize=24)
        xlim([-nbdry, nbdry])
        ylim([-nbdry, nbdry])
        axis('off')
        # FESOM
        ax = fig.add_subplot(2, 4, season+5, aspect='equal')
        img = PatchCollection(patches, cmap='jet')
        img.set_array(fesom_aice[season,:])
        img.set_clim(vmin=aice_bounds[0], vmax=aice_bounds[1])
        img.set_edgecolor('face')
        ax.add_collection(img)
        xlim([-nbdry, nbdry])
        ylim([-nbdry, nbdry])
        axis('off')
        if season == 0:
            text(-43, 0, 'FESOM', fontsize=24, ha='right')
    # Add a horizontal colorbar at the bottom
    cbaxes = fig.add_axes([0.25, 0.04, 0.5, 0.02])
    cbar = colorbar(img, orientation='horizontal', ticks=arange(aice_bounds[0],aice_bounds[1]+0.25,0.25), cax=cbaxes)
    cbar.ax.tick_params(labelsize=16)
    # Add the main title
    suptitle('Sea ice concentration', fontsize=30)
    # Decrease space between plots
    subplots_adjust(wspace=0.025,hspace=0.025)
    fig.savefig('aice_seasonal.png')

    print 'Plotting sea ice thickness'
    fig = figure(figsize=(20,9))
    for season in range(4):
        # MetROMS
        ax = fig.add_subplot(2, 4, season+1, aspect='equal')
        pcolor(cice_x, cice_y, cice_hi[season,:,:], vmin=hi_bounds[0], vmax=hi_bounds[1], cmap='jet')
        if season == 0:
            text(-43, 0, 'MetROMS', fontsize=24, ha='right')
        title(season_names[season], fontsize=24)
        xlim([-nbdry, nbdry])
        ylim([-nbdry, nbdry])
        axis('off')
        # FESOM
        ax = fig.add_subplot(2, 4, season+5, aspect='equal')
        img = PatchCollection(patches, cmap='jet')
        img.set_array(fesom_hi[season,:])
        img.set_clim(vmin=hi_bounds[0], vmax=hi_bounds[1])
        img.set_edgecolor('face')
        ax.add_collection(img)
        xlim([-nbdry, nbdry])
        ylim([-nbdry, nbdry])
        axis('off')
        if season == 0:
            text(-43, 0, 'FESOM', fontsize=24, ha='right')
    cbaxes = fig.add_axes([0.25, 0.04, 0.5, 0.02])
    cbar = colorbar(img, orientation='horizontal', ticks=arange(hi_bounds[0],hi_bounds[1]+0.5,0.5), cax=cbaxes, extend='max')
    cbar.ax.tick_params(labelsize=16)
    suptitle('Sea ice effective thickness (m)', fontsize=30)
    subplots_adjust(wspace=0.025,hspace=0.025)
    fig.savefig('hi_seasonal.png')


# Command-line interface
if __name__ == "__main__":

    cice_seasonal_file = raw_input("Path to CICE seasonal climatology file containing aice and hi: ")
    fesom_mesh_path = raw_input("Path to FESOM mesh directory: ")
    fesom_seasonal_file = raw_input("Path to FESOM seasonal climatology file containing area and hice: ")
    mip_aice_hi_seasonal(cice_seasonal_file, fesom_mesh_path, fesom_seasonal_file)
        

    
