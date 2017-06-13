from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from cartesian_grid_2d import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from patches import *

# Make a 2x1 plot showing horizontal grid resolution (square root of the area 
# of each rectanuglar grid cell or triangular element) in MetROMS (left) and
# FESOM (right), zoomed into the Antarctic continental shelf.
# Input:
# roms_grid_file = path to ROMS grid file
# fesom_mesh_dir = path to FESOM mesh directory
# save = optional boolean indicating to save the figure, rather than display it
#        on screen
# fig_name = if save=True, filename for figure
def mip_grid_res (roms_grid_file, fesom_mesh_dir, save=False, fig_name=None):

    # Spatial bounds on plot
    lat_max = -63 + 90
    # Bounds on colour scale (km)
    limits = [0, 20]
    # Degrees to radians conversion factor
    deg2rad = pi/180
    # FESOM plotting parameters
    circumpolar = True

    print 'Processing ROMS'
    # Read ROMS grid    
    id = Dataset(roms_grid_file, 'r')
    roms_lon = id.variables['lon_rho'][:,:]
    roms_lat = id.variables['lat_rho'][:,:]
    roms_mask = id.variables['mask_rho'][:,:]
    id.close()
    # Get differentials
    roms_dx, roms_dy = cartesian_grid_2d(roms_lon, roms_lat)
    # Calculate resolution: square root of the area, converted to km
    roms_res = sqrt(roms_dx*roms_dy)*1e-3
    # Apply land mask
    roms_res = ma.masked_where(roms_mask==0, roms_res)
    # Polar coordinates for plotting
    roms_x = -(roms_lat+90)*cos(roms_lon*deg2rad+pi/2)
    roms_y = (roms_lat+90)*sin(roms_lon*deg2rad+pi/2)

    print 'Processing FESOM'
    # Build triangular patches for each element
    elements, patches = make_patches(fesom_mesh_dir, circumpolar)
    # Calculate the resolution at each element
    fesom_res = []
    for elm in elements:
        fesom_res.append(sqrt(elm.area())*1e-3)

    print 'Plotting'
    fig = figure(figsize=(20,9))
    # ROMS
    ax1 = fig.add_subplot(1,2,1, aspect='equal')
    pcolor(roms_x, roms_y, roms_res, vmin=limits[0], vmax=limits[1], cmap='jet')
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    axis('off')
    title('MetROMS', fontsize=24)    
    # FESOM
    ax2 = fig.add_subplot(1,2,2, aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(array(fesom_res))
    img.set_clim(vmin=limits[0], vmax=limits[1])
    img.set_edgecolor('face')
    ax2.add_collection(img)
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    axis('off')
    title('FESOM', fontsize=24)
    cbaxes = fig.add_axes([0.92, 0.2, 0.01, 0.6])
    cbar = colorbar(img, cax=cbaxes, extend='max', ticks=arange(limits[0], limits[1]+5, 5))
    cbar.ax.tick_params(labelsize=20)
    suptitle('Horizontal grid resolution (km)', fontsize=30)
    subplots_adjust(wspace=0.05)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    roms_grid_file = raw_input("Path to ROMS grid file: ")
    fesom_mesh_dir = raw_input("Path to FESOM mesh directory: ")
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    mip_grid_res (roms_grid_file, fesom_mesh_dir, save, fig_name)
    
    
        
    
    
