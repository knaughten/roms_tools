from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from cartesian_grid_2d import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from patches import *

# Make a 3x1 plot showing horizontal grid resolution (square root of the area 
# of each rectanuglar grid cell or triangular element) in MetROMS, FESOM
# low-res, and FESOM high-res, zoomed into the Antarctic continental shelf.
# Input:
# roms_grid_file = path to ROMS grid file
# fesom_mesh_low, fesom_mesh_high = paths to FESOM low-res and high-res mesh
#                                   directories
# save = optional boolean indicating to save the figure, rather than display it
#        on screen
# fig_name = if save=True, filename for figure
def mip_grid_res (roms_grid_file, fesom_mesh_low, fesom_mesh_high, save=False, fig_name=None):

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

    print 'Processing FESOM low-res'
    # Build triangular patches for each element
    elements_low, patches_low = make_patches(fesom_mesh_low, circumpolar)
    # Calculate the resolution at each element
    fesom_res_low = []
    for elm in elements_low:
        fesom_res_low.append(sqrt(elm.area())*1e-3)

    print 'Processing FESOM high-res'
    # Build triangular patches for each element
    elements_high, patches_high = make_patches(fesom_mesh_high, circumpolar)
    # Calculate the resolution at each element
    fesom_res_high = []
    for elm in elements_high:
        fesom_res_high.append(sqrt(elm.area())*1e-3)

    print 'Plotting'
    fig = figure(figsize=(27,9))
    # ROMS
    ax1 = fig.add_subplot(1,3,1, aspect='equal')
    pcolor(roms_x, roms_y, roms_res, vmin=limits[0], vmax=limits[1], cmap='jet')
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    axis('off')
    title('a) MetROMS', fontsize=28)    
    # FESOM low-res
    ax2 = fig.add_subplot(1,3,2, aspect='equal')
    img_low = PatchCollection(patches_low, cmap='jet')
    img_low.set_array(array(fesom_res_low))
    img_low.set_clim(vmin=limits[0], vmax=limits[1])
    img_low.set_edgecolor('face')
    ax2.add_collection(img_low)
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    axis('off')
    title('b) FESOM low-res', fontsize=28)
    # FESOM high-res
    ax3 = fig.add_subplot(1,3,3, aspect='equal')
    img_high = PatchCollection(patches_high, cmap='jet')
    img_high.set_array(array(fesom_res_high))
    img_high.set_clim(vmin=limits[0], vmax=limits[1])
    img_high.set_edgecolor('face')
    ax3.add_collection(img_high)
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    axis('off')
    title('c) FESOM high-res', fontsize=28)
    cbaxes = fig.add_axes([0.92, 0.2, 0.01, 0.6])
    cbar = colorbar(img_high, cax=cbaxes, extend='max', ticks=arange(limits[0], limits[1]+5, 5))
    cbar.ax.tick_params(labelsize=24)
    suptitle('Horizontal grid resolution (km)', fontsize=36)
    subplots_adjust(wspace=0.05)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    roms_grid_file = raw_input("Path to ROMS grid file: ")
    fesom_mesh_low = raw_input("Path to FESOM low-res mesh directory: ")
    fesom_mesh_high = raw_input("Path to FESOM high-res mesh directory: ")
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    mip_grid_res (roms_grid_file, fesom_mesh_low, fesom_mesh_high, save, fig_name)
    
    
        
    
    
