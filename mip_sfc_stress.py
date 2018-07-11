from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.collections import PatchCollection
from rotate_vector_roms import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from patches import *
from unrotate_vector import *

def mip_sfc_stress ():

    # File paths
    roms_grid = '/short/m68/kaa561/metroms_iceshelf/apps/common/grid/circ30S_quarterdegree.nc'
    roms_file = '/short/m68/kaa561/metroms_iceshelf/tmproms/run/intercomparison/stress_firstyear.nc'  # Already averaged over first year
    fesom_mesh_path_lr = '/short/y99/kaa561/FESOM/mesh/meshA/'
    fesom_mesh_path_hr = '/short/y99/kaa561/FESOM/mesh/meshB/'
    fesom_file_lr = '/short/y99/kaa561/FESOM/intercomparison_lowres/output/MK44005.1992.forcing.diag.nc'
    fesom_file_hr = '/short/y99/kaa561/FESOM/intercomparison_highres/output/MK44005.1992.forcing.diag.nc'
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Northern boundaries for plots
    nbdry_acc = -30+90
    nbdry_shelf = -64+90
    # Bounds for colour scale
    colour_bound_acc = 0.25
    colour_bound_shelf = 0.25

    print 'Processing ROMS'
    # Read grid
    id = Dataset(roms_grid, 'r')
    roms_lat = id.variables['lat_rho'][:,:]
    roms_lon = id.variables['lon_rho'][:,:]
    angle = id.variables['angle'][:,:]
    zice = id.variables['zice'][:,:]
    id.close()
    # Read surface stress
    id = Dataset(roms_file, 'r')
    sustr_tmp = id.variables['sustr'][0,:,:]
    svstr_tmp = id.variables['svstr'][0,:,:]
    id.close()
    # Unrotate
    sustr, svstr = rotate_vector_roms(sustr_tmp, svstr_tmp, angle)
    # Get magnitude
    roms_stress = sqrt(sustr**2 + svstr**2)
    # Mask cavities
    roms_stress = ma.masked_where(zice<0, roms_stress)
    # Calculate polar projection
    roms_x = -(roms_lat+90)*cos(roms_lon*deg2rad+pi/2)
    roms_y = (roms_lat+90)*sin(roms_lon*deg2rad+pi/2)

    print 'Processing low-res FESOM'
    # Build mesh and patches
    elements_lr, patches_lr = make_patches(fesom_mesh_path_lr, circumpolar=True, mask_cavities=True)
    # Read rotated and and lon
    f = open(fesom_mesh_path_lr + 'nod2d.out', 'r')
    f.readline()
    rlon_lr = []
    rlat_lr = []
    for line in f:
        tmp = line.split()
        lon_tmp = float(tmp[1])
        if lon_tmp < -180:
            lon_tmp += 360
        elif lon_tmp > 180:
            lon_tmp -= 360
        rlon_lr.append(lon_tmp)
        rlat_lr.append(float(tmp[2]))
    f.close()
    rlon_lr = array(rlon_lr)
    rlat_lr = array(rlat_lr)
    # Read surface stress
    id = Dataset(fesom_file_lr, 'r')
    stress_x_tmp = mean(id.variables['stress_x'][:,:], axis=0)
    stress_y_tmp = mean(id.variables['stress_y'][:,:], axis=0)
    id.close()
    # Unrotate
    stress_x_lr, stress_y_lr = unrotate_vector(rlon_lr, rlat_lr, stress_x_tmp, stress_y_tmp)
    # Get magnitude
    fesom_stress_lr_nodes = sqrt(stress_x_lr**2 + stress_y_lr**2)
    # Average over elements
    fesom_stress_lr = []
    for elm in elements_lr:
        if not elm.cavity:
            fesom_stress_lr.append(mean([fesom_stress_lr_nodes[elm.nodes[0].id], fesom_stress_lr_nodes[elm.nodes[1].id], fesom_stress_lr_nodes[elm.nodes[2].id]]))

    print 'Processing high-res FESOM'
    elements_hr, patches_hr = make_patches(fesom_mesh_path_hr, circumpolar=True, mask_cavities=True)
    f = open(fesom_mesh_path_hr + 'nod2d.out', 'r')
    f.readline()
    rlon_hr = []
    rlat_hr = []
    for line in f:
        tmp = line.split()
        lon_tmp = float(tmp[1])
        if lon_tmp < -180:
            lon_tmp += 360
        elif lon_tmp > 180:
            lon_tmp -= 360
        rlon_hr.append(lon_tmp)
        rlat_hr.append(float(tmp[2]))
    f.close()
    rlon_hr = array(rlon_hr)
    rlat_hr = array(rlat_hr)
    id = Dataset(fesom_file_hr, 'r')
    stress_x_tmp = mean(id.variables['stress_x'][:,:], axis=0)
    stress_y_tmp = mean(id.variables['stress_y'][:,:], axis=0)
    id.close()
    stress_x_hr, stress_y_hr = unrotate_vector(rlon_hr, rlat_hr, stress_x_tmp, stress_y_tmp)
    fesom_stress_hr_nodes = sqrt(stress_x_hr**2 + stress_y_hr**2)
    fesom_stress_hr = []
    for elm in elements_hr:
        if not elm.cavity:
            fesom_stress_hr.append(mean([fesom_stress_hr_nodes[elm.nodes[0].id], fesom_stress_hr_nodes[elm.nodes[1].id], fesom_stress_hr_nodes[elm.nodes[2].id]]))

    print 'Plotting'

    # ACC
    fig = figure(figsize=(19,8))
    fig.patch.set_facecolor('white')
    gs = GridSpec(1,3)
    gs.update(left=0.05, right=0.95, bottom=0.1, top=0.85, wspace=0.05)
    # ROMS
    ax = subplot(gs[0,0], aspect='equal')
    ax.pcolor(roms_x, roms_y, roms_stress, vmin=0, vmax=colour_bound_acc, cmap='jet')
    xlim([-nbdry_acc, nbdry_acc])
    ylim([-nbdry_acc, nbdry_acc])
    title('a) MetROMS', fontsize=28)
    ax.set_xticks([])
    ax.set_yticks([])
    # FESOM (low-res)
    ax = subplot(gs[0,1], aspect='equal')
    img = PatchCollection(patches_lr, cmap='jet')
    img.set_array(array(fesom_stress_lr))
    img.set_clim(vmin=0, vmax=colour_bound_acc)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry_acc, nbdry_acc])
    ylim([-nbdry_acc, nbdry_acc])
    title('b) FESOM (low-res)', fontsize=28)
    ax.set_xticks([])
    ax.set_yticks([])
    # FESOM (high-res)
    ax = subplot(gs[0,2], aspect='equal')
    img = PatchCollection(patches_hr, cmap='jet')
    img.set_array(array(fesom_stress_hr))
    img.set_clim(vmin=0, vmax=colour_bound_acc)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry_acc, nbdry_acc])
    ylim([-nbdry_acc, nbdry_acc])
    title('c) FESOM (high-res)', fontsize=28)
    ax.set_xticks([])
    ax.set_yticks([])
    # Add a horizontal colourbar on the bottom
    cbaxes = fig.add_axes([0.3, 0.05, 0.4, 0.04])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='max', ticks=arange(0, colour_bound_acc+0.05, 0.05))
    cbar.ax.tick_params(labelsize=20)
    # Main title
    suptitle(r'Ocean surface stress (N/m$^2$), 1992 mean', fontsize=34)
    fig.show()
    fig.savefig('sfc_stress_acc.png')

    # Continental shelf
    fig = figure(figsize=(19,8))
    fig.patch.set_facecolor('white')
    gs = GridSpec(1,3)
    gs.update(left=0.05, right=0.95, bottom=0.1, top=0.85, wspace=0.05)
    # ROMS
    ax = subplot(gs[0,0], aspect='equal')
    ax.pcolor(roms_x, roms_y, roms_stress, vmin=0, vmax=colour_bound_shelf, cmap='jet')
    xlim([-nbdry_shelf, nbdry_shelf])
    ylim([-nbdry_shelf, nbdry_shelf])
    title('a) MetROMS', fontsize=28)
    ax.set_xticks([])
    ax.set_yticks([])
    # FESOM (low-res)
    ax = subplot(gs[0,1], aspect='equal')
    img = PatchCollection(patches_lr, cmap='jet')
    img.set_array(array(fesom_stress_lr))
    img.set_clim(vmin=0, vmax=colour_bound_shelf)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry_shelf, nbdry_shelf])
    ylim([-nbdry_shelf, nbdry_shelf])
    title('b) FESOM (low-res)', fontsize=28)
    ax.set_xticks([])
    ax.set_yticks([])
    # FESOM (high-res)
    ax = subplot(gs[0,2], aspect='equal')
    img = PatchCollection(patches_hr, cmap='jet')
    img.set_array(array(fesom_stress_hr))
    img.set_clim(vmin=0, vmax=colour_bound_shelf)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry_shelf, nbdry_shelf])
    ylim([-nbdry_shelf, nbdry_shelf])
    title('c) FESOM (high-res)', fontsize=28)
    ax.set_xticks([])
    ax.set_yticks([])
    # Add a horizontal colourbar on the bottom
    cbaxes = fig.add_axes([0.3, 0.05, 0.4, 0.04])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='max', ticks=arange(0, colour_bound_shelf+0.05, 0.05))
    cbar.ax.tick_params(labelsize=20)
    # Main title
    suptitle(r'Ocean surface stress (N/m$^2$), 1992 mean', fontsize=34)
    fig.show()
    fig.savefig('sfc_stress_shelf.png')


# Command-line interface
if __name__ == "__main__":

    mip_sfc_stress()
