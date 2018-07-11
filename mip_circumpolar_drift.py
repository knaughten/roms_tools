from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from scipy.interpolate import RegularGridInterpolator
from cartesian_grid_3d import *
from circumpolar_plot import average_btw_depths
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from fesom_grid import *
from unrotate_grid import *

# Needs python/2.7.6

def mip_circumpolar_drift ():

    # File paths
    # ECCO2 initial conditions file for temperature
    ecco2_ini_file = '/short/m68/kaa561/metroms_iceshelf/data/originals/ECCO2/THETA.1440x720x50.199201.nc'
    # ROMS grid file
    roms_grid = '/short/m68/kaa561/metroms_iceshelf/apps/common/grid/circ30S_quarterdegree.nc'
    # ROMS January 2016 mean temp
    roms_end_file = '/short/m68/kaa561/metroms_iceshelf/tmproms/run/intercomparison/temp_salt_jan2016.nc'
    # FESOM mesh paths
    fesom_mesh_path_lr = '/short/y99/kaa561/FESOM/mesh/meshA/'
    fesom_mesh_path_hr = '/short/y99/kaa561/FESOM/mesh/meshB/'
    # FESOM January 2016 mean temp
    fesom_end_file_lr = '/short/y99/kaa561/FESOM/intercomparison_lowres/output/temp_salt_jan2016.nc'
    fesom_end_file_hr = '/short/y99/kaa561/FESOM/intercomparison_highres/output/temp_salt_jan2016.nc'
    # Depth bounds to average between
    shallow_bound = 300
    deep_bound = 1000
    # ROMS grid parameters
    theta_s = 7.0
    theta_b = 2.0
    hc = 250
    N = 31
    deg2rad = pi/180
    # Bound for colour scale
    colour_bound = 3
    # Northern boundary for plot
    nbdry = -50+90

    print 'Processing ECCO2'
    id = Dataset(ecco2_ini_file, 'r')
    ecco_lon_tmp = id.variables['LONGITUDE_T'][:]
    ecco_lat = id.variables['LATITUDE_T'][:]
    ecco_depth = id.variables['DEPTH_T'][:]  # Depth is positive
    ecco_temp_3d_tmp = id.variables['THETA'][0,:,:,:]
    id.close()
    # Wrap periodic boundary
    ecco_lon = zeros(size(ecco_lon_tmp)+2)
    ecco_lon[0] = ecco_lon_tmp[-1]-360
    ecco_lon[1:-1] = ecco_lon_tmp
    ecco_lon[-1] = ecco_lon_tmp[0]+360
    ecco_temp_3d = ma.array(zeros((size(ecco_depth), size(ecco_lat), size(ecco_lon))))
    ecco_temp_3d[:,:,0] = ecco_temp_3d_tmp[:,:,-1]
    ecco_temp_3d[:,:,1:-1] = ecco_temp_3d_tmp
    ecco_temp_3d[:,:,-1] = ecco_temp_3d_tmp[:,:,0]
    # Calculate dz
    ecco_depth_edges = zeros(size(ecco_depth)+1)
    ecco_depth_edges[1:-1] = 0.5*(ecco_depth[:-1] + ecco_depth[1:])
    # Surface is zero
    # Extrapolate for bottom
    ecco_depth_edges[-1] = 2*ecco_depth[-1] - ecco_depth_edges[-2]
    ecco_dz = ecco_depth_edges[1:] - ecco_depth_edges[:-1]
    # Average between bounds
    # Find the first level below shallow_bound
    k_start = nonzero(ecco_depth > shallow_bound)[0][0]
    # Find the first level below deep_bound
    # Don't worry about regions where this hits the seafloor, as they will
    # get masked out in the final plot
    k_end = nonzero(ecco_depth > deep_bound)[0][0]
    # Integrate between
    ecco_temp = sum(ecco_temp_3d[k_start:k_end,:,:]*ecco_dz[k_start:k_end,None,None], axis=0)/sum(ecco_dz[k_start:k_end])
    # Fill land mask with zeros
    index = ecco_temp.mask
    ecco_temp = ecco_temp.data
    ecco_temp[index] = 0.0
    # Prepare interpolation function
    interp_function = RegularGridInterpolator((ecco_lat, ecco_lon), ecco_temp)

    print 'Processing MetROMS'
    # Read grid
    id = Dataset(roms_grid, 'r')
    roms_h = id.variables['h'][:,:]
    roms_zice = id.variables['zice'][:,:]
    roms_mask = id.variables['mask_rho'][:,:]
    roms_lon = id.variables['lon_rho'][:,:]
    roms_lat = id.variables['lat_rho'][:,:]
    num_lon = size(roms_lon,1)
    num_lat = size(roms_lat,0)
    id.close()
    # Interpolate ECCO2 depth-averaged values to the ROMS grid
    roms_temp_ini = interp_function((roms_lat, roms_lon))
    # Apply ROMS land mask
    roms_temp_ini = ma.masked_where(roms_mask==0, roms_temp_ini)
    # Read Jan 2016 values
    id = Dataset(roms_end_file, 'r')
    roms_temp_3d_end = id.variables['temp'][0,:,:,:]
    id.close()
    # Get z and dz
    roms_dx, roms_dy, roms_dz, roms_z = cartesian_grid_3d(roms_lon, roms_lat, roms_h, roms_zice, theta_s, theta_b, hc, N)
    # Vertically average between given depths
    roms_temp_end = average_btw_depths(roms_temp_3d_end, roms_z, roms_dz, [-1*shallow_bound, -1*deep_bound])
    # Mask regions shallower than 1000 m
    roms_temp_ini = ma.masked_where(roms_h < deep_bound, roms_temp_ini)
    roms_temp_end = ma.masked_where(roms_h < deep_bound, roms_temp_end)
    # Mask ice shelf cavities
    roms_temp_ini = ma.masked_where(roms_zice < 0, roms_temp_ini)
    roms_temp_end = ma.masked_where(roms_zice < 0, roms_temp_end)
    # Get difference
    roms_temp_drift = roms_temp_end - roms_temp_ini
    # Convert to spherical coordinates
    roms_x = -(roms_lat+90)*cos(roms_lon*deg2rad+pi/2)
    roms_y = (roms_lat+90)*sin(roms_lon*deg2rad+pi/2)

    print 'Processing low-res FESOM'
    print '...Building mesh'
    elements_lr = fesom_grid(fesom_mesh_path_lr, circumpolar=True)
    # Read rotated lat and lon for each 2D node
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
    # Unrotate grid
    fesom_lon_lr, fesom_lat_lr = unrotate_grid(rlon_lr, rlat_lr)
    # Get longitude in the range (-180, 180) to match ECCO
    index = fesom_lon_lr < 0
    fesom_lon_lr[index] = fesom_lon_lr[index] + 360
    print '...Interpolating ECCO2'
    fesom_temp_nodes_ini_lr = interp_function((fesom_lat_lr, fesom_lon_lr))
    # Read January 2016 temp
    id = Dataset(fesom_end_file_lr, 'r')
    fesom_temp_3d_nodes_end_lr = id.variables['temp'][0,:]
    id.close()
    print '...Looping over elements'
    fesom_temp_ini_lr = []
    fesom_temp_end_lr = []
    patches_lr = []
    for elm in elements_lr:
        # Make sure we're not in an ice shelf cavity, or shallower than deep_bound
        if not elm.cavity:
            if all(array([elm.nodes[0].find_bottom().depth, elm.nodes[1].find_bottom().depth, elm.nodes[2].find_bottom().depth]) > deep_bound):
                # Add a new patch
                coord = transpose(vstack((elm.x, elm.y)))
                patches_lr.append(Polygon(coord, True, linewidth=0.))
                # Average initial temp over element
                fesom_temp_ini_lr.append(mean([fesom_temp_nodes_ini_lr[elm.nodes[0].id], fesom_temp_nodes_ini_lr[elm.nodes[1].id], fesom_temp_nodes_ini_lr[elm.nodes[2].id]]))
                # Vertically integrate final temp for this element
                fesom_temp_end_lr.append(fesom_element_average_btw_depths(elm, shallow_bound, deep_bound, fesom_temp_3d_nodes_end_lr))
    fesom_temp_ini_lr = array(fesom_temp_ini_lr)
    fesom_temp_end_lr = array(fesom_temp_end_lr)
    # Get difference
    fesom_temp_drift_lr = fesom_temp_end_lr - fesom_temp_ini_lr

    print 'Processing high-res FESOM'
    print '...Building mesh'
    elements_hr = fesom_grid(fesom_mesh_path_hr, circumpolar=True)
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
    fesom_lon_hr, fesom_lat_hr = unrotate_grid(rlon_hr, rlat_hr)
    index = fesom_lon_hr < 0
    fesom_lon_hr[index] = fesom_lon_hr[index] + 360
    print '...Interpolating ECCO2'
    fesom_temp_nodes_ini_hr = interp_function((fesom_lat_hr, fesom_lon_hr))
    id = Dataset(fesom_end_file_hr, 'r')
    fesom_temp_3d_nodes_end_hr = id.variables['temp'][0,:]
    id.close()
    print '...Looping over elements'
    fesom_temp_ini_hr = []
    fesom_temp_end_hr = []
    patches_hr = []
    for elm in elements_hr:
        if not elm.cavity:
            if all(array([elm.nodes[0].find_bottom().depth, elm.nodes[1].find_bottom().depth, elm.nodes[2].find_bottom().depth]) > deep_bound):
                coord = transpose(vstack((elm.x, elm.y)))
                patches_hr.append(Polygon(coord, True, linewidth=0.))
                fesom_temp_ini_hr.append(mean([fesom_temp_nodes_ini_hr[elm.nodes[0].id], fesom_temp_nodes_ini_hr[elm.nodes[1].id], fesom_temp_nodes_ini_hr[elm.nodes[2].id]]))
                fesom_temp_end_hr.append(fesom_element_average_btw_depths(elm, shallow_bound, deep_bound, fesom_temp_3d_nodes_end_hr))
    fesom_temp_ini_hr = array(fesom_temp_ini_hr)
    fesom_temp_end_hr = array(fesom_temp_end_hr)
    fesom_temp_drift_hr = fesom_temp_end_hr - fesom_temp_ini_hr

    print 'Plotting'
    fig = figure(figsize=(19,8))
    fig.patch.set_facecolor('white')
    gs = GridSpec(1,3)
    gs.update(left=0.05, right=0.95, bottom=0.1, top=0.85, wspace=0.05)
    # ROMS
    ax = subplot(gs[0,0], aspect='equal')
    ax.pcolor(roms_x, roms_y, roms_temp_drift, vmin=-colour_bound, vmax=colour_bound, cmap='RdBu_r')
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    title('a) MetROMS', fontsize=28)
    ax.set_xticks([])
    ax.set_yticks([])
    # FESOM (low-res)
    ax = subplot(gs[0,1], aspect='equal')
    img = PatchCollection(patches_lr, cmap='RdBu_r')
    img.set_array(fesom_temp_drift_lr)
    img.set_clim(vmin=-colour_bound, vmax=colour_bound)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    title('b) FESOM (low-res)', fontsize=28)
    ax.set_xticks([])
    ax.set_yticks([])
    # FESOM (high-res)
    ax = subplot(gs[0,2], aspect='equal')
    img = PatchCollection(patches_hr, cmap='RdBu_r')
    img.set_array(fesom_temp_drift_hr)
    img.set_clim(vmin=-colour_bound, vmax=colour_bound)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    title('c) FESOM (high-res)', fontsize=28)
    ax.set_xticks([])
    ax.set_yticks([])
    # Add a horizontal colourbar on the bottom
    cbaxes = fig.add_axes([0.3, 0.05, 0.4, 0.04])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, ticks=arange(-colour_bound, colour_bound+1, 1), extend='both')
    cbar.ax.tick_params(labelsize=20)
    # Main title
    suptitle(r'Change in temperature from initial conditions ($^{\circ}$C), '+str(shallow_bound)+'-'+str(deep_bound)+' m average', fontsize=34)
    fig.show()
    fig.savefig('circumpolar_temp_drift.png')                


# Element is assumed to be not in an ice shelf cavity, with bathymetry deeper
# than deep_bound
def fesom_element_average_btw_depths (elm, shallow_bound, deep_bound, data):

    area = elm.area()
    nodes = [elm.nodes[0], elm.nodes[1], elm.nodes[2]]
    # Find the first 3D element which is entirely below shallow_bound
    while True:
        if nodes[0].depth > shallow_bound and nodes[1].depth > shallow_bound and nodes[2].depth > shallow_bound:
            # Save these nodes
            first_nodes = copy(nodes)
            break
        else:
            for i in range(3):
                nodes[i] = nodes[i].below
    # Integrate downward until one of the next nodes is deeper than deep_bound
    integral = 0.0
    volume = 0.0
    while True:
        if nodes[0].below.depth > deep_bound or nodes[1].below.depth > deep_bound or nodes[2].below.depth > deep_bound:
            # Save these nodes
            last_nodes = copy(nodes)
            break
        else:
            # Calculate mean of data at six corners of this triangular prism,
            # and mean depths at three edges
            values_tmp = []
            dz_tmp = []
            for i in range(3):
                values_tmp.append(data[nodes[i].id])
                values_tmp.append(data[nodes[i].below.id])
                dz_tmp.append(abs(nodes[i].depth - nodes[i].below.depth))
                # Get ready for next iteration of while loop
                nodes[i] = nodes[i].below
                # Integrand is mean of data at corners * area of upper face * mean of depths at edges
                integral += mean(array(values_tmp))*area*mean(array(dz_tmp))
                volume += mean(array(dz_tmp))*area
    # Now integrate from shallow_bound to first_nodes by linearly interpolating
    # each node to shallow_bound
    values_tmp = []
    dz_tmp = []
    for i in range(3):
        values_tmp.append(data[first_nodes[i].id])
        id1, id2, coeff1, coeff2 = elm.nodes[i].find_depth(shallow_bound)
        values_tmp.append(coeff1*data[id1] + coeff2*data[id2])
        dz_tmp.append(abs(first_nodes[i].depth - shallow_bound))
    integral += mean(array(values_tmp))*area*mean(array(dz_tmp))
    volume += mean(array(dz_tmp))*area
    # Now integrate from last_nodes to deep_bound by linearly interpolating
    # each node to deep_bound
    values_tmp = []
    dz_tmp = []
    for i in range(3):
        values_tmp.append(data[last_nodes[i].id])
        id1, id2, coeff1, coeff2 = elm.nodes[i].find_depth(deep_bound)
        values_tmp.append(coeff1*data[id1] + coeff2*data[id2])
        dz_tmp.append(abs(deep_bound - last_nodes[i].depth))
    integral += mean(array(values_tmp))*area*mean(array(dz_tmp))
    volume += mean(array(dz_tmp))*area
    # All done; divide integral by volume to get the average
    return integral/volume
