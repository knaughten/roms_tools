from netCDF4 import Dataset
from numpy import *
from numpy.ma import MaskedArray
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.colors import LinearSegmentedColormap
from rotate_vector_roms import *
from cartesian_grid_3d import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from fesom_grid import *
from patches import *
from unrotate_vector import *
from unrotate_grid import *

# For each of the 8 ice shelf regions analysed in this intercomparison paper,
# make one 3x1 figure for 5 different variables (ice shelf draft, ice shelf
# melt rate, bottom water temperature, bottom water salinity, vertically
# averaged velocity) where each variable is averaged over the years 2002-2016
# and zoomed into the region of interest. The 3 parts of each figure are
# MetROMS, FESOM low-res, and FESOM high-res. The velocity figures show
# absolute value (speed) in colours, and direction with vector arrows.
def mip_regions_1var ():

    # Path to ROMS grid file
    roms_grid = '/short/m68/kaa561/metroms_iceshelf/apps/common/grid/circ30S_quarterdegree.nc'
    # Path to ROMS time-averaged file
    roms_file = '/short/m68/kaa561/metroms_iceshelf/tmproms/run/intercomparison/2002_2016_avg.nc'
    # Path to FESOM mesh directories
    fesom_mesh_path_lr = '/short/y99/kaa561/FESOM/mesh/low_res/'
    fesom_mesh_path_hr = '/short/y99/kaa561/FESOM/mesh/high_res/'
    # Path to FESOM time-averaged ocean files (temp, salt, u, v)
    fesom_file_lr_o = '/short/y99/kaa561/FESOM/intercomparison_lowres/output/oce_2002_2016_avg.nc'
    fesom_file_hr_o = '/short/y99/kaa561/FESOM/intercomparison_highres/output/oce_2002_2016_avg.nc'
    # Path to FESOM time-averaged ice shelf files (wnet)
    fesom_file_lr_i = '/short/y99/kaa561/FESOM/intercomparison_lowres/output/wnet_2002_2016_avg.nc'
    fesom_file_hr_i = '/short/y99/kaa561/FESOM/intercomparison_highres/output/wnet_2002_2016_avg.nc'

    # Name of each region
    region_names = ['Filchner-Ronne Ice Shelf', 'Eastern Weddell Region', 'Amery Ice Shelf', 'Australian Sector', 'Ross Sea', 'Amundsen Sea', 'Bellingshausen Sea', 'Larsen Ice Shelves']
    num_regions = len(region_names)
    # Beginning of filenames for figures
    fig_heads = ['filchner_ronne', 'eweddell', 'amery', 'australian', 'ross', 'amundsen', 'bellingshausen', 'larsen']
    # Bounds for each region (using polar coordinate transformation as below)
    x_min = [-14, -8, 15.25, 12, -9.5, -15.5, -20.25, -22.5]
    x_max = [-4.5, 13, 20.5, 25.5, 4, -10.5, -15.5, -14.5]
    y_min = [1, 12, 4.75, -20, -13, -11.25, -4.5, 8.3]
    y_max = [10, 21, 8, 4, -4.75, -2.25, 7.6, 13]
    # Size of each plot in the y direction
    ysize = [8, 6, 7, 9, 7, 9, 10, 7]
    # Variables to process
    var_names = ['vel'] #['draft', 'melt', 'temp', 'salt', 'vel']
    # Constants
    sec_per_year = 365*24*3600
    deg2rad = pi/180.0
    # Parameters for missing circle in ROMS grid
    lon_c = 50
    lat_c = -83
    radius = 10.1
    nbdry = -63+90
    # ROMS vertical grid parameters
    theta_s = 7.0
    theta_b = 2.0
    hc = 250
    N = 31
    # Number of bins in each direction for vector overlay
    num_bins = 30

    print 'Reading ROMS grid'
    # Read the fields we need
    id = Dataset(roms_grid, 'r')
    roms_lon = id.variables['lon_rho'][:,:]
    roms_lat = id.variables['lat_rho'][:,:]
    roms_h = id.variables['h'][:,:]
    roms_mask = id.variables['mask_rho'][:,:]
    roms_zice = id.variables['zice'][:,:]
    roms_angle = id.variables['angle'][:,:]
    id.close()
    # Get land/zice mask
    open_ocn = copy(roms_mask)
    open_ocn[roms_zice!=0] = 0
    land_zice = ma.masked_where(open_ocn==1, open_ocn)
    # Convert grid to spherical coordinates
    roms_x = -(roms_lat+90)*cos(roms_lon*deg2rad+pi/2)
    roms_y = (roms_lat+90)*sin(roms_lon*deg2rad+pi/2)
    # Find centre in spherical coordinates
    x_c = -(lat_c+90)*cos(lon_c*deg2rad+pi/2)
    y_c = (lat_c+90)*sin(lon_c*deg2rad+pi/2)
    # Build a regular x-y grid and select the missing circle
    x_reg_roms, y_reg_roms = meshgrid(linspace(-nbdry, nbdry, num=1000), linspace(-nbdry, nbdry, num=1000))
    land_circle = zeros(shape(x_reg_roms))
    land_circle = ma.masked_where(sqrt((x_reg_roms-x_c)**2 + (y_reg_roms-y_c)**2) > radius, land_circle)

    print 'Building FESOM low-res mesh'
    # Mask open ocean
    elements_lr, mask_patches_lr = make_patches(fesom_mesh_path_lr, circumpolar=True, mask_cavities=True)
    # Unmask ice shelves
    patches_lr = iceshelf_mask(elements_lr)

    print 'Building FESOM high-res mesh'
    elements_hr, mask_patches_hr = make_patches(fesom_mesh_path_hr, circumpolar=True, mask_cavities=True)
    patches_hr = iceshelf_mask(elements_hr)

    for var in var_names:
        print 'Processing variable ' + var

        print 'Reading ROMS fields'
        if var == 'draft':
            # Swap sign on existing zice field; nothing more to read
            roms_data = -1*roms_zice
        else:
            id = Dataset(roms_file, 'r')
            if var == 'melt':
                # Convert from m/s to m/y
                roms_data = id.variables['m'][0,:,:]*sec_per_year
            elif var in ['temp', 'salt']:
                # Bottom layer
                roms_data = id.variables[var][0,0,:,:]
            elif var == 'vel':
                # Read full 3D u and v
                u_3d_tmp = id.variables['u'][0,:,:,:]
                v_3d_tmp = id.variables['v'][0,:,:,:]
                # Get integrands on 3D grid; we only care about dz
                dx, dy, dz, z = cartesian_grid_3d(roms_lon, roms_lat, roms_h, roms_zice, theta_s, theta_b, hc, N)
                # Unrotate each vertical level
                u_3d = ma.empty(shape(dz))
                v_3d = ma.empty(shape(dz))
                num_lat_u = size(u_3d_tmp,1)
                num_lon_u = size(u_3d_tmp,2)
                num_lat_v = size(v_3d_tmp,1)
                num_lon_v = size(v_3d_tmp,2)
                for k in range(N):
                    # Extend into land mask before interpolation to rho-grid so
                    # the land mask doesn't change in the final plot
                    for j in range(1,num_lat_u-1):
                        for i in range(1,num_lon_u-1):
                            # Check for masked points
                            if u_3d_tmp[k,j,i] is ma.masked:
                                # Look at 4 neighbours
                                neighbours = ma.array([u_3d_tmp[k,j-1,i], u_3d_tmp[k,j,i-1], u_3d_tmp[k,j+1,i], u_3d_tmp[k,j,i+1]])
                                # Find how many of them are unmasked
                                num_unmasked = MaskedArray.count(neighbours)
                                if num_unmasked > 0:
                                    # There is at least one unmasked neighbour;
                                    # set u_3d_tmp to their average
                                    u_3d_tmp[k,j,i] = sum(neighbours)/num_unmasked
                    # Repeat for v
                    for j in range(1,num_lat_v-1):
                        for i in range(1,num_lon_v-1):
                            if v_3d_tmp[k,j,i] is ma.masked:
                                neighbours = ma.array([v_3d_tmp[k,j-1,i], v_3d_tmp[k,j,i-1], v_3d_tmp[k,j+1,i], v_3d_tmp[k,j,i+1]])
                                num_unmasked = MaskedArray.count(neighbours)
                                if num_unmasked > 0:
                                    v_3d_tmp[k,j,i] = sum(neighbours)/num_unmasked
                    # Interpolate to rho grid and rotate
                    u_k, v_k = rotate_vector_roms(u_3d_tmp[k,:,:], v_3d_tmp[k,:,:], roms_angle)
                    u_3d[k,:,:] = u_k
                    v_3d[k,:,:] = v_k
                # Vertically average u and v
                u_rho = sum(u_3d*dz, axis=0)/sum(dz, axis=0)
                v_rho = sum(v_3d*dz, axis=0)/sum(dz, axis=0)    
                # Get speed
                roms_data = sqrt(u_rho**2 + v_rho**2)
            id.close()
        # Mask the open ocean and land out of the data field
        roms_data = ma.masked_where(roms_zice==0, roms_data)

        print 'Reading FESOM low-res fields'
        if var != 'draft':
            if var == 'melt':
                id = Dataset(fesom_file_lr_i, 'r')
                # Convert from m/s to m/y
                node_data_lr = id.variables['wnet'][0,:]*sec_per_year
            elif var in ['temp', 'salt']:
                id = Dataset(fesom_file_lr_o, 'r')
                # Read full 3D field for now
                node_data_lr = id.variables[var][0,:]
            elif var == 'vel':
                id = Dataset(fesom_file_lr_o, 'r')
                # The overlaid vectors are based on nodes not elements, so many
                # of the fesom_grid data structures fail to apply and we need to
                # read some of the FESOM grid files again.
                # Read the cavity flag for each 2D surface node
                fesom_cavity_lr = []
                f = open(fesom_mesh_path_lr + 'cavity_flag_nod2d.out', 'r')
                for line in f:
                    tmp = int(line)
                    if tmp == 1:
                        fesom_cavity_lr.append(True)
                    elif tmp == 0:
                        fesom_cavity_lr.append(False)
                    else:
                        print 'Problem'
                        #return
                f.close()
                # Save the number of 2D nodes
                fesom_n2d_lr = len(fesom_cavity_lr)
                # Read rotated lat and lon for each node, also depth
                f = open(fesom_mesh_path_lr + 'nod3d.out', 'r')
                f.readline()
                rlon_lr = []
                rlat_lr = []
                node_depth_lr = []
                for line in f:
                    tmp = line.split()
                    lon_tmp = float(tmp[1])
                    lat_tmp = float(tmp[2])
                    node_depth_tmp = -1*float(tmp[3])
                    if lon_tmp < -180:
                        lon_tmp += 360
                    elif lon_tmp > 180:
                        lon_tmp -= 360
                    rlon_lr.append(lon_tmp)
                    rlat_lr.append(lat_tmp)
                    node_depth_lr.append(node_depth_tmp)
                f.close()
                # For lat and lon, only care about the 2D nodes (the first
                # fesom_n2d indices)
                rlon_lr = array(rlon_lr[0:fesom_n2d_lr])
                rlat_lr = array(rlat_lr[0:fesom_n2d_lr])
                node_depth_lr = array(node_depth_lr)
                # Unrotate longitude
                fesom_lon_lr, fesom_lat_lr = unrotate_grid(rlon_lr, rlat_lr)
                # Calculate polar coordinates of each node
                fesom_x_lr = -(fesom_lat_lr+90)*cos(fesom_lon_lr*deg2rad+pi/2)
                fesom_y_lr = (fesom_lat_lr+90)*sin(fesom_lon_lr*deg2rad+pi/2)
                # Read lists of which nodes are directly below which
                f = open(fesom_mesh_path_lr + 'aux3d.out', 'r')
                max_num_layers_lr = int(f.readline())
                node_columns_lr = zeros([fesom_n2d_lr, max_num_layers_lr])
                for n in range(fesom_n2d_lr):
                    for k in range(max_num_layers_lr):
                        node_columns_lr[n,k] = int(f.readline())
                node_columns_lr = node_columns_lr.astype(int)
                f.close()
                # Now we can actually read the data
                # Read full 3D field for both u and v
                node_ur_3d_lr = id.variables['u'][0,:]
                node_vr_3d_lr = id.variables['v'][0,:]
                # Vertically average
                node_ur_lr = zeros(fesom_n2d_lr)
                node_vr_lr = zeros(fesom_n2d_lr)
                for n in range(fesom_n2d_lr):
                    # Integrate udz, vdz, and dz over this water column
                    udz_col = 0
                    vdz_col = 0
                    dz_col = 0
                    for k in range(max_num_layers_lr-1):
                        if node_columns_lr[n,k+1] == -999:
                            # Reached the bottom
                            break
                        # Trapezoidal rule
                        top_id = node_columns_lr[n,k]
                        bot_id = node_columns_lr[n,k+1]
                        dz_tmp = node_depth_lr[bot_id-1] - node_depth_lr[top_id-1]
                        udz_col += 0.5*(node_ur_3d_lr[top_id-1]+node_ur_3d_lr[bot_id-1])*dz_tmp
                        vdz_col += 0.5*(node_vr_3d_lr[top_id-1]+node_vr_3d_lr[bot_id-1])*dz_tmp
                        dz_col += dz_tmp
                    # Convert from integrals to averages
                    node_ur_lr[n] = udz_col/dz_col
                    node_vr_lr[n] = vdz_col/dz_col
                # Unrotate
                node_u_lr, node_v_lr = unrotate_vector(rlon_lr, rlat_lr, node_ur_lr, node_vr_lr)
                # Calculate speed
                node_data_lr = sqrt(node_u_lr**2 + node_v_lr**2)
            id.close()
        # Calculate given field at each element
        fesom_data_lr = []
        for elm in elements_lr:
            # For each element in an ice shelf cavity, append the mean value
            # for the 3 component Nodes
            if elm.cavity:
                if var == 'draft':
                    # Ice shelf draft is depth of surface layer
                    fesom_data_lr.append(mean([elm.nodes[0].depth, elm.nodes[1].depth, elm.nodes[2].depth]))
                elif var in ['melt', 'vel']:
                    # Surface nodes
                    fesom_data_lr.append(mean([node_data_lr[elm.nodes[0].id], node_data_lr[elm.nodes[1].id], node_data_lr[elm.nodes[2].id]]))
                elif var in ['temp', 'salt']:
                    # Bottom nodes
                    fesom_data_lr.append(mean([node_data_lr[elm.nodes[0].find_bottom().id], node_data_lr[elm.nodes[1].find_bottom().id], node_data_lr[elm.nodes[2].find_bottom().id]]))

        print 'Reading FESOM high-res fields'
        # As before
        if var != 'draft':
            if var == 'melt':
                id = Dataset(fesom_file_hr_i, 'r')
                node_data_hr = id.variables['wnet'][0,:]*sec_per_year
            elif var in ['temp', 'salt']:
                id = Dataset(fesom_file_hr_o, 'r')
                node_data_hr = id.variables[var][0,:]
            elif var == 'vel':
                id = Dataset(fesom_file_hr_o, 'r')
                fesom_cavity_hr = []
                f = open(fesom_mesh_path_hr + 'cavity_flag_nod2d.out', 'r')
                for line in f:
                    tmp = int(line)
                    if tmp == 1:
                        fesom_cavity_hr.append(True)
                    elif tmp == 0:
                        fesom_cavity_hr.append(False)
                    else:
                        print 'Problem'
                        #return
                f.close()
                fesom_n2d_hr = len(fesom_cavity_hr)
                f = open(fesom_mesh_path_hr + 'nod3d.out', 'r')
                f.readline()
                rlon_hr = []
                rlat_hr = []
                node_depth_hr = []
                for line in f:
                    tmp = line.split()
                    lon_tmp = float(tmp[1])
                    lat_tmp = float(tmp[2])
                    node_depth_tmp = -1*float(tmp[3])
                    if lon_tmp < -180:
                        lon_tmp += 360
                    elif lon_tmp > 180:
                        lon_tmp -= 360
                    rlon_hr.append(lon_tmp)
                    rlat_hr.append(lat_tmp)
                    node_depth_hr.append(node_depth_tmp)
                f.close()
                rlon_hr = array(rlon_hr[0:fesom_n2d_hr])
                rlat_hr = array(rlat_hr[0:fesom_n2d_hr])
                node_depth_hr = array(node_depth_hr)
                fesom_lon_hr, fesom_lat_hr = unrotate_grid(rlon_hr, rlat_hr)
                fesom_x_hr = -(fesom_lat_hr+90)*cos(fesom_lon_hr*deg2rad+pi/2)
                fesom_y_hr = (fesom_lat_hr+90)*sin(fesom_lon_hr*deg2rad+pi/2)
                f = open(fesom_mesh_path_hr + 'aux3d.out', 'r')
                max_num_layers_hr = int(f.readline())
                node_columns_hr = zeros([fesom_n2d_hr, max_num_layers_hr])
                for n in range(fesom_n2d_hr):
                    for k in range(max_num_layers_hr):
                        node_columns_hr[n,k] = int(f.readline())
                node_columns_hr = node_columns_hr.astype(int)
                f.close()
                node_ur_3d_hr = id.variables['u'][0,:]
                node_vr_3d_hr = id.variables['v'][0,:]
                node_ur_hr = zeros(fesom_n2d_hr)
                node_vr_hr = zeros(fesom_n2d_hr)
                for n in range(fesom_n2d_hr):
                    udz_col = 0
                    vdz_col = 0
                    dz_col = 0
                    for k in range(max_num_layers_hr-1):
                        if node_columns_hr[n,k+1] == -999:
                            break
                        top_id = node_columns_hr[n,k]
                        bot_id = node_columns_hr[n,k+1]
                        dz_tmp = node_depth_hr[bot_id-1] - node_depth_hr[top_id-1]
                        udz_col += 0.5*(node_ur_3d_hr[top_id-1]+node_ur_3d_hr[bot_id-1])*dz_tmp
                        vdz_col += 0.5*(node_vr_3d_hr[top_id-1]+node_vr_3d_hr[bot_id-1])*dz_tmp
                        dz_col += dz_tmp
                    node_ur_hr[n] = udz_col/dz_col
                    node_vr_hr[n] = vdz_col/dz_col
                node_u_hr, node_v_hr = unrotate_vector(rlon_hr, rlat_hr, node_ur_hr, node_vr_hr)
                node_data_hr = sqrt(node_u_hr**2 + node_v_hr**2)
            id.close()
        fesom_data_hr = []
        for elm in elements_hr:
            if elm.cavity:
                if var == 'draft':
                    fesom_data_hr.append(mean([elm.nodes[0].depth, elm.nodes[1].depth, elm.nodes[2].depth]))
                elif var in ['melt', 'vel']:
                    fesom_data_hr.append(mean([node_data_hr[elm.nodes[0].id], node_data_hr[elm.nodes[1].id], node_data_hr[elm.nodes[2].id]]))
                elif var in ['temp', 'salt']:
                    fesom_data_hr.append(mean([node_data_hr[elm.nodes[0].find_bottom().id], node_data_hr[elm.nodes[1].find_bottom().id], node_data_hr[elm.nodes[2].find_bottom().id]]))

    # Loop over regions
    for index in range(num_regions):
        print 'Processing ' + region_names[index]
        # Set up a grey square for FESOM to fill the background with land
        x_reg_fesom, y_reg_fesom = meshgrid(linspace(x_min[index], x_max[index], num=100), linspace(y_min[index], y_max[index], num=100))
        land_square = zeros(shape(x_reg_fesom))
        # Find bounds on variable in this region, for both ROMS and FESOM
        # Start with ROMS
        loc = (roms_x >= x_min[index])*(roms_x <= x_max[index])*(roms_y >= y_min[index])*(roms_y <= y_max[index])
        var_min = amin(roms_data[loc])
        var_max = amax(roms_data[loc])
        # Modify with FESOM
        # Low-res
        i = 0
        for elm in elements_lr:
            if elm.cavity:
                if any(elm.x >= x_min[index]) and any(elm.x <= x_max[index]) and any(elm.y >= y_min[index]) and any(elm.y <= y_max[index]):
                    if fesom_data_lr[i] < var_min:
                        var_min = fesom_data_lr[i]
                    if fesom_data_lr[i] > var_max:
                        var_max = fesom_data_lr[i]
                i += 1
        # High-res
        i = 0
        for elm in elements_hr:
            if elm.cavity:
                if any(elm.x >= x_min[index]) and any(elm.x <= x_max[index]) and any(elm.y >= y_min[index]) and any(elm.y <= y_max[index]):
                    if fesom_data_hr[i] < var_min:
                        var_min = fesom_data_hr[i]
                    if fesom_data_hr[i] > var_max:
                        var_max = fesom_data_hr[i]
                i += 1
        if var == 'melt':
            # Special colour map
            if var_min < 0:
                # There is refreezing here; include blue for elements < 0
                cmap_vals = array([var_min, 0, 0.25*var_max, 0.5*var_max, 0.75*var_max, var_max])
                cmap_colors = [(0.26, 0.45, 0.86), (1, 1, 1), (1, 0.9, 0.4), (0.99, 0.59, 0.18), (0.5, 0.0, 0.08), (0.96, 0.17, 0.89)]
                cmap_vals_norm = (cmap_vals - var_min)/(var_max - var_min)
                cmap_list = []
                for i in range(size(cmap_vals)):
                    cmap_list.append((cmap_vals_norm[i], cmap_colors[i]))
                mf_cmap = LinearSegmentedColormap.from_list('melt_freeze', cmap_list)
            else:
                # No refreezing
                cmap_vals = array([0, 0.25*var_max, 0.5*var_max, 0.75*var_max, var_max])
                cmap_colors = [(1, 1, 1), (1, 0.9, 0.4), (0.99, 0.59, 0.18), (0.5, 0.0, 0.08), (0.96, 0.17, 0.89)]
                cmap_vals_norm = cmap_vals/var_max
                cmap_list = []
                for i in range(size(cmap_vals)):
                    cmap_list.append((cmap_vals_norm[i], cmap_colors[i]))
                mf_cmap = LinearSegmentedColormap.from_list('melt_freeze', cmap_list)
            colour_map = mf_cmap            
        elif var == 'vel':
            colour_map = 'cool'
        else:
            colour_map = 'jet'
        if var == 'vel':
            # Make vectors for overlay
            # Set up bins (edges)
            x_bins = linspace(x_min[index], x_max[index], num=num_bins+1)
            y_bins = linspace(y_min[index], y_max[index], num=num_bins+1)
            # Calculate centres of bins (for plotting)
            x_centres = 0.5*(x_bins[:-1] + x_bins[1:])
            y_centres = 0.5*(y_bins[:-1] + y_bins[1:])
            # ROMS
            # First set up arrays to integrate velocity in each bin
            # Simple averaging of all the points inside each bin
            roms_u = zeros([size(y_centres), size(x_centres)])
            roms_v = zeros([size(y_centres), size(x_centres)])
            roms_num_pts = zeros([size(y_centres), size(x_centres)])
            # First convert to polar coordinates, rotate to account for
            # longitude in circumpolar projection, and convert back to vector
            # components
            theta_roms = arctan2(v_rho, u_rho)
            theta_circ_roms = theta_roms - roms_lon*deg2rad
            u_circ_roms = roms_data*cos(theta_circ_roms) # roms_data is speed
            v_circ_roms = roms_data*sin(theta_circ_roms)
            # Loop over all points (can't find a better way to do this)
            for j in range(size(roms_data,0)):
                for i in range(size(roms_data,1)):
                    # Make sure data isn't masked (i.e. land or open ocean)
                    if u_circ_roms[j,i] is not ma.masked:
                        # Check if we're in the region of interest
                        if roms_x[j,i] > x_min[index] and roms_x[j,i] < x_max[index] and roms_y[j,i] > y_min[index] and roms_y[j,i] < y_max[index]:
                            # Figure out which bins this falls into
                            x_index = nonzero(x_bins > roms_x[j,i])[0][0]-1
                            y_index = nonzero(y_bins > roms_y[j,i])[0][0]-1
                            # Integrate
                            roms_u[y_index, x_index] += u_circ_roms[j,i]
                            roms_v[y_index, x_index] += v_circ_roms[j,i]
                            roms_num_pts[y_index, x_index] += 1
            # Convert from sums to averages
            # First mask out points with no data
            roms_u = ma.masked_where(roms_num_pts==0, roms_u)
            roms_v = ma.masked_where(roms_num_pts==0, roms_v)
            # Divide everything else by the number of points
            flag = roms_num_pts > 0
            roms_u[flag] = roms_u[flag]/roms_num_pts[flag]
            roms_v[flag] = roms_v[flag]/roms_num_pts[flag]
            # FESOM low-res
            fesom_u_lr = zeros([size(y_centres), size(x_centres)])
            fesom_v_lr = zeros([size(y_centres), size(x_centres)])
            fesom_num_pts_lr = zeros([size(y_centres), size(x_centres)])
            theta_fesom_lr = arctan2(node_v_lr, node_u_lr)
            theta_circ_fesom_lr = theta_fesom_lr - fesom_lon_lr*deg2rad
            u_circ_fesom_lr = node_data_lr*cos(theta_circ_fesom_lr) # node_data is speed
            v_circ_fesom_lr = node_data_lr*sin(theta_circ_fesom_lr)
            # Loop over 2D nodes to fill in the velocity bins as before
            for n in range(fesom_n2d_lr):
                if fesom_cavity_lr[n]:
                    if fesom_x_lr[n] > x_min[index] and fesom_x_lr[n] < x_max[index] and fesom_y_lr[n] > y_min[index] and fesom_y_lr[n] < y_max[index]:
                        x_index = nonzero(x_bins > fesom_x_lr[n])[0][0]-1
                        y_index = nonzero(y_bins > fesom_y_lr[n])[0][0]-1
                        fesom_u_lr[y_index, x_index] += u_circ_fesom_lr[n]
                        fesom_v_lr[y_index, x_index] += v_circ_fesom_lr[n]
                        fesom_num_pts_lr[y_index, x_index] += 1
            fesom_u_lr = ma.masked_where(fesom_num_pts_lr==0, fesom_u_lr)
            fesom_v_lr = ma.masked_where(fesom_num_pts_lr==0, fesom_v_lr)
            flag = fesom_num_pts_lr > 0
            fesom_u_lr[flag] = fesom_u_lr[flag]/fesom_num_pts_lr[flag]
            fesom_v_lr[flag] = fesom_v_lr[flag]/fesom_num_pts_lr[flag]
            # FESOM high-res
            fesom_u_hr = zeros([size(y_centres), size(x_centres)])
            fesom_v_hr = zeros([size(y_centres), size(x_centres)])
            fesom_num_pts_hr = zeros([size(y_centres), size(x_centres)])
            theta_fesom_hr = arctan2(node_v_hr, node_u_hr)
            theta_circ_fesom_hr = theta_fesom_hr - fesom_lon_hr*deg2rad
            u_circ_fesom_hr = node_data_hr*cos(theta_circ_fesom_hr) # node_data is speed
            v_circ_fesom_hr = node_data_hr*sin(theta_circ_fesom_hr)
            # Loop over 2D nodes to fill in the velocity bins as before
            for n in range(fesom_n2d_hr):
                if fesom_cavity_hr[n]:
                    if fesom_x_hr[n] > x_min[index] and fesom_x_hr[n] < x_max[index] and fesom_y_hr[n] > y_min[index] and fesom_y_hr[n] < y_max[index]:
                        x_index = nonzero(x_bins > fesom_x_hr[n])[0][0]-1
                        y_index = nonzero(y_bins > fesom_y_hr[n])[0][0]-1
                        fesom_u_hr[y_index, x_index] += u_circ_fesom_hr[n]
                        fesom_v_hr[y_index, x_index] += v_circ_fesom_hr[n]
                        fesom_num_pts_hr[y_index, x_index] += 1
            fesom_u_hr = ma.masked_where(fesom_num_pts_hr==0, fesom_u_hr)
            fesom_v_hr = ma.masked_where(fesom_num_pts_hr==0, fesom_v_hr)
            flag = fesom_num_pts_hr > 0
            fesom_u_hr[flag] = fesom_u_hr[flag]/fesom_num_pts_hr[flag]
            fesom_v_hr[flag] = fesom_v_hr[flag]/fesom_num_pts_hr[flag]
        # Plot
        fig = figure(figsize=(20, ysize[index]))
        fig.patch.set_facecolor('white')
        # MetROMS
        ax = fig.add_subplot(1,3,1, aspect='equal')
        # First shade land and zice in grey
        contourf(roms_x, roms_y, land_zice, 1, colors=(('0.6', '0.6', '0.6')))
        # Fill in the missing circle
        contourf(x_reg_roms, y_reg_roms, land_circle, 1, colors=(('0.6', '0.6', '0.6')))
        # Now shade the data
        pcolor(roms_x, roms_y, roms_data, vmin=var_min, vmax=var_max, cmap=colour_map)
        if var == 'vel':
            # Overlay vectors
            quiver(x_centres, y_centres, roms_u, roms_v, scale=1.5, headwidth=6, headlength=7, color='black')
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        axis('off')
        title('MetROMS', fontsize=24)
        # FESOM low-res
        ax = fig.add_subplot(1,3,2, aspect='equal')
        # Start with land background
        contourf(x_reg_fesom, y_reg_fesom, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        # Add ice shelf elements
        img = PatchCollection(patches_lr, cmap=colour_map)
        img.set_array(array(fesom_data_lr))
        img.set_edgecolor('face')
        img.set_clim(vmin=var_min, vmax=var_max)
        ax.add_collection(img)
        # Mask out the open ocean in white
        overlay = PatchCollection(mask_patches_lr, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        if var == 'vel':
            # Overlay vectors
            quiver(x_centres, y_centres, fesom_u_lr, fesom_v_lr, scale=1.5, headwidth=6, headlength=7, color='black')
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        axis('off')
        title('FESOM (low-res)', fontsize=24)
        # FESOM high-res
        ax = fig.add_subplot(1,3,3, aspect='equal')
        contourf(x_reg_fesom, y_reg_fesom, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(patches_hr, cmap=colour_map)
        img.set_array(array(fesom_data_hr))
        img.set_edgecolor('face')
        img.set_clim(vmin=var_min, vmax=var_max)
        ax.add_collection(img)
        overlay = PatchCollection(mask_patches_hr, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax.add_collection(overlay)
        if var == 'vel':
            # Overlay vectors
            quiver(x_centres, y_centres, fesom_u_hr, fesom_v_hr, scale=1.5, headwidth=6, headlength=7, color='black')
        xlim([x_min[index], x_max[index]])
        ylim([y_min[index], y_max[index]])
        axis('off')
        title('FESOM (high-res)', fontsize=24)
        # Colourbar on the right
        cbaxes = fig.add_axes([0.92, 0.2, 0.01, 0.6])
        cbar = colorbar(img, cax=cbaxes)
        cbar.ax.tick_params(labelsize=20)
        # Main title
        if var == 'draft':
            title_string = ' draft (m)'
        elif var == 'melt':
            title_string = ' melt rate (m/y)'
        elif var == 'temp':
            title_string = r' bottom water temperature ($^{\circ}$C)'
        elif var == 'salt':
            title_string = ' bottom water salinity (psu)'
        elif var == 'vel':
            title_string = ' vertically averaged ocean velocity (m/s)'
        suptitle(region_names[index] + title_string, fontsize=30)
        subplots_adjust(wspace=0.05)
        #fig.show()
        fig.savefig(fig_heads[index] + '_' + var + '.png')


# Command-line interface
if __name__ == "__main__":

    mip_regions_1var ()
