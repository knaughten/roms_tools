from netCDF4 import Dataset
from numpy import *
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

# For each major ice shelf, make a 2x1 plot of the given field for MetROMS
# (left) and FESOM (right), zoomed into that region. Current options are
# ice shelf draft, ice shelf melt rate, bottom water temperature, bottom water
# salinity, surface velocity (with vectors overlaid), or vertically averaged
# velocity (with vectors overlaid).
# Input:
# var_name = 'draft', 'melt', 'temp', 'salt', 'vsfc', 'vavg'
# roms_grid = if var_name='draft', 'vsfc', or 'vavg', path to the ROMS grid
#             file; otherwise, None
# roms_file = if var_name='draft', None; otherwise, path to the ROMS file
#             containing the given field (m, temp, salt, u and v) time-averaged
#             over the desired period. Must also include the mask_rho and zice
#             fields. You can create this using NCO with (for the melt rate)
#             ncra -v m,mask_rho,zice <input_file_list> output_file.nc
# fesom_mesh_path = path to FESOM mesh directory
# fesom_file = if var_name='draft', None; otherwise, path to the FESOM file
#              containing the given field (wnet, temp, salt, u and v)
#              time-averaged over the same period as ROMS. You can create this
#              using NCO with (for the melt rate)
#              ncra -v wnet <input_file_list> output_file.nc
def mip_cavity_fields (var_name, roms_grid, roms_file, fesom_mesh_path, fesom_file):

    # Name of each ice shelf
    shelf_names = ['Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    # Beginnings of filenames for figures
    fig_heads = ['larsen_d', 'larsen_c', 'wilkins_georgevi_stange', 'ronne_filchner', 'abbot', 'pig', 'thwaites', 'dotson', 'getz', 'nickerson', 'sulzberger', 'mertz', 'totten_moscowuni', 'shackleton', 'west', 'amery', 'prince_harald', 'baudouin_borchgrevink', 'lazarev', 'nivl', 'fimbul_jelbart_ekstrom', 'brunt_riiser_larsen', 'ross']
    # Limits on longitude and latitude for each ice shelf
    # Note Ross crosses 180W=180E
    lon_min = [-62.67, -65.5, -79.17, -85, -104.17, -102.5, -108.33, -114.5, -135.67, -149.17, -155, 144, 115, 94.17, 80.83, 65, 33.83, 19, 12.9, 9.33, -10.05, -28.33, 158.33]
    lon_max = [-59.33, -60, -66.67, -28.33, -88.83, -99.17, -103.33, -111.5, -114.33, -140, -145, 146.62, 123.33, 102.5, 89.17, 75, 37.67, 33.33, 16.17, 12.88, 7.6, -10.33, -146.67]
    lat_min = [-73.03, -69.35, -74.17, -83.5, -73.28, -75.5, -75.5, -75.33, -74.9, -76.42, -78, -67.83, -67.17, -66.67, -67.83, -73.67, -69.83, -71.67, -70.5, -70.75, -71.83, -76.33, -85]
    lat_max = [-69.37, -66.13, -69.5, -74.67, -71.67, -74.17, -74.67, -73.67, -73, -75.17, -76.41, -66.67, -66.5, -64.83, -66.17, -68.33, -68.67, -68.33, -69.33, -69.83, -69.33, -71.5, -77]
    num_shelves = len(shelf_names)

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
    num_bins = 50

    print 'Reading ROMS fields'
    if var_name == 'draft':
        id = Dataset(roms_grid, 'r')
    else:
        id = Dataset(roms_file, 'r')
    roms_lon = id.variables['lon_rho'][:,:]
    roms_lat = id.variables['lat_rho'][:,:]
    roms_mask = id.variables['mask_rho'][:,:]
    roms_zice = id.variables['zice'][:,:]
    if var_name == 'draft':
        # Switch signs
        roms_data = -1*id.variables['zice'][:,:]
    elif var_name == 'melt':
        # Convert from m/s to m/y
        roms_data = id.variables['m'][0,:,:]*sec_per_year
    elif var_name == 'temp':
        # Bottom layer
        roms_data = id.variables['temp'][0,0,:,:]
    elif var_name == 'salt':
        # Bottom layer
        roms_data = id.variables['salt'][0,0,:,:]
    elif var_name in ['vsfc', 'vavg']:
        # Get angle from the grid file
        id2 = Dataset(roms_grid, 'r')
        angle = id2.variables['angle'][:,:]
        id2.close()
        if var_name == 'vsfc':
            # Read surface u and v
            u_tmp = id.variables['u'][0,0,:,:]
            v_tmp = id.variables['v'][0,0,:,:]
            # Interpolate to rho grid and unrotate
            u_rho, v_rho = rotate_vector_roms(u_tmp, v_tmp, angle)
        elif var_name == 'vavg':
            # Read full 3D u and v
            u_3d_tmp = id.variables['u'][0,:,:,:]
            v_3d_tmp = id.variables['v'][0,:,:,:]
            # Read bathymetry from grid file
            id2 = Dataset(roms_grid, 'r')
            roms_h = id2.variables['h'][:,:]
            id2.close()
            # Get integrands on 3D grid; we only care about dz
            dx, dy, dz, z = cartesian_grid_3d(roms_lon, roms_lat, roms_h, roms_zice, theta_s, theta_b, hc, N)
            # Unrotate each vertical level
            u_3d = ma.empty(shape(dz))
            v_3d = ma.empty(shape(dz))
            for k in range(N):
                u_k, v_k = rotate_vector_roms(u_3d_tmp[k,:,:], v_3d_tmp[k,:,:], angle)
                u_3d[k,:,:] = u_k
                v_3d[k,:,:] = v_k
            # Vertically average u and v
            u_rho = sum(u_3d*dz, axis=0)/sum(dz, axis=0)
            v_rho = sum(v_3d*dz, axis=0)/sum(dz, axis=0)          
        # Get speed
        roms_data = sqrt(u_rho**2 + v_rho**2)
    id.close()
    # Get land/zice mask
    open_ocn = copy(roms_mask)
    open_ocn[roms_zice!=0] = 0
    land_zice = ma.masked_where(open_ocn==1, open_ocn)
    # Mask the open ocean and land out of the data field
    roms_data = ma.masked_where(roms_zice==0, roms_data)
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

    print 'Reading FESOM fields'
    # Mask open ocean
    elements, mask_patches = make_patches(fesom_mesh_path, circumpolar=True, mask_cavities=True)
    # Unmask ice shelves
    patches = iceshelf_mask(elements)
    if var_name == 'draft':
        # Nothing more to read
        pass
    else:
        id = Dataset(fesom_file, 'r')
        if var_name == 'melt':
            # Convert from m/s to m/y
            node_data = id.variables['wnet'][0,:]*sec_per_year
        elif var_name == 'temp':
            # Read full 3D field for now
            node_data = id.variables['temp'][0,:]
        elif var_name == 'salt':
            # Read full 3D field for now
            node_data = id.variables['salt'][0,:]
        elif var_name in ['vsfc', 'vavg']:
            # The overlaid vectors are based on nodes not elements, so many
            # of the fesom_grid data structures fail to apply and we need to
            # read some of the FESOM grid files again.
            # Read the cavity flag for each 2D surface node
            fesom_cavity = []
            f = open(fesom_mesh_path + 'cavity_flag_nod2d.out', 'r')
            for line in f:
                tmp = int(line)
                if tmp == 1:
                    fesom_cavity.append(True)
                elif tmp == 0:
                    fesom_cavity.append(False)
                else:
                    print 'Problem'
                    #return
            f.close()
            # Save the number of 2D nodes
            fesom_n2d = len(fesom_cavity)
            # Reed rotated lat and lon for each node; also read depth which is
            # needed for vertically averaged velocity
            f = open(fesom_mesh_path + 'nod3d.out', 'r')
            f.readline()
            rlon = []
            rlat = []
            node_depth = []
            for line in f:
                tmp = line.split()
                lon_tmp = float(tmp[1])
                lat_tmp = float(tmp[2])
                node_depth_tmp = -1*float(tmp[3])
                if lon_tmp < -180:
                    lon_tmp += 360
                elif lon_tmp > 180:
                    lon_tmp -= 360
                rlon.append(lon_tmp)
                rlat.append(lat_tmp)
                node_depth.append(node_depth_tmp)
            f.close()
            # For lat and lon, only care about the 2D nodes (the first
            # fesom_n2d indices)
            rlon = array(rlon[0:fesom_n2d])
            rlat = array(rlat[0:fesom_n2d])
            node_depth = array(node_depth)
            # Unrotate longitude
            fesom_lon, fesom_lat = unrotate_grid(rlon, rlat)
            # Calculate polar coordinates of each node
            fesom_x = -(fesom_lat+90)*cos(fesom_lon*deg2rad+pi/2)
            fesom_y = (fesom_lat+90)*sin(fesom_lon*deg2rad+pi/2)
            if var_name == 'vavg':
                # Read lists of which nodes are directly below which
                f = open(fesom_mesh_path + 'aux3d.out', 'r')
                max_num_layers = int(f.readline())
                node_columns = zeros([fesom_n2d, max_num_layers])
                for n in range(fesom_n2d):
                    for k in range(max_num_layers):
                        node_columns[n,k] = int(f.readline())
                node_columns = node_columns.astype(int)
                f.close()
            # Now we can actually read the data
            # Read full 3D field for both u and v
            node_ur_3d = id.variables['u'][0,:]
            node_vr_3d = id.variables['v'][0,:]
            if var_name == 'vsfc':
                # Only care about the first fesom_n2d nodes (surface)
                node_ur = node_ur_3d[0:fesom_n2d]
                node_vr = node_vr_3d[0:fesom_n2d]
            elif var_name == 'vavg':
                # Vertically average
                node_ur = zeros(fesom_n2d)
                node_vr = zeros(fesom_n2d)
                for n in range(fesom_n2d):
                    # Integrate udz, vdz, and dz over this water column
                    udz_col = 0
                    vdz_col = 0
                    dz_col = 0
                    for k in range(max_num_layers-1):
                        if node_columns[n,k+1] == -999:
                            # Reached the bottom
                            break
                        # Trapezoidal rule
                        top_id = node_columns[n,k]
                        bot_id = node_columns[n,k+1]
                        dz_tmp = node_depth[bot_id-1] - node_depth[top_id-1]
                        udz_col += 0.5*(node_ur_3d[top_id-1]+node_ur_3d[bot_id-1])*dz_tmp
                        vdz_col += 0.5*(node_vr_3d[top_id-1]+node_vr_3d[bot_id-1])*dz_tmp
                        dz_col += dz_tmp
                    # Convert from integrals to averages
                    node_ur[n] = udz_col/dz_col
                    node_vr[n] = vdz_col/dz_col
            # Unrotate
            node_u, node_v = unrotate_vector(rlon, rlat, node_ur, node_vr)
            # Calculate speed
            node_data = sqrt(node_u**2 + node_v**2)
        id.close()
    # Calculate given field at each element
    fesom_data = []
    for elm in elements:
        # For each element in an ice shelf cavity, append the mean value
        # for the 3 component Nodes
        if elm.cavity:
            if var_name == 'draft':
                # Ice shelf draft is depth of surface layer
                fesom_data.append(mean([elm.nodes[0].depth, elm.nodes[1].depth, elm.nodes[2].depth]))
            elif var_name in ['melt', 'vsfc', 'vavg']:
                # Surface nodes (or 2D in the case of vavg)
                fesom_data.append(mean([node_data[elm.nodes[0].id], node_data[elm.nodes[1].id], node_data[elm.nodes[2].id]]))
            elif var_name in ['temp', 'salt']:
                # Bottom nodes
                fesom_data.append(mean([node_data[elm.nodes[0].find_bottom().id], node_data[elm.nodes[1].find_bottom().id], node_data[elm.nodes[2].find_bottom().id]]))

    # Loop over ice shelves
    for index in range(num_shelves):
        print 'Processing ' + shelf_names[index]
        # Convert lat/lon bounds to polar coordinates for plotting
        x1 = -(lat_min[index]+90)*cos(lon_min[index]*deg2rad+pi/2)
        y1 = (lat_min[index]+90)*sin(lon_min[index]*deg2rad+pi/2)
        x2 = -(lat_min[index]+90)*cos(lon_max[index]*deg2rad+pi/2)
        y2 = (lat_min[index]+90)*sin(lon_max[index]*deg2rad+pi/2)
        x3 = -(lat_max[index]+90)*cos(lon_min[index]*deg2rad+pi/2)
        y3 = (lat_max[index]+90)*sin(lon_min[index]*deg2rad+pi/2)
        x4 = -(lat_max[index]+90)*cos(lon_max[index]*deg2rad+pi/2)
        y4 = (lat_max[index]+90)*sin(lon_max[index]*deg2rad+pi/2)
        # Find the new bounds on x and y
        x_min = amin(array([x1, x2, x3, x4]))
        x_max = amax(array([x1, x2, x3, x4]))
        y_min = amin(array([y1, y2, y3, y4]))
        y_max = amax(array([y1, y2, y3, y4]))
        # Now make the plot square: enlarge the smaller of delta_x and delta_y
        # so they are equal
        delta_x = x_max - x_min
        delta_y = y_max - y_min
        if delta_x > delta_y:
            diff = 0.5*(delta_x - delta_y)
            y_min -= diff
            y_max += diff
        elif delta_y > delta_x:
            diff = 0.5*(delta_y - delta_x)
            x_min -= diff
            x_max += diff
        # Set up a grey square for FESOM to fill the background with land
        x_reg_fesom, y_reg_fesom = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
        land_square = zeros(shape(x_reg_fesom))
        # Find bounds on variable in this region, for both ROMS and FESOM
        # Start with ROMS
        loc = (roms_x >= x_min)*(roms_x <= x_max)*(roms_y >= y_min)*(roms_y <= y_max)
        var_min = amin(roms_data[loc])
        var_max = amax(roms_data[loc])
        # Modify with FESOM
        i = 0
        for elm in elements:
            if elm.cavity:
                if any(elm.x >= x_min) and any(elm.x <= x_max) and any(elm.y >= y_min) and any(elm.y <= y_max):
                    if fesom_data[i] < var_min:
                        var_min = fesom_data[i]
                    if fesom_data[i] > var_max:
                        var_max = fesom_data[i]
                i += 1
        if var_name == 'melt':
            # Special colour map
            if var_min < 0:
                # There is refreezing here; include blue for elements below 0
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
        else:
            colour_map = 'jet'
        if var_name in ['vsfc', 'vavg']:
            # Make vectors for overlay
            # Set up bins (edges)
            x_bins = linspace(x_min, x_max, num=num_bins+1)
            y_bins = linspace(y_min, y_max, num=num_bins+1)
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
                    # Make sure data isn't masked
                    if u_circ_roms[j,i] is not ma.masked:
                        # Check if we're in the region of interest
                        if roms_x[j,i] > x_min and roms_x[j,i] < x_max and roms_y[j,i] > y_min and roms_y[j,i] < y_max:
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
            # FESOM
            fesom_u = zeros([size(y_centres), size(x_centres)])
            fesom_v = zeros([size(y_centres), size(x_centres)])
            fesom_num_pts = zeros([size(y_centres), size(x_centres)])
            theta_fesom = arctan2(node_v, node_u)
            theta_circ_fesom = theta_fesom - fesom_lon*deg2rad
            u_circ_fesom = node_data*cos(theta_circ_fesom) # node_data is speed
            v_circ_fesom = node_data*sin(theta_circ_fesom)
            # Loop over 2D nodes to fill in the velocity bins as before
            for n in range(fesom_n2d):
                if fesom_cavity[n]:
                    if fesom_x[n] > x_min and fesom_x[n] < x_max and fesom_y[n] > y_min and fesom_y[n] < y_max:
                        x_index = nonzero(x_bins > fesom_x[n])[0][0]-1
                        y_index = nonzero(y_bins > fesom_y[n])[0][0]-1
                        fesom_u[y_index, x_index] += u_circ_fesom[n]
                        fesom_v[y_index, x_index] += v_circ_fesom[n]
                        fesom_num_pts[y_index, x_index] += 1
            fesom_u = ma.masked_where(fesom_num_pts==0, fesom_u)
            fesom_v = ma.masked_where(fesom_num_pts==0, fesom_v)
            flag = fesom_num_pts > 0
            fesom_u[flag] = fesom_u[flag]/fesom_num_pts[flag]
            fesom_v[flag] = fesom_v[flag]/fesom_num_pts[flag]            
        # Plot
        fig = figure(figsize=(30,12))
        fig.patch.set_facecolor('white')
        # ROMS
        ax1 = fig.add_subplot(1,2,1, aspect='equal')
        # First shade land and zice in grey
        contourf(roms_x, roms_y, land_zice, 1, colors=(('0.6', '0.6', '0.6')))
        # Fill in the missing circle
        contourf(x_reg_roms, y_reg_roms, land_circle, 1, colors=(('0.6', '0.6', '0.6')))
        # Now shade the data
        pcolor(roms_x, roms_y, roms_data, vmin=var_min, vmax=var_max, cmap=colour_map)
        if var_name in ['vsfc', 'vavg']:
            # Overlay vectors
            quiver(x_centres, y_centres, roms_u, roms_v, scale=1.5, color='black')
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        axis('off')
        title('MetROMS', fontsize=24)
        # FESOM
        ax2 = fig.add_subplot(1,2,2, aspect='equal')
        # Start with land background
        contourf(x_reg_fesom, y_reg_fesom, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        # Add ice shelf elements
        img = PatchCollection(patches, cmap=colour_map)
        img.set_array(array(fesom_data))
        img.set_edgecolor('face')
        img.set_clim(vmin=var_min, vmax=var_max)
        ax2.add_collection(img)
        # Mask out the open ocean in white
        overlay = PatchCollection(mask_patches, facecolor=(1,1,1))
        overlay.set_edgecolor('face')
        ax2.add_collection(overlay)
        if var_name in ['vsfc', 'vavg']:
            quiver(x_centres, y_centres, fesom_u, fesom_v, scale=1.5, color='black')
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        axis('off')
        title('FESOM', fontsize=24)
        # Colourbar on the right
        cbaxes = fig.add_axes([0.92, 0.2, 0.01, 0.6])
        cbar = colorbar(img, cax=cbaxes)
        cbar.ax.tick_params(labelsize=20)
        # Main title
        if var_name == 'draft':
            title_string = ' draft (m)'
        elif var_name == 'melt':
            title_string = ' melt rate (m/y)'
        elif var_name == 'temp':
            title_string = r' bottom water temperature ($^{\circ}$C)'
        elif var_name == 'salt':
            title_string = ' bottom water salinity (psu)'
        elif var_name == 'vsfc':
            title_string = ' surface velocity (m/s)'
        elif var_name == 'vavg':
            title_string = ' vertically averaged velocity (m/s)'
        suptitle(shelf_names[index] + title_string, fontsize=30)
        subplots_adjust(wspace=0.05)
        #fig.show()
        fig.savefig(fig_heads[index] + '_' + var_name + '.png')


# Command-line interface
if __name__ == "__main__":

    var_key = int(raw_input("Ice shelf draft (1), melt rate (2), bottom water temperature (3), bottom water salinity (4), surface velocity (5), or vertically averaged velocity (6)? "))
    if var_key == 1:
        var_name = 'draft'
    elif var_key == 2:
        var_name = 'melt'
    elif var_key == 3:
        var_name = 'temp'
    elif var_key == 4:
        var_name = 'salt'
    elif var_key == 5:
        var_name = 'vsfc'
    elif var_key == 6:
        var_name = 'vavg'
    if var_name in ['draft', 'vsfc', 'vavg']:
        roms_grid = raw_input("Path to ROMS grid file: ")
    else:
        roms_grid = None
    if var_name == 'draft':
        roms_file = None
        fesom_file = None
    else:
        roms_file = raw_input("Path to ROMS time-averaged file: ")
        fesom_file = raw_input("Path to FESOM time-averaged file: ")
    fesom_mesh_path = raw_input("Path to FESOM mesh directory: ")
    mip_cavity_fields(var_name, roms_grid, roms_file, fesom_mesh_path, fesom_file)
        
            
                
     
    
