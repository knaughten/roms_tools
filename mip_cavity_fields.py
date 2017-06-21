from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.colors import LinearSegmentedColormap
from rotate_vector_roms import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from patches import *

# For each major ice shelf, make a 2x1 plot of the given field for MetROMS
# (left) and FESOM (right), zoomed into that region. Current options are
# ice shelf draft, ice shelf melt rate, bottom water temperature, bottom water
# salinity, or surface velocity.
# Input:
# var_name = 'draft', 'melt', 'temp', 'salt', or 'vel'
# roms_file = if var_name='draft', path to the ROMS grid file; otherwise, path
#             to the ROMS file containing the given field (m, temp, salt, u and
#             v) time-averaged over the desired period. Must also include the
#             mask_rho and zice fields. You can create this using NCO with (for
#             the melt rate)
#             ncra -v m,mask_rho,zice <input_file_list> output_file.nc
# fesom_file = if var_name='draft', None; otherwise, path to the FESOM file
#              containing the given field (wnet, temp, salt, u and v)
#              time-averaged over the same period as ROMS. You can create this
#              using NCO with (for the melt rate)
#              ncra -v wnet <input_file_list> output_file.nc
# fesom_mesh_path = path to FESOM mesh directory
def mip_cavity_fields (var_name, roms_file, fesom_file, fesom_mesh_path):

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

    print 'Reading ROMS fields'
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
    elif var_name == 'vel':
        # Read surface u and v
        u_tmp = id.variables['u'][0,0,:,:]
        v_tmp = id.variables['v'][0,0,:,:]
        # Interpolate both to the rho-grid
        w_bdry_u = 0.5*(u_tmp[:,0] + u_tmp[:,-1])
        middle_u = 0.5*(u_tmp[:,0:-1] + u_tmp[:,1:])
        e_bdry_u = w_bdry_u[:]
        u_rho = ma.concatenate((w_bdry_u[:,None], middle_u, e_bdry_u[:,None]), axis=1)
        s_bdry_v = v_tmp[0,:]
        middle_v = 0.5*(v_tmp[0:-1,:] + v_tmp[1:,:])
        n_bdry_v = v_tmp[-1,:]
        v_rho = ma.concatenate((s_bdry_v[None,:], middle_v, n_bdry_v[None,:]), axis=0)
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
        elif var_name == 'vel':
            # Read full 3D field for both u and v
            node_u = id.variables['u'][0,:]
            node_v = id.variables['v'][0,:]
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
            elif var_name in ['melt', 'vel']:
                # Surface nodes
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
        elif var_name == 'vel':
            title_string = ' surface velocity (m/s)'
        suptitle(shelf_names[index] + title_string, fontsize=30)
        subplots_adjust(wspace=0.05)
        #fig.show()
        fig.savefig(fig_heads[index] + '_' + var_name + '.png')


# Command-line interface
if __name__ == "__main__":

    var_key = int(raw_input("Ice shelf draft (1), melt rate (2), bottom water temperature (3), bottom water salinity (4), or surface velocity (5)? "))
    if var_key == 1:
        var_name = 'draft'
        roms_file = raw_input("Path to ROMS grid file: ")
        fesom_file = None
    else:
        roms_file = raw_input("Path to ROMS time-averaged file: ")
        fesom_file = raw_input("Path to FESOM time-averaged file: ")
        if var_key == 2:
            var_name = 'melt'
        elif var_key == 3:
            var_name = 'temp'
        elif var_key == 4:
            var_name = 'salt'
        elif var_key == 5:
            var_name = 'vel'    
    fesom_mesh_path = raw_input("Path to FESOM mesh directory: ")
    mip_cavity_fields(var_name, roms_file, fesom_file, fesom_mesh_path)
        
            
                
     
    
