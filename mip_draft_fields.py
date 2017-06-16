from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from patches import *

# For each major ice shelf, make a 2x1 plot of the ice shelf draft for MetROMS
# (left) and FESOM (right), zoomed into that region. 
# Input:
# roms_grid_file = path to ROMS grid file
# fesom_mesh_path = path to FESOM mesh directory
def mip_draft_fields (roms_grid_file, fesom_mesh_path):

    # Name of each ice shelf
    shelf_names = ['Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    # Filenames for figures
    fig_names = ['larsen_d_draft.png', 'larsen_c_draft.png', 'wilkins_georgevi_stange_draft.png', 'ronne_filchner_draft.png', 'abbot_draft.png', 'pig_draft.png', 'thwaites_draft.png', 'dotson_draft.png', 'getz_draft.png', 'nickerson_draft.png', 'sulzberger_draft.png', 'mertz_draft.png', 'totten_moscowuni_draft.png', 'shackleton_draft.png', 'west_draft.png', 'amery_draft.png', 'prince_harald_draft.png', 'baudouin_borchgrevink_draft.png', 'lazarev_draft.png', 'nivl_draft.png', 'fimbul_jelbart_ekstrom_draft.png', 'brunt_riiser_larsen_draft.png', 'ross_draft.png']
    # Limits on longitude and latitude for each ice shelf
    # Note Ross crosses 180W=180E
    lon_min = [-62.67, -65.5, -79.17, -85, -104.17, -102.5, -108.33, -114.5, -135.67, -149.17, -155, 144, 115, 94.17, 80.83, 65, 33.83, 19, 12.9, 9.33, -10.05, -28.33, 158.33]
    lon_max = [-59.33, -60, -66.67, -28.33, -88.83, -99.17, -103.33, -111.5, -114.33, -140, -145, 146.62, 123.33, 102.5, 89.17, 75, 37.67, 33.33, 16.17, 12.88, 7.6, -10.33, -146.67]
    lat_min = [-73.03, -69.35, -74.17, -83.5, -73.28, -75.5, -75.5, -75.33, -74.9, -76.42, -78, -67.83, -67.17, -66.67, -67.83, -73.67, -69.83, -71.67, -70.5, -70.75, -71.83, -76.33, -85]
    lat_max = [-69.37, -66.13, -69.5, -74.67, -71.67, -74.17, -74.67, -73.67, -73, -75.17, -76.41, -66.67, -66.5, -64.83, -66.17, -68.33, -68.67, -68.33, -69.33, -69.83, -69.33, -71.5, -77]
    num_shelves = len(shelf_names)

    # Constants
    deg2rad = pi/180.0
    # Parameters for missing circle in ROMS grid
    lon_c = 50
    lat_c = -83
    radius = 10.1
    nbdry = -63+90

    print 'Processing ROMS'
    id = Dataset(roms_grid_file, 'r')
    roms_lon = id.variables['lon_rho'][:,:]
    roms_lat = id.variables['lat_rho'][:,:]
    roms_mask = id.variables['mask_rho'][:,:]
    roms_draft = -1*id.variables['zice'][:,:]
    id.close()
    # Get land/zice mask
    open_ocn = copy(roms_mask)
    open_ocn[roms_draft!=0] = 0
    land_zice = ma.masked_where(open_ocn==1, open_ocn)
    # Mask the open ocean and land out of the ice shelf draft
    roms_draft = ma.masked_where(roms_draft==0, roms_draft)
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

    print 'Processing FESOM'
    # Mask open ocean
    elements, mask_patches = make_patches(fesom_mesh_path, circumpolar=True, mask_cavities=True)
    # Unmask ice shelves
    patches = iceshelf_mask(elements)
    # Calculate ice shelf draft
    fesom_draft = []
    for elm in elements:
        # For each element in an ice shelf cavity, append the mean value
        # for the 3 component Nodes
        if elm.cavity:
            # Ice shelf draft is depth of surface layer
            fesom_draft.append(mean([elm.nodes[0].depth, elm.nodes[1].depth, elm.nodes[2].depth]))

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
        # Find bounds on ice shelf draft in this region, for both ROMS and FESOM
        # Start with ROMS
        loc = (roms_x >= x_min)*(roms_x <= x_max)*(roms_y >= y_min)*(roms_y <= y_max)
        var_min = amin(roms_draft[loc])
        var_max = amax(roms_draft[loc])
        # Modify with FESOM
        i = 0
        for elm in elements:
            if elm.cavity:
                if any(elm.x >= x_min) and any(elm.x <= x_max) and any(elm.y >= y_min) and any(elm.y <= y_max):
                    if fesom_draft[i] < var_min:
                        var_min = fesom_draft[i]
                    if fesom_draft[i] > var_max:
                        var_max = fesom_draft[i]
                i += 1
        # Plot
        fig = figure(figsize=(30,12))
        fig.patch.set_facecolor('white')
        # ROMS
        ax1 = fig.add_subplot(1,2,1, aspect='equal')
        contourf(roms_x, roms_y, land_zice, 1, colors=(('0.6', '0.6', '0.6')))
        contourf(x_reg_roms, y_reg_roms, land_circle, 1, colors=(('0.6', '0.6', '0.6')))
        pcolor(roms_x, roms_y, roms_draft, vmin=var_min, vmax=var_max, cmap='jet')
        xlim([x_min, x_max])
        ylim([y_min, y_max])
        axis('off')
        title('MetROMS', fontsize=24)
        # FESOM
        ax2 = fig.add_subplot(1,2,2, aspect='equal')
        contourf(x_reg_fesom, y_reg_fesom, land_square, 1, colors=(('0.6', '0.6', '0.6')))
        img = PatchCollection(patches, cmap='jet')
        img.set_array(array(fesom_draft))
        img.set_edgecolor('face')
        img.set_clim(vmin=var_min, vmax=var_max)
        ax2.add_collection(img)
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
        suptitle(shelf_names[index] + ' ice shelf draft (m)', fontsize=30)
        subplots_adjust(wspace=0.05)
        #fig.show()
        fig.savefig(fig_names[index])


# Command-line interface
if __name__ == "__main__":

    roms_grid_file = raw_input("Path to ROMS grid file: ")
    fesom_mesh_path = raw_input("Path to FESOM mesh directory: ")
    mip_draft_fields(roms_grid_file, fesom_mesh_path)
        
            
                
     
    
