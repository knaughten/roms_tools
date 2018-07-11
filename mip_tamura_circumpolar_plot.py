from numpy import *
from netCDF4 import Dataset
from matplotlib.pyplot import *
from matplotlib.collections import PatchCollection
#from matplotlib.colors import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from patches import *

def mip_tamura_circumpolar_plot ():

    # File paths
    # FESOM mesh paths
    fesom_mesh_path_lr = '/short/y99/kaa561/FESOM/mesh/meshA/'
    fesom_mesh_path_hr = '/short/y99/kaa561/FESOM/mesh/meshB/'
    # CICE 1992-2013 mean ice production files (precomputed in calc_ice_prod.py)
    cice_file = '/short/m68/kaa561/metroms_iceshelf/tmproms/run/intercomparison/ice_prod_1992_2013.nc'
    # FESOM 1992-2013 mean ice production files (precomputed in calc_annual_ice_prod.py in fesomtools)
    fesom_lr_file = '/short/y99/kaa561/FESOM/intercomparison_lowres/output/ice_prod_1992_2013.nc'
    fesom_hr_file = '/short/y99/kaa561/FESOM/intercomparison_highres/output/ice_prod_1992_2013.nc'
    # Tamura's 1992-2013 mean ice production (precomputed on desktop with Matlab)
    tamura_file = '/short/m68/kaa561/tamura_1992_2013_monthly_climatology.nc'
    # FESOM plotting parameters
    circumpolar = True
    mask_cavities = True
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Northern boundary for plot: 64S
    nbdry = -64 + 90
    # Maximum value to plot
    plot_bound = 30    

    print 'Processing Tamura obs'
    id = Dataset(tamura_file, 'r')
    tamura_lon = id.variables['longitude'][:,:]
    tamura_lat = id.variables['latitude'][:,:]
    # Read precomputed sea ice formation (already in m/y)
    tamura_data = id.variables['ice_prod'][:,:]
    id.close()
    # Polar coordinates for plotting
    tamura_x = -(tamura_lat+90)*cos(tamura_lon*deg2rad+pi/2)
    tamura_y = (tamura_lat+90)*sin(tamura_lon*deg2rad+pi/2)

    print 'Processing MetROMS'
    id = Dataset(cice_file, 'r')
    cice_lon_tmp = id.variables['TLON'][:,:]
    cice_lat_tmp = id.variables['TLAT'][:,:]
    # Read precomputed sea ice formation (already in m/y)
    cice_data_tmp = id.variables['ice_prod'][:,:]
    id.close()
    # Wrap the periodic boundary by 1 cell
    cice_lon = ma.empty([size(cice_lon_tmp,0), size(cice_lon_tmp,1)+1])
    cice_lat = ma.empty([size(cice_lat_tmp,0), size(cice_lat_tmp,1)+1])
    cice_data = ma.empty([size(cice_data_tmp,0), size(cice_data_tmp,1)+1])
    cice_lon[:,:-1] = cice_lon_tmp
    cice_lon[:,-1] = cice_lon_tmp[:,0]
    cice_lat[:,:-1] = cice_lat_tmp
    cice_lat[:,-1] = cice_lat_tmp[:,0]
    cice_data[:,:-1] = cice_data_tmp
    cice_data[:,-1] = cice_data_tmp[:,0]
    # Polar coordinates for plotting
    cice_x = -(cice_lat+90)*cos(cice_lon*deg2rad+pi/2)
    cice_y = (cice_lat+90)*sin(cice_lon*deg2rad+pi/2)

    print 'Processing low-res FESOM'
    # Build mesh
    elements_lr, patches_lr = make_patches(fesom_mesh_path_lr, circumpolar, mask_cavities)
    # Read precomputed sea ice formation (already in m/y)
    id = Dataset(fesom_lr_file, 'r')
    fesom_data_nodes_lr = id.variables['ice_prod'][:]
    id.close()
    # Find element-averages
    fesom_data_lr = []
    for elm in elements_lr:
        if not elm.cavity:
            fesom_data_lr.append(mean(array([fesom_data_nodes_lr[elm.nodes[0].id], fesom_data_nodes_lr[elm.nodes[1].id], fesom_data_nodes_lr[elm.nodes[2].id]])))
    fesom_data_lr = array(fesom_data_lr)

    print 'Processing high-res FESOM'
    elements_hr, patches_hr = make_patches(fesom_mesh_path_hr, circumpolar, mask_cavities)
    id = Dataset(fesom_hr_file, 'r')
    fesom_data_nodes_hr = id.variables['ice_prod'][:]
    id.close()
    fesom_data_hr = []
    for elm in elements_hr:
        if not elm.cavity:
            fesom_data_hr.append(mean(array([fesom_data_nodes_hr[elm.nodes[0].id], fesom_data_nodes_hr[elm.nodes[1].id], fesom_data_nodes_hr[elm.nodes[2].id]])))
    fesom_data_hr = array(fesom_data_hr)

    #bounds = linspace(0, plot_bound**(1.0/3), num=100)**3
    #norm = BoundaryNorm(boundaries=bounds, ncolors=256)

    print 'Plotting'
    fig = figure(figsize=(10,12))
    gs = GridSpec(2,2)
    gs.update(left=0.05, right=0.95, bottom=0.09, top=0.9, wspace=0.02, hspace=0.1)
    # Tamura obs
    ax = subplot(gs[0,0], aspect='equal')
    img = pcolor(tamura_x, tamura_y, tamura_data, vmin=0, vmax=plot_bound, cmap='jet')
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    ax.set_xticks([])
    ax.set_yticks([])
    title('Observations', fontsize=24)
    # MetROMS
    ax = subplot(gs[0,1], aspect='equal')
    img = pcolor(cice_x, cice_y, cice_data, vmin=0, vmax=plot_bound, cmap='jet') #, norm=norm)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    ax.set_xticks([])
    ax.set_yticks([])
    title('MetROMS', fontsize=24)
    # FESOM low-res
    ax = subplot(gs[1,0], aspect='equal')
    img = PatchCollection(patches_lr, cmap='jet') #, norm=norm)
    img.set_array(fesom_data_lr)
    img.set_clim(0, plot_bound)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    ax.set_xticks([])
    ax.set_yticks([])
    title('FESOM (low-res)', fontsize=24)
    # FESOM high-res
    ax = subplot(gs[1,1], aspect='equal')
    img = PatchCollection(patches_hr, cmap='jet') #, norm=norm)
    img.set_array(fesom_data_hr)
    img.set_clim(0, plot_bound)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    ax.set_xticks([])
    ax.set_yticks([])
    title('FESOM (high-res)', fontsize=24)
    # Colourbar below
    cbaxes = fig.add_axes([0.3, 0.03, 0.4, 0.03])
    cbar = colorbar(img, cax=cbaxes, extend='max', orientation='horizontal')
    cbar.ax.tick_params(labelsize=16)
    # Main title
    suptitle('Sea ice production (m/y), 1992-2013 mean', fontsize=28)
    fig.show()
    fig.savefig('seaice_tamura.png')
    

# Command-line interface
if __name__ == "__main__":

    mip_tamura_circumpolar_plot()
