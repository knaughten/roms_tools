from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from monthly_avg_cice import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from patches import *
from monthly_avg import *

# Make a 3x2 plot showing February (top) and August (bottom) sea ice
# concentration (1992-2015 average) for NSIDC, MetROMS, and FESOM.
# Input:
# cice_file = path to CICE output file containing 5-day averages for the entire
#             simulation
# fesom_mesh_path = path to FESOM mesh directory
# fesom_output_dir = path to FESOM output directory containing one ice.mean.nc
#                    file for each year (5-day averages)
def mip_aice_minmax_nsidc (cice_file, fesom_mesh_path, fesom_output_dir):

    # Range of years to process (exclude 2016 because no NSIDC data)
    start_year = 1992
    end_year = 2015
    # Beginning of FESOM output filenames
    fesom_file_head = 'MK44005'
    # Paths to reconstruct NSIDC files for each month
    nsidc_head1 = '/short/m68/kaa561/nsidc_aice/seaice_conc_monthly_sh_f11_'
    nsidc_head2 = '/short/m68/kaa561/nsidc_aice/seaice_conc_monthly_sh_f13_'
    nsidc_head3 = '/short/m68/kaa561/nsidc_aice/seaice_conc_monthly_sh_f17_'
    nsidc_tail = '_v02r00.nc'
    # FESOM plotting parameters
    circumpolar = True
    mask_cavities = True
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    num_years = end_year - start_year + 1

    print 'Processing NSIDC'
    # First read the grid
    id = Dataset(nsidc_head1 + '199201' + nsidc_tail, 'r')
    nsidc_lon = id.variables['longitude'][:,:]
    nsidc_lat = id.variables['latitude'][:,:]
    id.close()
    # Read February and August concentration for each year
    nsidc_feb = ma.empty([num_years, size(nsidc_lon,0), size(nsidc_lat,1)])
    nsidc_aug = ma.empty([num_years, size(nsidc_lon,0), size(nsidc_lat,1)])
    # Loop over years
    for year in range(start_year, end_year+1):
        # Reconstruct file paths
        if year < 1996:
            feb_file = nsidc_head1 + str(year) + '02' + nsidc_tail
            aug_file = nsidc_head1 + str(year) + '08' + nsidc_tail
        elif year < 2008:
            feb_file = nsidc_head2 + str(year) + '02' + nsidc_tail
            aug_file = nsidc_head2 + str(year) + '08' + nsidc_tail
        else:
            feb_file = nsidc_head3 + str(year) + '02' + nsidc_tail
            aug_file = nsidc_head3 + str(year) + '08' + nsidc_tail
        # Read February data and mask
        id = Dataset(feb_file, 'r')
        feb_data_tmp = id.variables['seaice_conc_monthly_cdr'][0,:,:]
        # (This variable is masked but sea ice concentration isn't)
        nsidc_mask = id.variables['stdev_of_seaice_conc_monthly_cdr'][0,:,:]
        id.close()
        # Apply mask
        feb_data = ma.empty(shape(feb_data_tmp))
        feb_data[:,:] = 0.0
        feb_data[~nsidc_mask.mask] = feb_data_tmp[~nsidc_mask.mask]
        feb_data[nsidc_mask.mask] = ma.masked
        # Save result
        nsidc_feb[year-start_year,:,:] = feb_data[:,:]
        # Repeat for August
        id = Dataset(aug_file, 'r')
        aug_data_tmp = id.variables['seaice_conc_monthly_cdr'][0,:,:]
        nsidc_mask = id.variables['stdev_of_seaice_conc_monthly_cdr'][0,:,:]
        id.close()
        aug_data = ma.empty(shape(aug_data_tmp))
        aug_data[:,:] = 0.0
        aug_data[~nsidc_mask.mask] = aug_data_tmp[~nsidc_mask.mask]
        aug_data[nsidc_mask.mask] = ma.masked
        nsidc_aug[year-start_year,:,:] = aug_data[:,:]
    # Average over years
    nsidc_feb = mean(nsidc_feb, axis=0)
    nsidc_aug = mean(nsidc_aug, axis=0)
    # Make sure mask is still there
    nsidc_feb[nsidc_mask.mask] = ma.masked
    nsidc_aug[nsidc_mask.mask] = ma.masked
    # Polar coordinates for plotting
    nsidc_x = -(nsidc_lat+90)*cos(nsidc_lon*deg2rad+pi/2)
    nsidc_y = (nsidc_lat+90)*sin(nsidc_lon*deg2rad+pi/2)
    # Choose boundaries based on extent of NSIDC grid
    bdry1 = amax(nsidc_x[:,0])
    bdry2 = amin(nsidc_x[:,-1])
    bdry3 = amin(nsidc_y[:,0])
    bdry4 = amax(nsidc_y[:,-1])

    print 'Processing MetROMS'
    # First read the grid
    id = Dataset(cice_file, 'r')
    cice_lon_tmp = id.variables['TLON'][:,:]
    cice_lat_tmp = id.variables['TLAT'][:,:]
    id.close()
    # Wrap the periodic boundary by 1 cell
    cice_lon = ma.empty([size(cice_lon_tmp,0), size(cice_lon_tmp,1)+1])
    cice_lat = ma.empty([size(cice_lat_tmp,0), size(cice_lat_tmp,1)+1])
    cice_lon[:,:-1] = cice_lon_tmp
    cice_lon[:,-1] = cice_lon_tmp[:,0]
    cice_lat[:,:-1] = cice_lat_tmp
    cice_lat[:,-1] = cice_lat_tmp[:,0]
    # Get averages for February and August
    # Start with first year just to initialise the arrays with the right size
    print '...monthly average for ' + str(start_year)
    cice_feb_tmp = monthly_avg_cice(cice_file, 'aice', shape(cice_lon_tmp), 1, instance=1)
    cice_aug_tmp = monthly_avg_cice(cice_file, 'aice', shape(cice_lon_tmp), 7, instance=1)
    # Loop over the rest of the years
    for year in range(start_year+1, end_year+1):
        print '...monthly average for ' + str(year)
        cice_feb_tmp = cice_feb_tmp + monthly_avg_cice(cice_file, 'aice', shape(cice_lon_tmp), 1, instance=year-start_year+1)
        cice_aug_tmp = cice_aug_tmp + monthly_avg_cice(cice_file, 'aice', shape(cice_lon_tmp), 7, instance=year-start_year+1)
    # Convert from integrals to averages
    cice_feb_tmp = cice_feb_tmp/num_years
    cice_aug_tmp = cice_aug_tmp/num_years
    # Wrap the periodic boundary
    cice_feb = ma.empty(shape(cice_lon))
    cice_aug = ma.empty(shape(cice_lon))
    cice_feb[:,:-1] = cice_feb_tmp
    cice_aug[:,:-1] = cice_aug_tmp
    cice_feb[:,-1] = cice_feb_tmp[:,0]
    cice_aug[:,-1] = cice_aug_tmp[:,0]
    # Polar coordinates for plotting
    cice_x = -(cice_lat+90)*cos(cice_lon*deg2rad+pi/2)
    cice_y = (cice_lat+90)*sin(cice_lon*deg2rad+pi/2)

    print 'Processing FESOM'
    # First build the grid
    elements, patches = make_patches(fesom_mesh_path, circumpolar, mask_cavities)
    # Get averages for February and August
    # Start with first year just to initialise the arrays with the right size
    print '...monthly average for ' + str(start_year)
    fesom_feb_nodes = monthly_avg(fesom_output_dir + fesom_file_head + '.' + str(start_year) + '.ice.mean.nc', 'area', 1)
    fesom_aug_nodes = monthly_avg(fesom_output_dir + fesom_file_head + '.' + str(start_year) + '.ice.mean.nc', 'area', 7)
    # Loop over the rest of the years
    for year in range(start_year+1, end_year+1):
        print '...monthly average for ' + str(year)
        fesom_feb_nodes = fesom_feb_nodes + monthly_avg(fesom_output_dir + fesom_file_head + '.' + str(year) + '.ice.mean.nc', 'area', 1)
        fesom_aug_nodes = fesom_aug_nodes + monthly_avg(fesom_output_dir + fesom_file_head + '.' + str(year) + '.ice.mean.nc', 'area', 7)
    # Convert from integrals to averages
    fesom_feb_nodes = fesom_feb_nodes/num_years
    fesom_aug_nodes = fesom_aug_nodes/num_years
    # Find element-averages
    fesom_feb = []
    fesom_aug = []
    for elm in elements:
        if not elm.cavity:
            # Average over 3 component nodes
            fesom_feb.append(mean(array([fesom_feb_nodes[elm.nodes[0].id], fesom_feb_nodes[elm.nodes[1].id], fesom_feb_nodes[elm.nodes[2].id]])))
            fesom_aug.append(mean(array([fesom_aug_nodes[elm.nodes[0].id], fesom_aug_nodes[elm.nodes[1].id], fesom_aug_nodes[elm.nodes[2].id]])))
    fesom_feb = array(fesom_feb)
    fesom_aug = array(fesom_aug)

    print 'Plotting'
    fig = figure(figsize=(16,10))
    # NSIDC, February
    ax = fig.add_subplot(2, 3, 1, aspect='equal')
    img = pcolor(nsidc_x, nsidc_y, nsidc_feb, vmin=0, vmax=1, cmap='jet')
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    title('NSIDC', fontsize=24)
    text(-39, 0, 'February', fontsize=24, ha='right')
    # MetROMS, February
    ax = fig.add_subplot(2, 3, 2, aspect='equal')
    img = pcolor(cice_x, cice_y, cice_feb, vmin=0, vmax=1, cmap='jet')
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    title('MetROMS', fontsize=24)
    # FESOM, February
    ax = fig.add_subplot(2, 3, 3, aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(fesom_feb)
    img.set_clim(vmin=0, vmax=1)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    title('FESOM', fontsize=24)
    # NSIDC, August
    ax = fig.add_subplot(2, 3, 4, aspect='equal')
    img = pcolor(nsidc_x, nsidc_y, nsidc_aug, vmin=0, vmax=1, cmap='jet')
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    text(-39, 0, 'August', fontsize=24, ha='right')
    # MetROMS, August
    ax = fig.add_subplot(2, 3, 5, aspect='equal')
    img = pcolor(cice_x, cice_y, cice_aug, vmin=0, vmax=1, cmap='jet')
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')    
    # FESOM, August
    ax = fig.add_subplot(2, 3, 6, aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(fesom_aug)
    img.set_clim(vmin=0, vmax=1)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    # Add a colourbar at the bottom
    cbaxes = fig.add_axes([0.35, 0.04, 0.3, 0.04])
    cbar = colorbar(img, orientation='horizontal', ticks=arange(0,1+0.25,0.25), cax=cbaxes)
    cbar.ax.tick_params(labelsize=16)
    # Main title
    suptitle('Sea ice concentration ('+str(start_year)+'-'+str(end_year)+' average)', fontsize=30)
    # Make panels closer together
    subplots_adjust(wspace=0.05, hspace=0.05)
    #fig.show()
    fig.savefig('aice_minmax_nsidc.png')


# Command-line interface
if __name__ == "__main__":

    cice_file = raw_input("Path to CICE output file containing data for the entire simulation: ")
    fesom_mesh_path = raw_input("Path to FESOM mesh directory: ")
    fesom_output_dir = raw_input("Path to FESOM output directory containing one ice.mean.nc file for each year: ")
    mip_aice_minmax_nsidc (cice_file, fesom_mesh_path, fesom_output_dir)
    
         
    
        
        
