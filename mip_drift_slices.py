from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from calc_z import *
from interp_lon_roms import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from fesom_grid import *
from fesom_sidegrid import *

# Make a 4x2 plot showing zonal slices through 0E of temperature (left) and
# salinity (right), in MetROMS (top) and FESOM (bottom). The top 2x2 plots show
# the average over the last year of simulation, while the bottom 2x2 plots show
# anomalies with respect to the first year.
# Input:
# roms_grid = path to ROMS grid file
# roms_file_first, roms_file_last = paths to files containing ROMS temperature
#                                   and salinity averaged over the first year
#                                   and the last year respectively
# fesom_mesh_path = path to FESOM mesh directory
# fesom_file_first, fesom_file_last = paths to files containing FESOM
#                                     temperature and salinity averaged over
#                                     the first year and the last year
#                                     respectively
def mip_drift_slices (roms_grid, roms_file_first, roms_file_last, fesom_mesh_path, fesom_file_first, fesom_file_last):

    # Longitude to plot (0E)
    lon0 = 0
    # Bounds on plot
    lat_min = -73
    lat_max = -30
    depth_min = -5300
    depth_max = 0
    # ROMS grid parameters
    theta_s = 7.0
    theta_b = 2.0
    hc = 250
    N = 31
    # Bounds on colour scales for temperature and salinity
    temp_min = -2
    temp_max = 6
    temp_anom = 2
    salt_min = 33.9
    salt_max = 34.9
    salt_anom = 0.5

     # Get longitude for the title
    if lon0 < 0:
        lon_string = str(int(round(-lon0))) + r'$^{\circ}$W'
    else:
        lon_string = str(int(round(lon0))) + r'$^{\circ}$E'

    print 'Processing ROMS'
    # Read grid variables we need
    id = Dataset(roms_grid, 'r')
    roms_lon_2d = id.variables['lon_rho'][:,:]
    roms_lat_2d = id.variables['lat_rho'][:,:]
    roms_h = id.variables['h'][:,:]
    roms_zice = id.variables['zice'][:,:]
    id.close()
    # Read temperature and salinity
    id = Dataset(roms_file_first, 'r')
    roms_temp_first_3d = id.variables['temp'][0,:,:,:]
    roms_salt_first_3d = id.variables['salt'][0,:,:,:]
    id.close()
    id = Dataset(roms_file_last, 'r')
    roms_temp_last_3d = id.variables['temp'][0,:,:,:]
    roms_salt_last_3d = id.variables['salt'][0,:,:,:]
    id.close()
    # Get anomalies for last year
    roms_temp_anom_3d = roms_temp_last_3d - roms_temp_first_3d
    roms_salt_anom_3d = roms_salt_last_3d - roms_salt_first_3d
    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    roms_z_3d, sc_r, Cs_r = calc_z(roms_h, roms_zice, theta_s, theta_b, hc, N)
    # Make sure we are in the range 0-360
    if lon0 < 0:
        lon0 += 360
    # Interpolate to lon0
    roms_temp, roms_z, roms_lat = interp_lon_roms(roms_temp_last_3d, roms_z_3d, roms_lat_2d, roms_lon_2d, lon0)
    roms_salt, roms_z, roms_lat = interp_lon_roms(roms_salt_last_3d, roms_z_3d, roms_lat_2d, roms_lon_2d, lon0)
    roms_temp_anom, roms_z, roms_lat = interp_lon_roms(roms_temp_anom_3d, roms_z_3d, roms_lat_2d, roms_lon_2d, lon0)
    roms_salt_anom, roms_z, roms_lat = interp_lon_roms(roms_salt_anom_3d, roms_z_3d, roms_lat_2d, roms_lon_2d, lon0)
    # Switch back to range -180-180
    if lon0 > 180:
        lon0 -= 360

    print 'Processing FESOM'
    # Build regular elements
    elements = fesom_grid(fesom_mesh_path)
    # Read temperature and salinity
    id = Dataset(fesom_file_first, 'r')
    fesom_temp_first_nodes = id.variables['temp'][0,:]
    fesom_salt_first_nodes = id.variables['salt'][0,:]
    id.close()
    id = Dataset(fesom_file_last, 'r')
    fesom_temp_last_nodes = id.variables['temp'][0,:]
    fesom_salt_last_nodes = id.variables['salt'][0,:]
    id.close()
    # Get anomalies for last year
    fesom_temp_anom_nodes = fesom_temp_last_nodes - fesom_temp_first_nodes
    fesom_salt_anom_nodes = fesom_salt_last_nodes - fesom_salt_first_nodes
    # Make SideElements
    selements_temp = fesom_sidegrid(elements, fesom_temp_last_nodes, lon0, lat_max)
    selements_salt = fesom_sidegrid(elements, fesom_salt_last_nodes, lon0, lat_max)
    selements_temp_anom = fesom_sidegrid(elements, fesom_temp_anom_nodes, lon0, lat_max)
    selements_salt_anom = fesom_sidegrid(elements, fesom_salt_anom_nodes, lon0, lat_max)
    # Build an array of quadrilateral patches for the plot, and of data values
    # corresponding to each SideElement
    patches = []
    fesom_temp = []
    for selm in selements_temp:
        # Make patch
        coord = transpose(vstack((selm.y, selm.z)))
        patches.append(Polygon(coord, True, linewidth=0.))
        # Save data value
        fesom_temp.append(selm.var)
    # Repeat for the other variables
    fesom_salt = []
    for selm in selements_salt:
        fesom_salt.append(selm.var)
    fesom_temp_anom = []
    for selm in selements_temp_anom:
        fesom_temp_anom.append(selm.var)
    fesom_salt_anom = []
    for selm in selements_salt_anom:
        fesom_salt_anom.append(selm.var)

    # Set up axis labels the way we want them
    lat_ticks = arange(lat_min+3, lat_max+10, 10)
    lat_labels = []
    for val in lat_ticks:
        lat_labels.append(str(int(round(-val))) + r'$^{\circ}$S')
    depth_ticks = range(depth_min+300, 0+1000, 1000)
    depth_labels = []
    for val in depth_ticks:
        depth_labels.append(str(int(round(-val))))
    
    print 'Plotting'
    fig = figure(figsize=(16,25))
    # Absolute values
    gs1 = GridSpec(2,2)
    gs1.update(left=0.1, right=0.95, bottom=0.57, top=0.925, wspace=0.08, hspace=0.2)
    # MetROMS temperature
    ax = subplot(gs1[0,0])
    pcolor(roms_lat, roms_z, roms_temp, vmin=temp_min, vmax=temp_max, cmap='jet')
    title(r'MetROMS temperature ($^{\circ}$)', fontsize=24)    
    ylabel('Depth (m)', fontsize=18)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels([])
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels, fontsize=16)
    text(-38, 1000, 'Last year, ' + lon_string, fontsize=30)
    # MetROMS salinity
    ax = subplot(gs1[0,1])
    pcolor(roms_lat, roms_z, roms_salt, vmin=salt_min, vmax=salt_max, cmap='jet')
    title('MetROMS salinity (psu)', fontsize=24)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels([])
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels([])
    # FESOM temperature
    ax = subplot(gs1[1,0])
    img = PatchCollection(patches, cmap='jet')
    img.set_array(array(fesom_temp))
    img.set_edgecolor('face')
    img.set_clim(vmin=temp_min, vmax=temp_max)
    ax.add_collection(img)
    title(r'FESOM (high-res) temperature ($^{\circ}$C)', fontsize=24)
    xlabel('Latitude', fontsize=18)
    ylabel('Depth (m)', fontsize=18)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels, fontsize=16)
    # Add a colorbar for temperature
    cbaxes = fig.add_axes([0.17, 0.51, 0.3, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='both', ticks=arange(temp_min, temp_max+2, 2))
    cbar.ax.tick_params(labelsize=16)
    # FESOM salinity
    ax = subplot(gs1[1,1])
    img = PatchCollection(patches, cmap='jet')
    img.set_array(array(fesom_salt))
    img.set_edgecolor('face')
    img.set_clim(vmin=salt_min, vmax=salt_max)
    ax.add_collection(img)
    title('FESOM (high-res) salinity (psu)', fontsize=24)
    xlabel('Latitude', fontsize=18)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels([])
    # Add a colorbar for salinity
    cbaxes = fig.add_axes([0.6, 0.51, 0.3, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='both', ticks=arange(salt_min+0.1, salt_max+0.1, 0.2))
    cbar.ax.tick_params(labelsize=16)
    # Anomalies for last year
    gs2 = GridSpec(2,2)
    gs2.update(left=0.1, right=0.95, bottom=0.075, top=0.43, wspace=0.08, hspace=0.2)
    # MetROMS temperature
    ax = subplot(gs2[0,0])
    pcolor(roms_lat, roms_z, roms_temp_anom, vmin=-temp_anom, vmax=temp_anom, cmap='RdBu_r')
    title(r'MetROMS temperature ($^{\circ}$)', fontsize=24)    
    ylabel('Depth (m)', fontsize=18)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels([])
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels, fontsize=16)
    text(-47, 1000, 'Last year minus first year, ' + lon_string, fontsize=30)
    # MetROMS salinity
    ax = subplot(gs2[0,1])
    pcolor(roms_lat, roms_z, roms_salt_anom, vmin=-salt_anom, vmax=salt_anom, cmap='RdBu_r')
    title('MetROMS salinity (psu)', fontsize=24)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels([])
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels([])
    # FESOM temperature
    ax = subplot(gs2[1,0])
    img = PatchCollection(patches, cmap='RdBu_r')
    img.set_array(array(fesom_temp_anom))
    img.set_edgecolor('face')
    img.set_clim(vmin=-temp_anom, vmax=temp_anom)
    ax.add_collection(img)
    title(r'FESOM (high-res) temperature ($^{\circ}$C)', fontsize=24)
    xlabel('Latitude', fontsize=18)
    ylabel('Depth (m)', fontsize=18)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels, fontsize=16)
    # Add a colorbar for temperature
    cbaxes = fig.add_axes([0.17, 0.02, 0.3, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='both', ticks=arange(-temp_anom, temp_anom+1, 1))
    cbar.ax.tick_params(labelsize=16)
    # FESOM salinity
    ax = subplot(gs2[1,1])
    img = PatchCollection(patches, cmap='RdBu_r')
    img.set_array(array(fesom_salt_anom))
    img.set_edgecolor('face')
    img.set_clim(vmin=-salt_anom, vmax=salt_anom)
    ax.add_collection(img)
    title('FESOM (high-res) salinity (psu)', fontsize=24)
    xlabel('Latitude', fontsize=18)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels([])
    # Add a colorbar for salinity
    cbaxes = fig.add_axes([0.6, 0.02, 0.3, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='both', ticks=arange(-salt_anom, salt_anom+0.25, 0.25))
    cbar.ax.tick_params(labelsize=16)
    #fig.show()
    fig.savefig('ts_slices.png')


# Command-line interface
if __name__ == "__main__":

    roms_grid = raw_input("Path to ROMS grid file: ")
    roms_file_first = raw_input("Path to ROMS temperature/salinity file averaged over first year: ")
    roms_file_last = raw_input("Path to ROMS temperature/salinity file averaged over last year: ")
    fesom_mesh_path = raw_input("Path to FESOM mesh directory: ")
    fesom_file_first = raw_input("Path to FESOM temperature/salinity file averaged over first year: ")
    fesom_file_last = raw_input("Path to FESOM temperature/salinity file averaged over last year: ")
    mip_drift_slices(roms_grid, roms_file_first, roms_file_last, fesom_mesh_path, fesom_file_first, fesom_file_last)

    
    

    
