from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.mlab import griddata
from calc_z import *
from interp_lon_roms import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from fesom_grid import *
from fesom_sidegrid import *

# Make a 3x2 plot of temperature (left) and salinity (right) through 0E.
# The top row is the initial conditions from ECCO2. The middle and bottom rows
# are the end of the simulation (last 5 day average) from MetROMS and FESOM
# respectively.
# Input:
# roms_grid = path to ROMS grid file
# roms_file = path to file containing Jan 2016 monthly average of temperature
#             and salinity in ROMS
# fesom_mesh_path = path to FESOM mesh directory
# fesom_file = path to file containing Jan 2016 monthly average of temperature
#              and salinity in FESOM
def mip_drift_slices (roms_grid, roms_file, fesom_mesh_path, fesom_file):

    # Paths to ECCO2 files with initial conditions for temp and salt
    ecco_temp_file = '/short/m68/kaa561/metroms_iceshelf/data/originals/ECCO2/THETA.1440x720x50.199201.nc'
    ecco_salt_file = '/short/m68/kaa561/metroms_iceshelf/data/originals/ECCO2/SALT.1440x720x50.199201.nc'
    # Longitude to interpolate to (OE)
    lon0 = 0
    # Bounds on plot
    lat_min = -73
    lat_max = -30
    depth_min = -6000
    depth_max = 0
    # ROMS grid parameters
    theta_s = 7.0
    theta_b = 2.0
    hc = 250
    N = 31
    # Bounds on colour scales for temperature and salinity
    temp_min = -2
    temp_max = 6
    salt_min = 33.9
    salt_max = 34.9

    # Get longitude for the title
    if lon0 < 0:
        lon_string = str(int(round(-lon0))) + r'$^{\circ}$W'
    else:
        lon_string = str(int(round(lon0))) + r'$^{\circ}$E'

    print 'Processing ECCO2'
    id = Dataset(ecco_temp_file, 'r')
    # Read grid variables
    ecco_lat = id.variables['LATITUDE_T'][:]
    ecco_depth = -1*id.variables['DEPTH_T'][:]
    if lon0 == 0:
        # Hard-coded lon0 = 0E: average between the first (0.125 E) and last
        # (359.875 E = -0.125 W) indices in the regular ECCO2 grid
        ecco_temp = 0.5*(id.variables['THETA'][0,:,:,0] + id.variables['THETA'][0,:,:,-1])
        id.close()
        id = Dataset(ecco_salt_file, 'r')
        ecco_salt = 0.5*(id.variables['SALT'][0,:,:,0] + id.variables['SALT'][0,:,:,-1])
        id.close()
    else:
        print 'lon0 is only coded for 0E at this time'
        return

    print 'Processing ROMS'
    # Read grid variables we need
    id = Dataset(roms_grid, 'r')
    roms_lon_2d = id.variables['lon_rho'][:,:]
    roms_lat_2d = id.variables['lat_rho'][:,:]
    roms_h = id.variables['h'][:,:]
    roms_zice = id.variables['zice'][:,:]
    id.close()
    # Read temperature and salinity
    id = Dataset(roms_file, 'r')
    roms_temp_3d = id.variables['temp'][0,:,:,:]
    roms_salt_3d = id.variables['salt'][0,:,:,:]
    id.close()
    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    roms_z_3d, sc_r, Cs_r = calc_z(roms_h, roms_zice, theta_s, theta_b, hc, N)
    # Make sure we are in the range 0-360
    if lon0 < 0:
        lon0 += 360
    # Interpolate to lon0
    roms_temp, roms_z, roms_lat = interp_lon_roms(roms_temp_3d, roms_z_3d, roms_lat_2d, roms_lon_2d, lon0)
    roms_salt, roms_z, roms_lat = interp_lon_roms(roms_salt_3d, roms_z_3d, roms_lat_2d, roms_lon_2d, lon0)
    # Switch back to range -180-180
    if lon0 > 180:
        lon0 -= 360

    print 'Processing FESOM'
    # Build regular elements
    elements = fesom_grid(fesom_mesh_path)
    # Read temperature and salinity
    id = Dataset(fesom_file, 'r')
    fesom_temp_nodes = id.variables['temp'][0,:]
    fesom_salt_nodes = id.variables['salt'][0,:]
    id.close()
    # Make SideElements
    selements_temp = fesom_sidegrid(elements, fesom_temp_nodes, lon0, lat_max)
    selements_salt = fesom_sidegrid(elements, fesom_salt_nodes, lon0, lat_max)
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

    # Set up axis labels the way we want them
    lat_ticks = arange(lat_min+3, lat_max+10, 10)
    lat_labels = []
    for val in lat_ticks:
        lat_labels.append(str(int(round(-val))) + r'$^{\circ}$S')
    depth_ticks = range(depth_min+1000, 0+1000, 1000)
    depth_labels = []
    for val in depth_ticks:
        depth_labels.append(str(int(round(-val))))

    print 'Plotting'
    fig = figure(figsize=(14,18))
    # ECCO2
    gs1 = GridSpec(1,2)
    gs1.update(left=0.1, right=0.95, bottom=0.69, top=0.93, wspace=0.08)
    # Temperature
    ax = subplot(gs1[0,0])
    pcolor(ecco_lat, ecco_depth, ecco_temp, vmin=temp_min, vmax=temp_max, cmap='jet')
    title(r'Temperature ($^{\circ}$C)', fontsize=24)
    ylabel('Depth (m)', fontsize=18)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels, fontsize=16)
    text(-64, 1000, 'a) ECCO2 initial conditions at ' + lon_string + ', January 1992', fontsize=28)
    # Salinity
    ax = subplot(gs1[0,1])
    pcolor(ecco_lat, ecco_depth, ecco_salt, vmin=salt_min, vmax=salt_max, cmap='jet')
    title('Salinity (psu)', fontsize=24)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels([])
    # MetROMS
    gs2 = GridSpec(1,2)
    gs2.update(left=0.1, right=0.95, bottom=0.38, top=0.62, wspace=0.08)
    # Temperature
    ax = subplot(gs2[0,0])
    pcolor(roms_lat, roms_z, roms_temp, vmin=temp_min, vmax=temp_max, cmap='jet')
    ylabel('Depth (m)', fontsize=18)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels, fontsize=16)
    text(-49, 300, 'b) MetROMS, January 2016', fontsize=28)
    # Salinity
    ax = subplot(gs2[0,1])
    pcolor(roms_lat, roms_z, roms_salt, vmin=salt_min, vmax=salt_max, cmap='jet')
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels([])
    # FESOM
    gs3 = GridSpec(1,2)
    gs3.update(left=0.1, right=0.95, bottom=0.07, top=0.31, wspace=0.08)
    # Temperature
    ax = subplot(gs3[0,0])
    img = PatchCollection(patches, cmap='jet')
    img.set_array(array(fesom_temp))
    img.set_edgecolor('face')
    img.set_clim(vmin=temp_min, vmax=temp_max)
    ax.add_collection(img)
    ylabel('Depth (m)', fontsize=18)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels, fontsize=16)
    text(-53, 300, 'c) FESOM (high-res), January 2016', fontsize=28)
    # Add a colorbar for temperature
    cbaxes = fig.add_axes([0.17, 0.02, 0.3, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='both', ticks=arange(temp_min, temp_max+2, 2))
    cbar.ax.tick_params(labelsize=16)
    # Salinity
    ax = subplot(gs3[0,1])
    img = PatchCollection(patches, cmap='jet')
    img.set_array(array(fesom_salt))
    img.set_edgecolor('face')
    img.set_clim(vmin=salt_min, vmax=salt_max)
    ax.add_collection(img)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels([])
    # Add a colorbar for salinity
    cbaxes = fig.add_axes([0.6, 0.02, 0.3, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='both', ticks=arange(salt_min+0.1, salt_max+0.1, 0.2))
    cbar.ax.tick_params(labelsize=16)
    fig.show()
    fig.savefig('ts_drift.png')


# Command-line interface
if __name__ == "__main__":

    roms_grid = raw_input("Path to ROMS grid file: ")
    roms_file = raw_input("Path to last ocean_avg file output by ROMS: ")
    fesom_mesh_path = raw_input("Path to FESOM mesh directory: ")
    fesom_file = raw_input("Path to last oce.mean file output by FESOM: ")
    mip_drift_slices(roms_grid, roms_file, fesom_mesh_path, fesom_file)
    
    
