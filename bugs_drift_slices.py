from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from calc_z import *
from interp_lon_roms import *

def bugs_drift_slices (grid_file, upwind_file, akima_file, split_file):

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
    # Contours to overlay
    temp_contour = 0.75
    salt_contour = 34.5

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

    print 'Building ROMS grid'
    id = Dataset(grid_file, 'r')
    lon_2d = id.variables['lon_rho'][:,:]
    lat_2d = id.variables['lat_rho'][:,:]
    h = id.variables['h'][:,:]
    zice = id.variables['zice'][:,:]
    id.close()
    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    z_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N)
    # Make sure we are in the range 0-360
    if lon0 < 0:
        lon0 += 360
    print 'Processing upwind advection'
    id = Dataset(upwind_file, 'r')
    upwind_temp_3d = id.variables['temp'][0,:,:,:]
    upwind_salt_3d = id.variables['salt'][0,:,:,:]
    id.close()
    # Interpolate to lon0
    upwind_temp, z, lat = interp_lon_roms(upwind_temp_3d, z_3d, lat_2d, lon_2d, lon0)
    upwind_salt, z, lat = interp_lon_roms(upwind_salt_3d, z_3d, lat_2d, lon_2d, lon0)
    print 'Processing Akima advection'
    id = Dataset(akima_file, 'r')
    akima_temp_3d = id.variables['temp'][0,:,:,:]
    akima_salt_3d = id.variables['salt'][0,:,:,:]
    id.close()
    akima_temp, z, lat = interp_lon_roms(akima_temp_3d, z_3d, lat_2d, lon_2d, lon0)
    akima_salt, z, lat = interp_lon_roms(akima_salt_3d, z_3d, lat_2d, lon_2d, lon0)    
    print 'Processing split advection'
    id = Dataset(split_file, 'r')
    split_temp_3d = id.variables['temp'][0,:,:,:]
    split_salt_3d = id.variables['salt'][0,:,:,:]
    id.close()
    split_temp, z, lat = interp_lon_roms(split_temp_3d, z_3d, lat_2d, lon_2d, lon0)
    split_salt, z, lat = interp_lon_roms(split_salt_3d, z_3d, lat_2d, lon_2d, lon0)
    # Switch back to range -180-180
    if lon0 > 180:
        lon0 -= 360

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
    fig = figure(figsize=(14,24))
    # ECCO2
    gs1 = GridSpec(1,2)
    gs1.update(left=0.1, right=0.95, bottom=0.7575, top=0.94, wspace=0.08)
    # Temperature
    ax = subplot(gs1[0,0])
    pcolor(ecco_lat, ecco_depth, ecco_temp, vmin=temp_min, vmax=temp_max, cmap='jet')
    # Overlay contour
    contour(ecco_lat, ecco_depth, ecco_temp, levels=[temp_contour], color='black')
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
    contour(ecco_lat, ecco_depth, ecco_salt, levels=[salt_contour], color='black')
    title('Salinity (psu)', fontsize=24)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels([])
    # Upwind advection
    gs2 = GridSpec(1,2)
    gs2.update(left=0.1, right=0.95, bottom=0.525, top=0.7075, wspace=0.08)
    # Temperature
    ax = subplot(gs2[0,0])
    pcolor(lat, z, upwind_temp, vmin=temp_min, vmax=temp_max, cmap='jet')
    contour(lat, z, upwind_temp, levels=[temp_contour], color='black')
    ylabel('Depth (m)', fontsize=18)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels, fontsize=16)
    text(-60, 300, 'b) Upwind third-order advection, January 2016', fontsize=28)
    # Salinity
    ax = subplot(gs2[0,1])
    pcolor(lat, z, upwind_salt, vmin=salt_min, vmax=salt_max, cmap='jet')
    contour(lat, z, upwind_salt, levels=[salt_contour], color='black')
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels([])
    # Akima advection
    gs3 = GridSpec(1,2)
    gs3.update(left=0.1, right=0.95, bottom=0.2925, top=0.475, wspace=0.08)
    # Temperature
    ax = subplot(gs3[0,0])
    pcolor(lat, z, akima_temp, vmin=temp_min, vmax=temp_max, cmap='jet')
    contour(lat, z, akima_temp, levels=[temp_contour], color='black')
    ylabel('Depth (m)', fontsize=18)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels, fontsize=16)
    text(-52, 300, 'c) Akima advection, January 2016', fontsize=28)
    # Salinity
    ax = subplot(gs3[0,1])
    pcolor(lat, z, akima_salt, vmin=salt_min, vmax=salt_max, cmap='jet')
    contour(lat, z, akima_salt, levels=[salt_contour], color='black')
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels([])
    # Split advection
    gs4 = GridSpec(1,2)
    gs4.update(left=0.1, right=0.95, bottom=0.06, top=0.2425, wspace=0.08)
    # Temperature
    ax = subplot(gs4[0,0])
    img = pcolor(lat, z, split_temp, vmin=temp_min, vmax=temp_max, cmap='jet')
    contour(lat, z, split_temp, levels=[temp_contour], color='black')
    ylabel('Depth (m)', fontsize=18)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels, fontsize=16)
    text(-52, 300, 'd) RSUP3 advection, January 2016', fontsize=28)
    # Add a colorbar for temperature
    cbaxes = fig.add_axes([0.17, 0.015, 0.3, 0.015])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='both', ticks=arange(temp_min, temp_max+2, 2))
    cbar.ax.tick_params(labelsize=16)
    # Salinity
    ax = subplot(gs4[0,1])
    img = pcolor(lat, z, split_salt, vmin=salt_min, vmax=salt_max, cmap='jet')
    contour(lat, z, split_salt, levels=[salt_contour], color='black')
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels([])
    # Add a colorbar for salinity
    cbaxes = fig.add_axes([0.6, 0.015, 0.3, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='both', ticks=arange(salt_min+0.1, salt_max+0.1, 0.2))
    cbar.ax.tick_params(labelsize=16)
    fig.show()
    fig.savefig('bugs_drift.png')


# Command-line interface
if __name__ == "__main__":

    grid_file = raw_input("Path to ROMS grid file: ")
    upwind_file = raw_input("Path to upwind advection file containing monthly averaged temperature and salinity for January 2016: ")
    akima_file = raw_input("Path to Akima advection file containing monthly averaged temperature and salinity for January 2016: ")
    split_file = raw_input("Path to split advection file containing monthly averaged temperature and salinity for January 2016: ")
    bugs_drift_slices(grid_file, upwind_file, akima_file, split_file)
