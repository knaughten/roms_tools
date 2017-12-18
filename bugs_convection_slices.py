from netCDF4 import Dataset, num2date
from numpy import *
from matplotlib.pyplot import *
from calc_z import *
from interp_lon_roms import *

def bugs_convection_slices ():

    # Files and timesteps to plot
    file_beg = '/short/m68/kaa561/metroms_iceshelf/tmproms/run/bug_chapter/no_restoring/ocean_avg_0003.nc'
    tstep_beg = 18
    file_end = '/short/m68/kaa561/metroms_iceshelf/tmproms/run/bug_chapter/no_restoring/ocean_avg_0012.nc'
    tstep_end = 2
    # Longitude to interpolate to
    lon0 = -13
    # Deepest depth to plot
    depth_min = -250
    depth_ticks = arange(depth_min, 0+50, 50)
    depth_labels = []
    for depth in depth_ticks:
        depth_labels.append(str(int(-depth)))
    # Latitudes to plot
    lat_min = -74
    lat_max = -55
    lat_ticks = arange(lat_min+4, lat_max+5, 5)
    lat_labels = []
    for lat in lat_ticks:
        lat_labels.append(str(int(-lat)) + r'$^{\circ}$S')
    # Bounds on colour scales
    var_min = [-2, 33.9]
    var_max = [2, 34.7]
    var_tick = [1, 0.2]
    lev1 = linspace(var_min[0], var_max[0], num=50)
    lev2 = linspace(var_min[1], var_max[1], num=50)
    # Grid parameters
    theta_s = 7.0
    theta_b = 2.0
    hc = 250
    N = 31
    # Month names for titles
    month_names = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']

    # Get longitude for the title
    if lon0 < 0:
        lon_string = str(int(round(-lon0))) + r'$^{\circ}$W'
    else:
        lon_string = str(int(round(lon0))) + r'$^{\circ}$E'
    # Make sure we are in the range 0-360
    if lon0 < 0:
        lon0 += 360

    # Set up figure
    fig = figure(figsize=(18,12))
    gs = GridSpec(2,2)
    gs.update(left=0.13, right=0.9, bottom=0.05, top=0.9, wspace=0.05, hspace=0.28)

    # Read 3D temperature, salinity, and grid variables at the beginning
    id = Dataset(file_beg, 'r')
    temp_3d_beg = id.variables['temp'][tstep_beg-1,:,:,:]
    salt_3d_beg = id.variables['salt'][tstep_beg-1,:,:,:]
    zeta_beg = id.variables['zeta'][tstep_beg-1,:,:]
    h = id.variables['h'][:,:]
    zice = id.variables['zice'][:,:]
    lon_2d = id.variables['lon_rho'][:,:]
    lat_2d = id.variables['lat_rho'][:,:]
    # Read time axis and convert to Date objects
    time_id = id.variables['ocean_time']
    time_beg = num2date(time_id[tstep_beg-1], units=time_id.units, calendar=time_id.calendar.lower())
    id.close()
    # Get a 3D array of z-coordinates
    z_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N, zeta_beg)
    # Get the date for the title
    date_string_beg = str(time_beg.day) + ' ' + month_names[time_beg.month-1] + ' ' + str(time_beg.year)
    # Interpolate to lon0
    temp_beg, z, lat = interp_lon_roms(temp_3d_beg, z_3d, lat_2d, lon_2d, lon0)
    salt_beg, z, lat = interp_lon_roms(salt_3d_beg, z_3d, lat_2d, lon_2d, lon0)
    # Plot temperature
    ax = subplot(gs[0,0])
    img = contourf(lat, z, temp_beg, lev1, extend='both', cmap='jet')
    title(r'Temperature ($^{\circ}$C)', fontsize=24)
    ax.set_xlim([lat_min, lat_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels([])
    ax.set_ylim([depth_min, 0])
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels, fontsize=16)
    ylabel('Depth (m)', fontsize=18)
    # Plot salinity
    ax = subplot(gs[0,1])
    img = contourf(lat, z, salt_beg, lev2, extend='both', cmap='jet')
    title('Salinity (psu)', fontsize=24)
    ax.set_xlim([lat_min, lat_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels([])
    ax.set_ylim([depth_min, 0])
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels([])
    # Label the timestep
    text(0.5, 0.95, date_string_beg + ', ' + lon_string, ha='center', transform=fig.transFigure, fontsize=28)

    # Read 3D temperature, salinity, and sea surface height at the end
    id = Dataset(file_end, 'r')
    temp_3d_end = id.variables['temp'][tstep_end-1,:,:,:]
    salt_3d_end = id.variables['salt'][tstep_end-1,:,:,:]
    zeta_end = id.variables['zeta'][tstep_end-1,:,:]
    # Read time axis and convert to Date objects
    time_id = id.variables['ocean_time']
    time_end = num2date(time_id[tstep_end-1], units=time_id.units, calendar=time_id.calendar.lower())
    id.close()
    # Get a 3D array of z-coordinates
    z_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N, zeta_end)
    # Get the date for the title
    date_string_end = str(time_end.day) + ' ' + month_names[time_end.month-1] + ' ' + str(time_end.year)
    # Interpolate to lon0
    temp_end, z, lat = interp_lon_roms(temp_3d_end, z_3d, lat_2d, lon_2d, lon0)
    salt_end, z, lat = interp_lon_roms(salt_3d_end, z_3d, lat_2d, lon_2d, lon0)
    # Plot temperature
    ax = subplot(gs[1,0])
    img = contourf(lat, z, temp_end, lev1, extend='both', cmap='jet')
    title(r'Temperature $^{\circ}$C)', fontsize=24)
    ax.set_xlim([lat_min, lat_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    xlabel('Latitude', fontsize=18)
    ax.set_ylim([depth_min, 0])
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels, fontsize=16)
    ylabel('Depth (m)', fontsize=18)
    # Colourbar
    cbaxes = fig.add_axes([0.03, 0.3, 0.02, 0.4])
    cbar = colorbar(img, cax=cbaxes, ticks=arange(var_min[0], var_max[0]+var_tick[0], var_tick[0]))
    cbar.ax.tick_params(labelsize=16)
    # Plot salinity
    ax = subplot(gs[1,1])
    img = contourf(lat, z, salt_end, lev2, extend='both', cmap='jet')
    title('Salinity (psu)', fontsize=24)
    ax.set_xlim([lat_min, lat_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    xlabel('Latitude', fontsize=18)
    ax.set_ylim([depth_min, 0])
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels([])
    # Colourbar
    cbaxes = fig.add_axes([0.93, 0.3, 0.02, 0.4])
    cbar = colorbar(img, cax=cbaxes, ticks=arange(var_min[1], var_max[1]+var_tick[1], var_tick[1]))
    cbar.ax.tick_params(labelsize=16)
    # Label the timestep
    text(0.5, 0.47, date_string_end + ', ' + lon_string, ha='center', transform=fig.transFigure, fontsize=28)
    fig.show()
    fig.savefig('bugs_convection_slices.png')


# Command-line interface
if __name__ == "__main__":

    bugs_convection_slices()
    
    
    
