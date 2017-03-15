from netCDF4 import Dataset, num2date
from numpy import *
from matplotlib.pyplot import *
from calc_z import *
from interp_lon_roms import *

def temp_salt_slice (file_path, tstep, lon0, depth_min, save=False, fig_name=None):

    # Grid parameters
    theta_s = 4.0
    theta_b = 0.9
    hc = 40
    N = 31

    month_names = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']

    var_min = [-2, 33.8]
    var_max = [3, 34.8]
    var_tick = [1, 0.2]

    id = Dataset(file_path, 'r')
    temp_3d = id.variables['temp'][tstep-1,:,:-15,:]
    salt_3d = id.variables['salt'][tstep-1,:,:-15,:]
    zeta = id.variables['zeta'][tstep-1,:-15,:]
    h = id.variables['h'][:-15,:]
    zice = id.variables['zice'][:-15,:]
    lon_2d = id.variables['lon_rho'][:-15,:]
    lat_2d = id.variables['lat_rho'][:-15,:]
    time_id = id.variables['ocean_time']
    time = num2date(time_id[tstep-1], units=time_id.units, calendar=time_id.calendar.lower())
    id.close()

    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    z_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N, zeta)

    date_string = str(time.day) + ' ' + month_names[time.month-1] + ' ' + str(time.year)

    if lon0 < 0:
        lon_string = str(int(round(-lon0))) + r'$^{\circ}$W'
    else:
        lon_string = str(int(round(lon0))) + r'$^{\circ}$E'

    if lon0 < 0:
        lon0 += 360

    temp, z, lat = interp_lon_roms(temp_3d, z_3d, lat_2d, lon_2d, lon0)
    salt, z, lat = interp_lon_roms(salt_3d, z_3d, lat_2d, lon_2d, lon0)

    # Choose latitude bounds based on land mask
    temp_sum = sum(temp, axis=0)    
    # Find southernmost and northernmost unmasked j-indices
    edges = ma.flatnotmasked_edges(temp_sum)
    j_min = edges[0]
    j_max = edges[1]
    if j_min == 0:
        # There are ocean points right to the southern boundary
        # Don't do anything special
        lat_min = min(lat[:,j_min])
    else:
        # There is land everywhere at the southern boundary
        # Show the last 2 degrees of this land mask
        lat_min = min(lat[:,j_min]) - 2
    if j_max == size(temp_sum) - 1:
        # There are ocean points right to the northern boundary
        # Don't do anything special
        lat_max = max(lat[:,j_max])
    else:
        # There is land everywhere at the northern boundary
        # Show the first 2 degrees of this land mask
        lat_max = max(lat[:,j_max]) + 2    

    lev1 = linspace(var_min[0], var_max[0], num=50)
    lev2 = linspace(var_min[1], var_max[1], num=50)

    fig = figure(figsize=(24,6))
    ax = fig.add_subplot(1,2,1)
    img1 = contourf(lat, z, temp, lev1, extend='both')
    xlim([lat_min, lat_max])
    ylim([depth_min, 0])
    xlabel('Latitude')
    ylabel('Depth (m)')
    title(r'Temperature ($^{\circ}$C)', fontsize=20)
    cbar1 = colorbar(img1, ticks=arange(var_min[0], var_max[0]+var_tick[0], var_tick[0]))
    cbar1.ax.tick_params(labelsize=16)
    ax = fig.add_subplot(1,2,2)
    img2 = contourf(lat, z, salt, lev2, extend='both')
    xlim([lat_min, lat_max])
    ylim([depth_min, 0])
    xlabel('Latitude')
    ylabel('Depth (m)')
    title('Salinity (psu)', fontsize=20)
    cbar2 = colorbar(img2, ticks=arange(var_min[1], var_max[1]+var_tick[1], var_tick[1]))
    cbar2.ax.tick_params(labelsize=16)
    suptitle(date_string + ', ' + lon_string, fontsize=24)
    subplots_adjust(wspace=0.025)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()

    if lon0 > 180:
        lon0 -= 360


if __name__ == "__main__":

    file_path = raw_input("Path to ocean averages file: ")
    tstep = int(raw_input("Time index to plot (starting at 1): "))
    lon0 = float(raw_input("Longitude to plot (-180 to 180): "))
    depth_min = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    temp_salt_slice(file_path, tstep, lon0, depth_min, save, fig_name)

    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                changes = raw_input("Enter a parameter to change: (1) file path, (2) time index, (3) longitude, (4) deepest depth, (5) save/display; or enter to continue: ")
                if len(changes) == 0:
                    break
                else:
                    if int(changes) == 1:
                        file_path = raw_input("Path to ocean averages file: ")
                    elif int(changes) == 2:
                        tstep = int(raw_input("Time index to plot (starting at 1): "))
                    elif int(changes) == 3:
                        lon0 = float(raw_input("Longitude to plot (-180 to 180): "))
                    elif int(changes) == 4:
                        depth_min = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
                    elif int(changes) == 5:
                        save = not save
            if save:
                fig_name = raw_input("File name for figure: ")
            temp_salt_slice(file_path, tstep, lon0, depth_min, save, fig_name)
        else:
            break
                

    
    
