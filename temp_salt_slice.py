from netCDF4 import Dataset, num2date
from numpy import *
from matplotlib.pyplot import *
from calc_z import *
from interp_lon_roms import *

# Create a 2x1 plot showing zonal slices (depth vs latitude) of temperature and
# salinity interpolated to the given longitude, at the given timestep.
# Input:
# file_path = path to ROMS history or averages file
# tstep = time index in file_path to plot (1-indexed)
# lon0 = longitude to interpolate to (-180 to 180)
# depth_min = deepest depth to plot (negative, in metres)
# save = optional boolean flag; if True, the figure will be saved with file name
#        fig_name, if False, the figure will display on the screen
# fig_name = optional string containing filename for figure, if save=True
def temp_salt_slice (file_path, tstep, lon0, depth_min, save=False, fig_name=None):

    # Grid parameters
    theta_s = 7.0
    theta_b = 2.0
    hc = 250
    N = 31

    # Month names for titles
    month_names = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']

    # Bounds on colour scales for temperature and salinity
    var_min = [-2, 33.8]
    var_max = [3, 34.8]
    var_tick = [1, 0.2]

    # Read temperature, salinity, and grid variables
    id = Dataset(file_path, 'r')
    temp_3d = id.variables['temp'][tstep-1,:,:-15,:]
    salt_3d = id.variables['salt'][tstep-1,:,:-15,:]
    zeta = id.variables['zeta'][tstep-1,:-15,:]
    h = id.variables['h'][:-15,:]
    zice = id.variables['zice'][:-15,:]
    lon_2d = id.variables['lon_rho'][:-15,:]
    lat_2d = id.variables['lat_rho'][:-15,:]
    # Read time axis and convert to Date objects
    time_id = id.variables['ocean_time']
    time = num2date(time_id[tstep-1], units=time_id.units, calendar=time_id.calendar.lower())
    id.close()

    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    z_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N, zeta)

    # Get the date for the title
    date_string = str(time.day) + ' ' + month_names[time.month-1] + ' ' + str(time.year)

    # Get longitude for the title
    if lon0 < 0:
        lon_string = str(int(round(-lon0))) + r'$^{\circ}$W'
    else:
        lon_string = str(int(round(lon0))) + r'$^{\circ}$E'

    # Make sure we are in the range 0-360
    if lon0 < 0:
        lon0 += 360

    # Interpolate temperature, salinity, z, and latitude to lon0
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

    # Colour levels
    lev1 = linspace(var_min[0], var_max[0], num=50)
    lev2 = linspace(var_min[1], var_max[1], num=50)

    # Plot
    fig = figure(figsize=(24,6))
    # Temperature
    ax = fig.add_subplot(1,2,1)
    img1 = contourf(lat, z, temp, lev1, extend='both')
    xlim([lat_min, lat_max])
    ylim([depth_min, 0])
    xlabel('Latitude')
    ylabel('Depth (m)')
    title(r'Temperature ($^{\circ}$C)', fontsize=20)
    cbar1 = colorbar(img1, ticks=arange(var_min[0], var_max[0]+var_tick[0], var_tick[0]))
    cbar1.ax.tick_params(labelsize=16)
    # Salinity
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

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()

    # Convert back to the range -180 to 180, in case this script repeats
    if lon0 > 180:
        lon0 -= 360


# Command-line interface
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
    # Make the plot
    temp_salt_slice(file_path, tstep, lon0, depth_min, save, fig_name)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            # Ask for changes to parameters until the user is done
            while True:
                changes = raw_input("Enter a parameter to change: (1) file path, (2) time index, (3) longitude, (4) deepest depth, (5) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
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
                # Get new figure name
                fig_name = raw_input("File name for figure: ")
            # Make the plot
            temp_salt_slice(file_path, tstep, lon0, depth_min, save, fig_name)
        else:
            # No more figures
            break
                

    
    
