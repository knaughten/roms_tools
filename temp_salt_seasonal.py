from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from calc_z import *
from seasonal_avg_roms import *
from interp_lon_roms import *

# Make a 4x2 plot showing lat vs. depth slices of seasonally averaged 
# temperature (top) and salinity (bottom) at the given longitude, over the
# last year of simulation.
# Input:
# file_path = path to ROMS output file, containing at least one complete 
#             Dec-Nov period (if there are multiple such periods, the last one
#             will be used for seasonal averages)
# lon0 = the specific longitude to plot (between -180 and 180)
# depth_bdry = deepest depth to plot (negative, in m)
# save = optional boolean flag; if True, the figure will be saved with file name
#        fig_name; if False, the figure will display on the screen
# fig_name = optional string containing filename for figure, if save=True
def temp_salt_seasonal (file_path, lon0, depth_bdry, save=False, fig_name=None):

    # Grid parameters
    theta_s = 4.0
    theta_b = 0.9
    hc = 40
    N = 31

    # Season names for titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']

    # Bounds on colour scale
    temp_min = -2.5
    temp_max = 7.5
    temp_ticks = 2
    salt_min = 33.8
    salt_max = 34.8
    salt_ticks = 0.2

    # Choose what to write on the title about longitude
    if lon0 < 0:
        lon_string = 'T/S slices at ' + str(int(round(-lon0))) + r'$^{\circ}$W'
    else:
        lon_string = 'T/S slices at ' + str(int(round(lon0))) + r'$^{\circ}$E'
    # Edit longitude bounds to be from 0 to 360, to fit with ROMS convention
    if lon0 < 0:
        lon0 += 360

    # Read grid
    id = Dataset(file_path, 'r')
    h = id.variables['h'][:-15,:]
    zice = id.variables['zice'][:-15,:]
    lon_roms_2d = id.variables['lon_rho'][:-15,:]
    lat_roms_2d = id.variables['lat_rho'][:-15,:]
    num_lon = id.variables['lon_rho'].shape[1]
    num_lat = id.variables['lat_rho'].shape[0]
    id.close()

    # Calculate seasonal averages of temperature and salinity
    print 'Reading temperature'
    temp_3d_roms = seasonal_avg_roms(file_path, 'temp', [N, num_lat, num_lon])
    print 'Reading salinity'
    salt_3d_roms = seasonal_avg_roms(file_path, 'salt', [N, num_lat, num_lon])
    # Chop off northern boundary
    temp_3d_roms = temp_3d_roms[:,:,:-15,:]
    salt_3d_roms = salt_3d_roms[:,:,:-15,:]

    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    z_roms_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N)

    # Calculate zonal slices for each season
    temp_roms = ma.empty([4, N, size(lat_roms_2d,0)])
    temp_roms[:,:,:] = 0.0
    salt_roms = ma.empty([4, N, size(lat_roms_2d,0)])
    salt_roms[:,:,:] = 0.0
    for season in range(4):
        print 'Calculating zonal slices for ' + season_names[season]
        temp_tmp, z_roms, lat_roms = interp_lon_roms(temp_3d_roms[season,:,:,:], z_roms_3d, lat_roms_2d, lon_roms_2d, lon0)
        temp_roms[season,:,:] = temp_tmp
        salt_tmp, z_roms, lat_roms = interp_lon_roms(salt_3d_roms[season,:,:,:], z_roms_3d, lat_roms_2d, lon_roms_2d, lon0)
        salt_roms[season,:,:] = salt_tmp

    # Set colour levels
    lev1 = linspace(temp_min, temp_max, num=50)
    lev2 = linspace(salt_min, salt_max, num=50)

    # Choose boundaries based on extent of ROMS grid
    sbdry = amin(lat_roms)
    nbdry = amax(lat_roms)

    # Plot
    print 'Plotting'
    fig = figure(figsize=(20,9))
    # Loop over seasons
    for season in range(4):
        # Temperature
        fig.add_subplot(2, 4, season+1)
        img = contourf(lat_roms, z_roms, temp_roms[season,:,:], lev1, cmap='jet', extend='both')
        xlim([sbdry, nbdry])
        ylim([depth_bdry, 0])
        title('Temperature (' + season_names[season] + ')', fontsize=24)
        if season == 0:
            ylabel('depth (m)', fontsize=18)
        if season == 3:
            cbaxes1 = fig.add_axes([0.92, 0.55, 0.01, 0.3])
            cbar1 = colorbar(img, cax=cbaxes1, ticks=arange(temp_min, temp_max+temp_ticks, temp_ticks))
            cbar1.ax.tick_params(labelsize=16)
        # Salinity
        fig.add_subplot(2, 4, season+5)
        img = contourf(lat_roms, z_roms, salt_roms[season,:,:], lev2, cmap='jet', extend='both')
        xlim([sbdry, nbdry])
        ylim([depth_bdry, 0])
        title('Salinity (' + season_names[season] + ')', fontsize=24)
        if season == 0:
            ylabel('depth (m)', fontsize=18)
        xlabel('Latitude', fontsize=18)
        if season == 3:
            cbaxes2 = fig.add_axes([0.92, 0.15, 0.01, 0.3])
            cbar2 = colorbar(img, cax=cbaxes2, ticks=arange(salt_min, salt_max+salt_ticks, salt_ticks))
            cbar2.ax.tick_params(labelsize=16)
    # Add the main title
    suptitle(lon_string, fontsize=30)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()    


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to ocean averages file, containing at least one complete Dec-Nov period: ")
    lon0 = float(raw_input("Enter longitude (-180 to 180): "))
    depth_bdry = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    temp_salt_seasonal(file_path, lon0, depth_bdry, save, fig_name)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the input parameters; repeat until the user is finished
                changes = raw_input("Enter a parameter to change: (1) file path, (2) longitude, (3) deepest depth, (4) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New file path
                        file_path = raw_input("Path to ocean averages file, containing at least one complete Dec-Nov period: ")
                    elif int(changes) == 2:
                        # New longitude
                        lon0 = float(raw_input("Enter longitude (-180 to 180): "))
                    elif int(changes) == 3:
                        # New depth bound
                        depth_bdry = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
                    elif int(changes) == 4:
                        # Change from save to display, or vice versa
                        save = not save
            if save:
                # Get file name for figure
                fig_name = raw_input("File name for figure: ")

            # Make the plot
            temp_salt_seasonal(file_path, lon0, depth_bdry, save, fig_name)
        else:
            break
        
                

    
