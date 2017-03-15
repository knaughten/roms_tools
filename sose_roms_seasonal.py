from netCDF4 import Dataset, num2date
from numpy import *
from matplotlib.pyplot import *
from calc_z import *
from seasonal_avg_roms import *
from interp_lon_roms import *
from interp_lon_sose import *

# Make a 4x2 plot comparing lat vs. depth slices of seasonally averaged 
# temperature or salinity at the given longitude, between ROMS (last year of 
# simulation) and SOSE (2005-2010 climatology).
# Input:
# file_path = path to ROMS output file, containing at least one complete 
#             Dec-Nov period (if there are multiple such periods, the last one
#             will be used for seasonal averages)
# lon0 = the specific longitude to plot (between -180 and 180)
# depth_bdry = deepest depth to plot (negative, in m)
# var_name = 'temp' for temperature or 'salt' for salinity
# save = optional boolean flag; if True, the figure will be saved with file name
#        fig_name; if False, the figure will display on the screen
# fig_name = optional string containing filename for figure, if save=True
def sose_roms_seasonal (file_path, var_name, lon0, depth_bdry, save=False, fig_name=None):

    # Path to SOSE seasonal climatology file
    sose_file = '../SOSE_seasonal_climatology.nc'

    # Grid parameters
    theta_s = 4.0
    theta_b = 0.9
    hc = 40
    N = 31

    # Season names for titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']

    # Bounds on colour scale
    if var_name == 'temp':
        var_min = -2.5
        var_max = 7.5
        var_ticks = 1
    elif var_name == 'salt':
        var_min = 33.6 #33.8 #33.6
        var_max = 35.0 #34.8 #35.0
        var_ticks = 0.4 #0.2 #0.4
    else:
        print 'Unknown variable ' + var_name
        return

    # Choose what to write on the title about the variable
    if var_name == 'temp':
        var_string = r'Temperature ($^{\circ}$C)'
    elif var_name == 'salt':
        var_string = 'Salinity (psu)'
    # Choose what to write on the title about longitude
    if lon0 < 0:
        lon_string = ' at ' + str(int(round(-lon0))) + r'$^{\circ}$W'
    else:
        lon_string = ' at ' + str(int(round(lon0))) + r'$^{\circ}$E'
    # Edit longitude bounds to be from 0 to 360, to fit with ROMS convention
    if lon0 < 0:
        lon0 += 360

    print 'Processing ROMS data'

    # Read grid
    id = Dataset(file_path, 'r')
    h = id.variables['h'][:-15,:]
    zice = id.variables['zice'][:-15,:]
    lon_roms_2d = id.variables['lon_rho'][:-15,:]
    lat_roms_2d = id.variables['lat_rho'][:-15,:]
    num_lon = id.variables['lon_rho'].shape[1]
    num_lat = id.variables['lat_rho'].shape[0]
    id.close()

    # Calculate seasonal averages of ROMS data
    var_3d_roms = seasonal_avg_roms(file_path, var_name, [N, num_lat, num_lon])
    # Chop off northern boundary
    var_3d_roms = var_3d_roms[:,:,:-15,:]

    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    z_roms_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N)

    # Calculate zonal slices for each season
    var_roms = ma.empty([4, N, size(lat_roms_2d,0)])
    var_roms[:,:,:] = 0.0
    for season in range(4):
        print 'Calculating zonal slices for ' + season_names[season]
        var_tmp, z_roms, lat_roms = interp_lon_roms(var_3d_roms[season,:,:,:], z_roms_3d, lat_roms_2d, lon_roms_2d, lon0)
        var_roms[season,:,:] = var_tmp

    print 'Processing SOSE data'

    # Read grid and 3D data (already seasonally averaged)
    id = Dataset(sose_file, 'r')
    lon_sose = id.variables['longitude'][0,:]
    lat_sose = id.variables['latitude'][:,0]
    z_sose = id.variables['depth'][:]
    var_3d_sose = id.variables[var_name][:,:,:,:]

    # Calculate zonal slices for each season
    var_sose = ma.empty([4, size(z_sose), size(lat_sose,0)])
    var_sose[:,:,:] = 0.0
    for season in range(4):
        print 'Calculating zonal slices for ' + season_names[season]
        var_sose[season,:,:] = interp_lon_sose(var_3d_sose[season,:,:,:], lon_sose, lon0)

    # Set colour levels
    lev = linspace(var_min, var_max, num=50)

    # Choose southern boundary based on extent of SOSE grid
    sbdry = amin(lat_sose)
    # Choose northern boundary based on extent of ROMS grid
    nbdry = amax(lat_roms)

    # Plot
    print 'Plotting'
    fig = figure(figsize=(20,9))
    # Loop over seasons
    for season in range(4):
        # ROMS
        fig.add_subplot(2, 4, season+1)
        img = contourf(lat_roms, z_roms, var_roms[season,:,:], lev, cmap='jet', extend='both')
        xlim([sbdry, nbdry])
        ylim([depth_bdry, 0])
        title('ROMS (' + season_names[season] + ')', fontsize=24)
        if season == 0:
            ylabel('depth (m)', fontsize=18)
        # SOSE
        fig.add_subplot(2, 4, season+5)
        contourf(lat_sose, z_sose, var_sose[season,:,:], lev, cmap='jet', extend='both')
        xlim([sbdry, nbdry])
        ylim([depth_bdry, 0])
        title('SOSE (' + season_names[season] + ')', fontsize=24)
        xlabel('Latitude', fontsize=18)
        if season == 0:
            ylabel('depth (m)', fontsize=18)
    # Add colourbar
    cbaxes = fig.add_axes([0.93, 0.2, 0.015, 0.6])
    cbar = colorbar(img, cax=cbaxes, ticks=arange(var_min, var_max+var_ticks, var_ticks))
    cbar.ax.tick_params(labelsize=16)
    # Add the main title
    suptitle(var_string + lon_string, fontsize=30)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()    


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to ocean averages file, containing at least one complete Dec-Nov period: ")
    var_key = raw_input("Temperature (t) or salinity (s)? ")
    if var_key == 't':
        var_name = 'temp'
    elif var_key == 's':
        var_name = 'salt'
    lon0 = float(raw_input("Enter longitude (-180 to 180): "))
    depth_bdry = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    sose_roms_seasonal(file_path, var_name, lon0, depth_bdry, save, fig_name)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the input parameters; repeat until the user is finished
                changes = raw_input("Enter a parameter to change: (1) file path, (2) temperature/salinity, (3) longitude, (4) deepest depth, (5) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New file path
                        file_path = raw_input("Path to ocean averages file, containing at least one complete Dec-Nov period: ")
                    elif int(changes) == 2:
                        # Switch from temperature to salinity or vice versa
                        if var_name == 'temp':
                            var_name = 'salt'
                        else:
                            var_name = 'temp'
                    elif int(changes) == 3:
                        # New longitude
                        lon0 = float(raw_input("Enter longitude (-180 to 180): "))
                    elif int(changes) == 4:
                        # New depth bound
                        depth_bdry = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
                    elif int(changes) == 5:
                        # Change from save to display, or vice versa
                        save = not save
            if save:
                # Get file name for figure
                fig_name = raw_input("File name for figure: ")

            # Make the plot
            sose_roms_seasonal(file_path, var_name, lon0, depth_bdry, save, fig_name)
        else:
            break
        
                

    
