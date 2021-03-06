from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from calc_z import *

# Plot the difference from the freezing temperature for a single i-slice
# (depth versus y) with no time-averaging, spatial averaging, or interpolation.
# This helps to show where there is spurious supercooling.
# Input:
# file_path = path to ROMS history file (NOT averages)
# tstep = timestep to plot (starting at 1)
# i_val = i-index to plot (starting at 1)
# depth_min = deepest depth to plot (negative, in metres)
# colour_bounds = optional bounds on colour scale, stored as an array of size
#                 2 with the lower bound first. If colour_bounds = None, then
#                 determine colour scale bounds automatically.
# save = optional boolean flag; if True, the figure will be saved with file name
#        fig_name, if False, the figure will display on the screen
# fig_name = optional string containing filename for figure, if save=True
def freezingpt_slice (file_path, tstep, i_val, depth_min, colour_bounds=None, save=False, fig_name=None):

    # Grid parameters
    theta_s = 7.0
    theta_b = 2.0
    hc = 250
    N = 31

    # Read temperature, salinity, and grid variables
    id = Dataset(file_path, 'r')
    temp = id.variables['temp'][tstep-1,:,:-15,i_val-1]
    salt = id.variables['salt'][tstep-1,:,:-15,i_val-1]
    h = id.variables['h'][:-15,:]
    zice = id.variables['zice'][:-15,:]
    # Sea surface height is time-dependent
    zeta = id.variables['zeta'][tstep-1,:-15,:]
    lon_2d = id.variables['lon_rho'][:-15,:]
    lat_2d = id.variables['lat_rho'][:-15,:]
    id.close()

    # Calculate freezing point as seen by supercooling code
    tfr = salt/(-18.48 + salt*18.48/1000.0) #-0.054*salt
    # Calculate difference from freezing point
    deltat = temp - tfr

    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    z_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N, zeta)
    # Select depth and latitude at the given i-index
    z = z_3d[:,:,i_val-1]
    lat = tile(lat_2d[:,i_val-1], (N,1))

    # Determine colour bounds
    if colour_bounds is not None:
        # Specified by user
        scale_min = colour_bounds[0]
        scale_max = colour_bounds[1]
        if scale_min == -scale_max:
            # Centered on zero; use a red-yellow-blue colour scale
            colour_map = 'RdBu_r'
        else:
            # Use a rainbow colour scale
            colour_map = 'jet'
    else:
        # Determine automatically
        scale_min = amin(deltat)
        scale_max = amax(deltat)
        colour_map = 'jet'

    # Plot (pcolor not contour to show what each individual cell is doing)
    fig = figure(figsize=(18,6))
    pcolor(lat, z, deltat, vmin=scale_min, vmax=scale_max, cmap=colour_map)
    colorbar()
    title('Difference from freezing point (K) at i=' + str(i_val))
    xlabel('Latitude')
    ylabel('Depth (m)')

    # Determine bounds on latitude
    data_sum = sum(deltat, axis=0)
    edges = ma.flatnotmasked_edges(data_sum)
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
    if j_max == size(data_sum) - 1:
        # There are ocean points right to the northern boundary
        # Don't do anything special
        lat_max = max(lat[:,j_max])
    else:
        # There is land everywhere at the northern boundary
        # Show the first 2 degrees of this land mask
        lat_max = max(lat[:,j_max]) + 2
    lat_max = -65
    xlim([lat_min, lat_max])
    ylim([depth_min, 0])

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to ocean history file: ")
    tstep = int(raw_input("Timestep number (starting at 1): "))
    i_val = int(raw_input("i-index to plot (1 to 1443): "))
    depth_min = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
    # Get colour bounds if necessary
    colour_bounds = None
    get_bounds = raw_input("Set bounds on colour scale (y/n)? ")
    if get_bounds == 'y':
        lower_bound = float(raw_input("Lower bound: "))
        upper_bound = float(raw_input("Upper bound: "))
        colour_bounds = [lower_bound, upper_bound]

    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None

    # Make the plot
    freezingpt_slice(file_path, tstep, i_val, depth_min, colour_bounds, save, fig_name)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the input parameters; repeat until the user
                # is finished
                changes = raw_input("Enter a parameter to change: (1) file path, (2) timestep number, (3) i-index, (4) deepest depth, (5) colour bounds, (6) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New file path
                        file_path = raw_input("Path to ocean history file: ")
                    elif int(changes) == 2:
                        # New timestep
                        tstep = int(raw_input("Timestep number (starting at 1): "))
                    elif int(changes) == 3:
                        # New i-index
                        i_val = int(raw_input("i-index to plot (1 to 1443): "))
                    elif int(changes) == 4:
                        # New depth bound
                        depth_min = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
                    elif int(changes) == 5:
                        # Get colour bounds if necessary
                        colour_bounds = None
                        get_bounds = raw_input("Set bounds on colour scale (y/n)? ")
                        if get_bounds == 'y':
                            lower_bound = float(raw_input("Lower bound: "))
                            upper_bound = float(raw_input("Upper bound: "))
                            colour_bounds = [lower_bound, upper_bound]
                    elif int(changes) == 6:
                        # Change from display to save, or vice versa
                        save = not save
            if save:
                # Get file name for figure
                fig_name = raw_input("File name for figure: ")
            # Make the plot
            freezingpt_slice(file_path, tstep, i_val, depth_min, colour_bounds, save, fig_name)
        else:
            break
            
            
    

    

    
