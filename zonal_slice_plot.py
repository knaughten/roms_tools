from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from calc_z import *

# Create a zonal slice plot (i.e. depth vs latitude, single i-index) of the
# given variable.
# Input:
# grid_path = path to ROMS grid file
# file_path = path to ocean history or output file
# var_name = variable name in file_path; must be a variable of dimension
#            time x depth x latitude x longitude, and on the rho, u, or v grid
# tstep = timestep in file_path to plot (1-indexed)
# i_val = i-index of zonal slice to plot (1-indexed)
# depth_min = deepest depth to plot (negative, metres)
# save = optional boolean flag; if True, the figure will be saved with file name
#        fig_name, if False, the figure will display on the screen
# fig_name = optional string containing filename for figure, if save=True
def zonal_slice_plot (grid_path, file_path, var_name, tstep, i_val, depth_min, save=False, fig_name=None):

    # Grid parameters
    theta_s = 0.9
    theta_b = 4.0
    hc = 40
    N = 31

    # Read the variable and figure out which grid it's on
    id = Dataset(file_path, 'r')
    data = id.variables[var_name][tstep-1,:,:,:]
    # Also read sea surface height
    zeta = id.variables['zeta'][tstep-1,:,:]
    if var_name == 'salt':
        units = 'psu'
    else:
        units = id.variables[var_name].units
    grid_string = id.variables[var_name].coordinates
    if grid_string.startswith('lon_rho'):
        grid_name = 'rho'
        lon_name = 'lon_rho'
        lat_name = 'lat_rho'
    elif grid_string.startswith('lon_u'):
        grid_name = 'u'
        lon_name = 'lon_u'
        lat_name = 'lat_u'
    elif grid_string.startswith('lon_v'):
        grid_name = 'v'
        lon_name = 'lon_v'
        lat_name = 'lat_v'
    else:
        print 'Grid type ' + grid_string + ' not supported'
        id.close()
        return
    id.close()

    # Read grid variables
    id = Dataset(grid_path, 'r')
    h = id.variables['h'][:,:]
    zice = id.variables['zice'][:,:]
    # h, zice, and zeta are on the rho-grid; interpolate if necessary
    if grid_name == 'u':
        h = 0.5*(h[:,0:-1] + h[:,1:])
        zice = 0.5*(zice[:,0:-1] + zice[:,1:])
        zeta = 0.5*(zeta[:,0:-1] + zeta[:,1:])
    elif grid_name == 'v':
        h = 0.5*(h[0:-1,:] + h[1:,:])
        zice = 0.5*(zice[0:-1,:] + zice[1:,:])
        zeta = 0.5*(zeta[0:-1,:] + zeta[1:,:])
    # Read the correct lat and lon for this grid (determined previously)
    lon = id.variables[lon_name][:,:]
    lat = id.variables[lat_name][:,:]
    id.close()

    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    z, sc_r, Cs_r = calc_z(h, zice, lon, lat, theta_s, theta_b, hc, N, zeta)

    # Select the correct slice of data, lat, and z
    data_slice = data[:,:,i_val-1]
    z_slice = z[:,:,i_val-1]
    lat_1D = lat[:,i_val-1]
    # Copy latitude into a 2D grid (depth x lat)
    lat_slice = tile(lat_1D, (N,1))
    # Select the mean longitude of this slice
    lon_mean = mean(lon[:,i_val-1])
    if lon_mean > 180:
        lon_mean -= 360
    elif lon_mean < -180:
        lon_mean += 360

    # Center levels on 0 for u and v, with a blue-to-red colourmap
    if var_name in ['u', 'v']:
        max_val = amax(abs(data_slice))
        lev = linspace(-max_val, max_val, num=40)
        colour_map = 'RdYlBu_r'
    else:
        lev = linspace(amin(data_slice), amax(data_slice), num=40)
        colour_map = 'jet'

    # Plot
    figure(figsize=(18,6))
    contourf(lat_slice, z_slice, data_slice, lev, cmap=colour_map)
    colorbar()

    # Configure plot
    if lon_mean < 0:
        lon_string = str(int(round(-lon_mean))) + r'$^{\circ}$W'
    else:
        lon_string = str(int(round(lon_mean))) + r'$^{\circ}$E'
    title(var_name + ' (' + units + ') at ' + lon_string)
    xlabel('Latitude')
    ylabel('Depth (m)')

    if save:
        savefig(fig_name)
    else:
        show()


# Command-line interface
if __name__ == "__main__":

    grid_path = raw_input("Path to ROMS grid file: ")
    file_path = raw_input("Path to ocean history/averages file: ")
    var_name = raw_input("Variable name: ")
    tstep = int(raw_input("Timestep number (starting at 1): "))
    i_val = int(raw_input("i-index of zonal slice (starting at 1): "))
    depth_min = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    zonal_slice_plot(grid_path, file_path, var_name, tstep, i_val, depth_min, save, fig_name)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the input parameters; repeat until the user is finished
                changes = raw_input("Enter a parameter to change: (1) grid path, (2) file path, (3) variable name, (4) timestep number, (5) i-index, (6) deepest depth, (7) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New grid path
                        grid_path = raw_input("Path to ROMS grid file: ")
                    elif int(changes) == 2:
                        # New file path
                        file_path = raw_input("Path to ocean history/averages file: ")
                    elif int(changes) == 3:
                        # New variable name
                        var_name = raw_input("Variable name: ")
                    elif int(changes) == 4:
                        # New timestep
                        tstep = int(raw_input("Timestep number (starting at 1): "))
                    elif int(changes) == 5:
                        # New i-index
                        i_val = int(raw_input("i-index of zonal slice (starting at 1): "))
                    elif int(changes) == 6:
                        # New depth bound
                        depth_min = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
                    elif int(changes) == 7:
                        # Change from display to save, or vice versa
                        save = not save
            if save:
                # Get file name for figure
                fig_name = raw_input("File name for figure: ")

            # Make the plot
            zonal_slice_plot(grid_path, file_path, var_name, tstep, i_val, depth_min, save, fig_name)

        else:
            break
                        

    

    

    
