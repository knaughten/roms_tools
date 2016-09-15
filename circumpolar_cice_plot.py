from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from rotate_vector_cice import *

# Make a circumpolar Antarctic plot of the given (horizontal) variable from CICE.
# Input:
# file_path = path to CICE history file
# var_name = name of variable in file_path to plot
# tstep = timestep in file_path to plot (1-indexed)
# colour_bounds = optional bounds on colour scale, stored as an array of size
#                 2 with the lower bound first. If colour_bounds = None, then
#                 determine colour scale bounds automatically.
# save = optional boolean flag indicating that the plot should be saved to a file rather than
#        displayed on the screen
# fig_name = if save=True, filename for figure
def circumpolar_cice_plot (file_path, var_name, tstep, colour_bounds=None, save=False, fig_name=None):

    deg2rad = pi/180

    # Read the variable
    id = Dataset(file_path, 'r')
    data = id.variables[var_name][tstep-1,:-15,:]
    if var_name == 'aice':
        units = 'fraction'
    else:
        units = id.variables[var_name].units

    # Check for vector variables that need to be rotated
    if var_name in ['uvel', 'vvel', 'uatm', 'vatm', 'uocn', 'vocn', 'strairx', 'strairy', 'strtltx', 'strtlty', 'strcorx', 'strcory', 'strocnx', 'strocny', 'strintx', 'strinty']:
        angle = id.variables['ANGLET'][:-15,:]
        if var_name in ['uvel', 'uatm', 'uocn', 'strairx', 'strtltx', 'strcorx', 'strocnx', 'strintx']:
            # u-variable
            u_data = data[:,:]
            if var_name[0] == 'u':
                v_data = id.variables[var_name.replace('u','v')][tstep-1,:-15,:]
            else:
                v_data = id.variables[var_name.replace('x','y')][tstep-1,:-15,:]
            u_data_lonlat, v_data_lonlat = rotate_vector_cice(u_data, v_data, angle)
            data = u_data_lonlat
        elif var_name in ['vvel', 'vatm', 'vocn', 'strairy', 'strtlty', 'strcory', 'strocny', 'strinty']:
            # v-variable
            v_data = data[:,:]
            if var_name[0] == 'v':
                u_data = id.variables[var_name.replace('v','u',1)][tstep-1,:-15,:]
            else:
                u_data = id.variables[var_name.replace('y','x')][tstep-1,:-15,:]
            u_data_lonlat, v_data_lonlat = rotate_vector_cice(u_data, v_data, angle)
            data = v_data_lonlat

    # Figure out which grid we're on
    grid_string = id.variables[var_name].coordinates
    if grid_string.startswith('ULON'):
        grid_name = 'u'
        lon_name = 'ULON'
        lat_name = 'ULAT'
    elif grid_string.startswith('TLON'):
        grid_name = 't'
        lon_name = 'TLON'
        lat_name = 'TLAT'
    else:
        print 'Grid type ' + grid_string + ' not supported'
        id.close()
        return

    # Read the correct lat and lon for this grid
    lon = id.variables[lon_name][:-15,:]
    lat = id.variables[lat_name][:-15,:]
    id.close()

    # Convert to spherical coordinates
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)

    if colour_bounds is not None:
        # User has set bounds on colour scale
        lev = linspace(colour_bounds[0], colour_bounds[1], num=40)
        if colour_bounds[0] == -colour_bounds[1]:
            # Bounds are centered on zero, so choose a blue-to-red colourmap
            # centered on yellow
            colour_map = 'RdYlBu_r'
        else:
            colour_map = 'jet'
    else:
        # Determine bounds automatically
        if var_name in ['uvel', 'vvel', 'uatm', 'vatm', 'uocn', 'vocn', 'fresh_ai', 'fsalt_ai', 'fhocn_ai', 'strairx', 'strairy', 'strtltx', 'strtlty', 'strcorx', 'strcory', 'strocnx', 'strocny', 'strintx', 'strinty']:
            # Center levels on 0 for certain variables, with a blue-to-red
            # colourmap
            max_val = amax(abs(data))
            lev = linspace(-max_val, max_val, num=40)
            colour_map = 'RdYlBu_r'
        else:
            lev = linspace(amin(data), amax(data), num=40)
            colour_map = 'jet'

    # Plot
    fig = figure(figsize=(16,12))
    fig.add_subplot(1,1,1, aspect='equal')
    contourf(x, y, data, lev, cmap=colour_map, extend='both')
    cbar = colorbar()
    cbar.ax.tick_params(labelsize=20)
    title(var_name+' ('+units+')', fontsize=30)
    axis('off')

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to CICE history file: ")
    var_name = raw_input("Variable name: ")
    tstep = int(raw_input("Timestep number (starting at 1): "))

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
    circumpolar_cice_plot(file_path, var_name, tstep, colour_bounds, save, fig_name)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the input parameters; repeat until the user is finished
                changes = raw_input("Enter a parameter to change: (1) file path, (2) variable name, (3) timestep number, (4) colour bounds, (5) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New file path
                        file_path = raw_input("Path to CICE history file: ")
                    elif int(changes) == 2:
                        # New variable name
                        var_name = raw_input("Variable name: ")
                    elif int(changes) == 3:
                        # New timestep number
                        tstep = int(raw_input("Timestep number (starting at 1): "))
                    elif int(changes) == 4:
                        # Get colour bounds if necessary
                        colour_bounds = None
                        get_bounds = raw_input("Set bounds on colour scale (y/n)? ")
                        if get_bounds == 'y':
                            lower_bound = float(raw_input("Lower bound: "))
                            upper_bound = float(raw_input("Upper bound: "))
                            colour_bounds = [lower_bound, upper_bound]
                    elif int(changes) == 5:
                        # Change from display to save, or vice versa
                        save = not save
            if save:
                # Get file name for figure
                fig_name = raw_input("File name for figure: ")

            # Make the plot
            circumpolar_cice_plot(file_path, var_name, tstep, colour_bounds, save, fig_name)

        else:
            break
                
    
        
        
        
        

    
    
