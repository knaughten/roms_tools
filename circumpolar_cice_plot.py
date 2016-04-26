from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *

# Make a circumpolar Antarctic plot of the given (horizontal) variable from CICE.
# Input:
# file_path = path to CICE history file
# var_name = name of variable in file_path to plot
# tstep = timestep in file_path to plot (1-indexed)
# save = optional boolean flag indicating that the plot should be saved to a file rather than
#        displayed on the screen
# fig_name = if save=True, filename for figure
def circumpolar_cice_plot (file_path, var_name, tstep, save=False, fig_name=None):

    deg2rad = pi/180

    # Read the variable and figure out which grid it's on
    id = Dataset(file_path, 'r')
    data = id.variables[var_name][tstep-1,:,:]
    if var_name == 'aice':
        units = 'fraction'
    else:
        units = id.variables[var_name].units
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
    lon = id.variables[lon_name][:,:]
    lat = id.variables[lat_name][:,:]
    id.close()

    # Convert to spherical coordinates
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)

    # Center levels on 0 for certain variables, with a blue-to-red colourmap
    if var_name in ['uvel', 'vvel', 'uatm', 'vatm', 'uocn', 'vocn']:
        max_val = amax(abs(data))
        lev = linspace(-max_val, max_val, num=40)
        colour_map = 'RdYlBu_r'
    else:
        lev = linspace(amin(data), amax(data), num=40)
        colour_map = 'jet'

    # Plot
    fig = figure(figsize=(16,12))
    fig.add_subplot(1,1,1, aspect='equal')
    contourf(x, y, data, lev, cmap=colour_map)
    cbar = colorbar()
    cbar.ax.tick_params(labelsize=20)
    title(var_name+' ('+units+')', fontsize=30)
    axis('off')

    if save:
        savefig(fig_name)
        close()
    else:
        show()


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to CICE history file: ")
    var_name = raw_input("Variable name: ")
    tstep = int(raw_input("Timestep number (starting at 1): "))
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    circumpolar_cice_plot(file_path, var_name, tstep, save, fig_name)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the input parameters; repeat until the user is finished
                changes = raw_input("Enter a parameter to change: (1) file path, (2) variable name, (3) timestep number, (4) save/display; or enter to continue: ")
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
                        # Change from display to save, or vice versa
                        save = not save
            if save:
                # Get file name for figure
                fig_name = raw_input("File name for figure: ")

            # Make the plot
            circumpolar_cice_plot(file_path, var_name, tstep, save, fig_name)

        else:
            break
                
    
        
        
        
        

    
    
