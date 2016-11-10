from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from rotate_vector_cice import *

# Make a circumpolar plot of effective sea ice thickness (area*height).
# Input:
# file_path = path to CICE history file
# tstep = timestep in file_path to plot (1-indexed)
# colour_bounds = optional bounds on colour scale, stored as an array of size
#                 2 with the lower bound first. If colour_bounds = None, then
#                 determine colour scale bounds automatically.
# save = optional boolean flag indicating that the plot should be saved to a file rather than
#        displayed on the screen
# fig_name = if save=True, filename for figure
def effective_thickness (file_path, tstep, colour_bounds=None, save=False, fig_name=None):

    deg2rad = pi/180

    # Read sea ice concentration and thickness
    id = Dataset(file_path, 'r')
    aice = id.variables['aice'][tstep-1,:-15,:]
    hi = id.variables['hi'][tstep-1,:-15,:]
    data_tmp = aice*hi

    # Read the correct lat and lon for this grid
    lon_tmp = id.variables['TLON'][:-15,:]
    lat_tmp = id.variables['TLAT'][:-15,:]
    id.close()

    # Wrap the periodic boundary by 1 cell
    lon = ma.empty([size(lon_tmp,0), size(lon_tmp,1)+1])
    lat = ma.empty([size(lat_tmp,0), size(lat_tmp,1)+1])
    data = ma.empty([size(data_tmp,0), size(data_tmp,1)+1])
    lon[:,:-1] = lon_tmp
    lon[:,-1] = lon_tmp[:,0]
    lat[:,:-1] = lat_tmp
    lat[:,-1] = lat_tmp[:,0]
    data[:,:-1] = data_tmp
    data[:,-1] = data_tmp[:,0]

    # Convert to spherical coordinates
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)

    if colour_bounds is not None:
        # User has set bounds on colour scale
        lev = linspace(colour_bounds[0], colour_bounds[1], num=40)
    else:
        lev = linspace(amin(data), amax(data), num=40)

    # Plot
    fig = figure(figsize=(16,12))
    fig.add_subplot(1,1,1, aspect='equal')
    contourf(x, y, data, lev, extend='both')
    cbar = colorbar()
    cbar.ax.tick_params(labelsize=20)
    title('Effective thickness (m)', fontsize=30)
    axis('off')

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to CICE history file: ")
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
    effective_thickness(file_path, tstep, colour_bounds, save, fig_name)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the input parameters; repeat until the user is finished
                changes = raw_input("Enter a parameter to change: (1) file path, (2) timestep number, (3) colour bounds, (4) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New file path
                        file_path = raw_input("Path to CICE history file: ")
                    elif int(changes) == 2:
                        # New timestep number
                        tstep = int(raw_input("Timestep number (starting at 1): "))
                    elif int(changes) == 3:
                        # Get colour bounds if necessary
                        colour_bounds = None
                        get_bounds = raw_input("Set bounds on colour scale (y/n)? ")
                        if get_bounds == 'y':
                            lower_bound = float(raw_input("Lower bound: "))
                            upper_bound = float(raw_input("Upper bound: "))
                            colour_bounds = [lower_bound, upper_bound]
                    elif int(changes) == 4:
                        # Change from display to save, or vice versa
                        save = not save
            if save:
                # Get file name for figure
                fig_name = raw_input("File name for figure: ")

            # Make the plot
            effective_thickness(file_path, tstep, colour_bounds, save, fig_name)

        else:
            break
                
    
        
        
        
        

    
    
