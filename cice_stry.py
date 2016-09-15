from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from rotate_vector_cice import *

def cice_stry (file_path, tstep, colour_bounds=None, save=False, fig_name=None):

    deg2rad = pi/180

    id = Dataset(file_path, 'r')
    lon = id.variables['ULON'][:-15,:]
    lat = id.variables['ULAT'][:-15,:]
    angle = id.variables['ANGLET'][:-15,:]
    strx = id.variables['strairx'][tstep-1,:-15,:] + id.variables['strocnx'][tstep-1,:-15,:] + id.variables['strtltx'][tstep-1,:-15,:] + id.variables['strcorx'][tstep-1,:-15,:] + id.variables['strintx'][tstep-1,:-15,:]
    stry = id.variables['strairy'][tstep-1,:-15,:] + id.variables['strocny'][tstep-1,:-15,:] + id.variables['strtlty'][tstep-1,:-15,:] + id.variables['strcory'][tstep-1,:-15,:] + id.variables['strinty'][tstep-1,:-15,:]
    id.close()

    strx_lonlat, stry_lonlat = rotate_vector_cice(strx, stry, angle)

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
        max_val = amax(abs(stry_lonlat))
        lev = linspace(-max_val, max_val, num=40)
        colour_map = 'RdYlBu_r'

    fig = figure(figsize=(16,12))
    fig.add_subplot(1,1,1, aspect='equal')
    contourf(x, y, stry_lonlat, lev, cmap=colour_map, extend='both')
    cbar = colorbar()
    cbar.ax.tick_params(labelsize=20)
    title('Sum of northward stresses N/m^2', fontsize=30)
    axis('off')

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


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

    cice_stry(file_path, tstep, colour_bounds, save, fig_name)

    
