from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *

def bw_diff_plot (file_path, var_name, colour_bounds=None, save=False, fig_name=None):

    init_path = '../ROMS-CICE-MCT/data/woa_ini.nc'
    deg2rad = pi/180
    lon_c = 50
    lat_c = -83
    radius = 10.1

    id = Dataset(init_path, 'r')
    data1 = id.variables[var_name][0,0,:-15,:-2]
    id.close()

    id = Dataset(file_path, 'r')
    lon = id.variables['lon_rho'][:-15,:-2]
    lat = id.variables['lat_rho'][:-15,:-2]
    data2 = id.variables[var_name][-1,0,:-15,:-2]
    id.close()

    data_diff = data2 - data1

    # Convert grid to spherical coordinates
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)
    # Find centre in spherical coordinates
    x_c = -(lat_c+90)*cos(lon_c*deg2rad+pi/2)
    y_c = (lat_c+90)*sin(lon_c*deg2rad+pi/2)
    # Build a regular x-y grid and select the missing circle
    x_reg, y_reg = meshgrid(linspace(amin(x), amax(x), num=1000), linspace(amin(y), amax(y), num=1000))
    land_circle = zeros(shape(x_reg))
    land_circle = ma.masked_where(sqrt((x_reg-x_c)**2 + (y_reg-y_c)**2) > radius, land_circle)

    if colour_bounds is not None:
        # User has set bounds on colour scale
        lev = linspace(colour_bounds[0], colour_bounds[1], num=50)
        if colour_bounds[0] == -colour_bounds[1]:
            # Bounds are centered on zero, so choose a blue-to-red colourmap
            # centered on yellow
            colour_map = 'RdYlBu_r'
        else:
            colour_map = 'jet'
    else:
        max_val = amax(abs(data_diff))
        lev = linspace(-max_val, max_val, num=50)
        colour_map = 'RdYlBu_r'

    # Plot
    fig = figure(figsize=(16,12))
    ax = fig.add_subplot(1,1,1,aspect='equal')
    fig.patch.set_facecolor('white')
    # First shade everything in grey
    contourf(x, y, ones(shape(data_diff)), 1, colors=(('0.6', '0.6', '0.6')))
    # Fill in the missing circle
    contourf(x_reg, y_reg, land_circle, 1, colors=(('0.6', '0.6', '0.6')))
    # Now shade the temperature (land mask will remain grey)
    contourf(x, y, data_diff, lev, cmap=colour_map, extend='both')
    cbar = colorbar()
    cbar.ax.tick_params(labelsize=20)
    if var_name == 'temp':
        title(r'Difference in bottom water temperature since initialisation ($^{\circ}$C)', fontsize=30)
    elif var_name == 'salt':
        title(r'Difference in bottom water salinity since initialisation (psu)', fontsize=30)
    axis('off')

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


if __name__ == "__main__":

    file_path = raw_input("Path to ocean history/averages file: ")
    tmp = raw_input("Temperature (t) or salinity (s)? ")
    if tmp == 't':
        var_name = 'temp'
    elif tmp == 's':
        var_name = 'salt'
    colour_bounds = None
    get_bounds = raw_input("Set bounds on colour scale (y/n)? ")
    if get_bounds == 'y':
        lower_bound = float(raw_input("Lower bound: "))
        upper_bound = float(raw_input("Upper bound: "))
        colour_bounds = [lower_bound, upper_bound]
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("Filename for figure: ")
    else:
        save = False
        fig_name = None

    bw_diff_plot(file_path, var_name, colour_bounds, save, fig_name)

    
