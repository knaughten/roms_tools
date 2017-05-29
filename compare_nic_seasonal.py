from numpy import *
from netCDF4 import Dataset
from matplotlib.pyplot import *
from rotate_vector_cice import *
from seasonal_avg_cice import *

# Make a 4x2 plot comparing seasonal averages of the given variable from my
# implementation of CICE in MetROMS (bottom) and Nic Hannah's implementation
# coupled to MOM (top).
# Input:
# cice_file = path to CICE file from MetROMS simulation, containing at least
#             one complete December-November period, in 5-day averages. If
#             there are multiple such periods, the last one will be used.
# var_name = variable name to plot
# colour_bounds = optional array of size 2 containing bounds on colour scale
# save = optional boolean indicating to save the file rather than display it
#        on screen
# fig_name = if save=True, filename for figure
def compare_nic_seasonal (cice_file, var_name, colour_bounds=None, save=False, fig_name=None):

    # Piece together the paths to Nic's monthly averaged output on raijin
    nic_dir_head = '/g/data/gh5/access_om_025-CORE_NYF/output'
    output_number = 137
    nic_dir_tail = '/ice/HISTORY/'
    nic_file_head = 'iceh.0'
    nic_year_number = 133
    # Maximum j-index to read in Nic's output
    max_j = 300

    # Look at the variable in my output
    id = Dataset(cice_file, 'r')
    # Save units
    units = id.variables[var_name].units
    # Check if this is a vector we need to rotate to lon-lat space
    if var_name in ['uvel', 'vvel', 'strairx', 'strairy', 'strocnx', 'strocny']:
        rotate = True
        # Read the angle of grid rotation
        angle = id.variables['ANGLE'][:,:]
        # Figure out whether this is an x or y component, and what the name of
        # the other component is
        if var_name == 'uvel':
            cmp_flag = 'x'
            other_name = 'vvel'
        elif var_name == 'vvel':
            cmp_flag = 'y'
            other_name = 'uvel'
        elif var_name in ['strairx', 'strocnx']:
            cmp_flag = 'x'
            other_name = var_name.replace('x', 'y')
        elif var_name in ['strairy', 'strocny']:
            cmp_flag = 'y'
            other_name = var_name.replace('x', 'y')            
    else:
        rotate = False
    # Read the correct grid (tracer or velocity)
    grid_string = id.variables[var_name].coordinates
    if grid_string.startswith('ULON'):
        lon_name = 'ULON'
        lat_name = 'ULAT'
    else:
        lon_name = 'TLON'
        lat_name = 'TLAT'
    id.close()

    # Number of days in each month (this is just for Nic's output)
    # Note Nic doesn't run with leap years
    ndays_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    # Season names for titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    # Degrees to radians conversion
    deg2rad = pi/180.0

    # Read my CICE grid
    id = Dataset(cice_file, 'r')
    cice_lon_tmp = id.variables[lon_name][:-15,:]
    cice_lat_tmp = id.variables[lat_name][:-15,:]
    num_lon = id.variables['TLON'].shape[1]
    num_lat = id.variables['TLAT'].shape[0]
    id.close()
    # Wrap the periodic boundary by 1 cell
    cice_lon = ma.empty([size(cice_lon_tmp,0), size(cice_lon_tmp,1)+1])
    cice_lat = ma.empty([size(cice_lat_tmp,0), size(cice_lat_tmp,1)+1])
    cice_lon[:,:-1] = cice_lon_tmp
    cice_lon[:,-1] = cice_lon_tmp[:,0]
    cice_lat[:,:-1] = cice_lat_tmp
    cice_lat[:,-1] = cice_lat_tmp[:,0]

    # Get seasonal averages of CICE data
    if rotate:
        # Average both components of the vector
        this_cmp = seasonal_avg_cice(cice_file, var_name, [num_lat, num_lon])
        other_cmp = seasonal_avg_cice(cice_file, other_name, [num_lat, num_lon])
        # Rotate to lon-lat space
        if cmp_flag == 'x':
            cice_data_tmp = ma.empty(shape(this_cmp))
            for season in range(4):
                tmp1, tmp2 = rotate_vector_cice(this_cmp[season,:,:], other_cmp[season,:,:], angle)
                cice_data_tmp[season,:,:] = tmp1
        elif cmp_flag == 'y':
            cice_data_tmp = ma.empty(shape(this_cmp))
            for season in range(4):
                tmp1, tmp2 = rotate_vector_cice(other_cmp[season,:,:], this_cmp[season,:,:], angle)
                cice_data_tmp[season,:,:] = tmp2
    else:
        cice_data_tmp = seasonal_avg_cice(cice_file, var_name, [num_lat, num_lon])

    # Chop off northern boundary
    cice_data_tmp = cice_data_tmp[:,:-15,:]
    # Wrap the periodic boundary
    cice_data = ma.empty([size(cice_data_tmp,0), size(cice_data_tmp,1), size(cice_data_tmp,2)+1])
    cice_data[:,:,:-1] = cice_data_tmp
    cice_data[:,:,-1] = cice_data_tmp[:,:,0]

    # Conversions to account for different thermodynamics schemes
    if var_name in ['frazil', 'snoice']:
        cice_data /= 3.6

    # Read Nic's grid from the January file
    id = Dataset(nic_dir_head + str(output_number) + nic_dir_tail + nic_file_head + str(nic_year_number) + '-01.nc', 'r')
    nic_lon = id.variables[lon_name][:max_j,:]
    nic_lat = id.variables[lat_name][:max_j,:]
    id.close()

    # Get seasonal averages of Nic's output
    nic_data = ma.empty([4, size(nic_lon,0), size(nic_lon,1)])
    nic_data[:,:,:] = 0.0
    # Loop over seasons
    for season in range(4):
        # Figure out what months we care about for this season
        if season == 0:
            months = [12, 1, 2]
        elif season == 1:
            months = [3, 4, 5]
        elif season == 2:
            months = [6, 7, 8]
        elif season == 3:
            months = [9, 10, 11]
        # Days in season so far
        season_days = 0
        # Loop over months
        for month in months:
            if month == 12:
                # Read December from the previous year
                filename = nic_dir_head + str(output_number-1) + nic_dir_tail + nic_file_head + str(nic_year_number-1) + '-' + str(month) + '.nc'
            else:
                if month < 10:
                    filename = nic_dir_head + str(output_number) + nic_dir_tail + nic_file_head + str(nic_year_number) + '-0' + str(month) + '.nc'
                else:
                    filename = nic_dir_head + str(output_number) + nic_dir_tail + nic_file_head + str(nic_year_number) + '-' + str(month) + '.nc'
            id = Dataset(filename, 'r')
            # Integrate over time
            nic_data[season,:,:] += id.variables[var_name][0,:max_j,:]*ndays_month[month-1]
            season_days += ndays_month[month-1]
        # Convert from integral to average
        nic_data[season,:,:] /= season_days

    # Convert both grids to spherical coordinates
    cice_x = -(cice_lat+90)*cos(cice_lon*deg2rad+pi/2)
    cice_y = (cice_lat+90)*sin(cice_lon*deg2rad+pi/2)
    nic_x = -(nic_lat+90)*cos(nic_lon*deg2rad+pi/2)
    nic_y = (nic_lat+90)*sin(nic_lon*deg2rad+pi/2)

    # Boundaries on plots
    bdry1 = -35
    bdry2 = 39
    bdry3 = -35
    bdry4 = 39

    if colour_bounds is not None:
        # User-defined colour bounds
        lev = linspace(colour_bounds[0], colour_bounds[1], num=50)
        if colour_bounds[0] == -colour_bounds[1]:
            # Centered on zero; use red-yellow-blue colourmap
            colour_map = 'RdYlBu_r'
        else:
            # Not centered on zero; go for a rainbow
            colour_map = 'jet'
    else:
        # Automatic colour bounds based on min/max of the data
        if var_name in ['uvel', 'vvel', 'strairx', 'strairy', 'strocnx', 'strocny']:
            # Center on zero and use red-yellow-blue colourmap
            max_val = max(amax(abs(nic_data)), amax(abs(cice_data)))
            lev = linspace(-max_val, max_val, num=50)
            colour_map = 'RdYlBu_r'
        else:
            # Not centered on zero; go for a rainbow
            min_val = min(amin(nic_data), amin(cice_data))
            max_val = max(amax(nic_data), amax(cice_data))
            lev = linspace(min_val, max_val, num=50)
            colour_map = 'jet'

    # Make the figure
    fig = figure(figsize=(20,9))
    # Loop over seasons
    for season in range(4):
        # Nic's output
        ax = fig.add_subplot(2, 4, season+1, aspect='equal')
        contourf(nic_x, nic_y, nic_data[season,:,:], lev, cmap=colour_map, extend='both')
        if season == 0:
            text(-39, 0, 'Nic', fontsize=24, ha='right')
        title(season_names[season], fontsize=24)
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
        # My output
        ax = fig.add_subplot(2, 4, season+5, aspect='equal')
        img = contourf(cice_x, cice_y, cice_data[season,:,:], lev, cmap=colour_map, extend='both')
        if season == 0:
            text(-39, 0, 'Me', fontsize=24, ha='right')
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
    # Add colourbar on the bottom
    cbaxes = fig.add_axes([0.25, 0.04, 0.5, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes)
    cbar.ax.tick_params(labelsize=16)
    # Main title with variable name and units
    suptitle(var_name + ' (' + units + ')', fontsize=30)
    subplots_adjust(wspace=0.025,hspace=0.025)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    cice_file = raw_input("Path to CICE file, containing at least one complete Dec-Nov period: ")
    var_name = raw_input("Variable name: ")
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
    compare_nic_seasonal(cice_file, var_name, colour_bounds, save, fig_name)

    # Repeat until the user is finished
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            # Ask for changes to input parameters until the user is finished
            while True:
                changes = raw_input("Enter a parameter to change: (1) file path, (2) variable name, (3) colour bounds, (4) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        cice_file = raw_input("Path to CICE file, containing at least one complete Dec-Nov period: ")
                    elif int(changes) == 2:
                        var_name = raw_input("Variable name: ")
                    elif int(changes) == 3:
                        colour_bounds = None
                        get_bounds = raw_input("Set bounds on colour scale (y/n)? ")
                        if get_bounds == 'y':
                            lower_bound = float(raw_input("Lower bound: "))
                            upper_bound = float(raw_input("Upper bound: "))
                            colour_bounds = [lower_bound, upper_bound]
                    elif int(changes) == 4:
                        save = not save
            if save:
                # Get new figure name
                fig_name = raw_input("File name for figure: ")
            # Make the plot
            compare_nic_seasonal(cice_file, var_name, colour_bounds, save, fig_name)
        else:
            # No more plots
            break


    
