from numpy import *
from netCDF4 import Dataset
from matplotlib.pyplot import *
from rotate_vector_cice import *
from seasonal_avg_cice import *

def compare_nic_seasonal (cice_file, var_name, colour_bounds=None, save=False, fig_name=None):

    nic_dir_head = '/g/data/gh5/access_om_025-CORE_NYF/output'
    output_number = 137
    nic_dir_tail = '/ice/HISTORY/'
    nic_file_head = 'iceh.0'
    nic_year_number = 133
    max_j = 300

    id = Dataset(cice_file, 'r')
    units = id.variables[var_name].units
    if var_name in ['uvel', 'vvel', 'strairx', 'strairy', 'strocnx', 'strocny']:
        rotate = True
        angle = id.variables['ANGLE'][:-15,:]
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
    grid_string = id.variables[var_name].coordinates
    if grid_string.startswith('ULON'):
        lon_name = 'ULON'
        lat_name = 'ULAT'
    else:
        lon_name = 'TLON'
        lat_name = 'TLAT'
    id.close()

    # Number of days in each month (this is just for Nic's output)
    ndays_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    # Season names for titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    # Degrees to radians conversion
    deg2rad = pi/180.0

    # Read the CICE grid
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
        this_cmp = seasonal_avg_cice(cice_file, var_name, [num_lat, num_lon])
        other_cmp = seasonal_avg_cice(cice_file, other_name, [num_lat, num_lon])
        if cmp_flag == 'x':
            cice_data_tmp, other_tmp = rotate_vector_cice(this_cmp, other_cmp, angle)
        elif cmp_flag == 'y':
            other_tmp, cice_data_tmp = rotate_vector_cice(other_cmp, this_cmp, angle)
    else:
        cice_data_tmp = seasonal_avg_cice(cice_file, var_name, [num_lat, num_lon])

    # Chop off northern boundary
    cice_data_tmp = cice_data_tmp[:,:-15,:]
    # Wrap the periodic boundary
    cice_data = ma.empty([size(cice_data_tmp,0), size(cice_data_tmp,1), size(cice_data_tmp,2)+1])
    cice_data[:,:,:-1] = cice_data_tmp
    cice_data[:,:,-1] = cice_data_tmp[:,:,0]

    # Read Nic's grid from the January file
    id = Dataset(nic_dir_head + str(output_number) + nic_dir_tail + nic_file_head + str(nic_year_number) + '-01.nc', 'r')
    nic_lon = id.variables[lon_name][:max_j,:]
    nic_lat = id.variables[lat_name][:max_j,:]
    id.close()

    nic_data = ma.empty([4, size(nic_lon,0), size(nic_lon,1)])
    nic_data[:,:,:] = 0.0
    for season in range(4):
        if season == 0:
            months = [12, 1, 2]
        elif season == 1:
            months = [3, 4, 5]
        elif season == 2:
            months = [6, 7, 8]
        elif season == 3:
            months = [9, 10, 11]
        season_days = 0

        for month in months:
            if month == 12:
                filename = nic_dir_head + str(output_number-1) + nic_dir_tail + nic_file_head + str(nic_year_number-1) + '-' + str(month) + '.nc'
            else:
                if month < 10:
                    filename = nic_dir_head + str(output_number) + nic_dir_tail + nic_file_head + str(nic_year_number) + '-0' + str(month) + '.nc'
                else:
                    filename = nic_dir_head + str(output_number) + nic_dir_tail + nic_file_head + str(nic_year_number) + '-' + str(month) + '.nc'
            id = Dataset(filename, 'r')
            nic_data[season,:,:] += id.variables[var_name][0,:max_j,:]*ndays_month[month-1]
            season_days += ndays_month[month-1]
        nic_data[season,:,:] /= season_days

    # Convert both grids to spherical coordinates
    cice_x = -(cice_lat+90)*cos(cice_lon*deg2rad+pi/2)
    cice_y = (cice_lat+90)*sin(cice_lon*deg2rad+pi/2)
    nic_x = -(nic_lat+90)*cos(nic_lon*deg2rad+pi/2)
    nic_y = (nic_lat+90)*sin(nic_lon*deg2rad+pi/2)

    bdry1 = -35
    bdry2 = 39
    bdry3 = -35
    bdry4 = 39

    if colour_bounds is not None:
        lev = linspace(colour_bounds[0], colour_bounds[1], num=50)
        if colour_bounds[0] == -colour_bounds[1]:
            colour_map = 'RdYlBu_r'
        else:
            colour_map = 'jet'
    else:
        if var_name in ['uvel', 'vvel', 'strairx', 'strairy', 'strocnx', 'strocny']:
            max_val = max(amax(abs(nic_data)), amax(abs(cice_data)))
            lev = linspace(-max_val, max_val, num=50)
            colour_map = 'RdYlBu_r'
        else:
            min_val = min(amin(nic_data), amin(cice_data))
            max_val = max(amax(nic_data), amax(cice_data))
            lev = linspace(min_val, max_val, num=50)
            colour_map = 'jet'

    fig = figure(figsize=(20,9))
    for season in range(4):
        ax = fig.add_subplot(2, 4, season+1, aspect='equal')
        contourf(nic_x, nic_y, nic_data[season,:,:], lev, cmap=colour_map, extend='both')
        if season == 0:
            text(-39, 0, 'Nic', fontsize=24, ha='right')
        title(season_names[season], fontsize=24)
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
        ax = fig.add_subplot(2, 4, season+5, aspect='equal')
        img = contourf(cice_x, cice_y, cice_data[season,:,:], lev, cmap=colour_map, extend='both')
        if season == 0:
            text(-39, 0, 'Me', fontsize=24, ha='right')
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
    cbaxes = fig.add_axes([0.25, 0.04, 0.5, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes)
    cbar.ax.tick_params(labelsize=16)
    suptitle(var_name + ' (' + units + ')', fontsize=30)
    subplots_adjust(wspace=0.025,hspace=0.025)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


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
    compare_nic_seasonal(cice_file, var_name, colour_bounds, save, fig_name)

    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                changes = raw_input("Enter a parameter to change: (1) file path, (2) variable name, (3) colour bounds, (4) save/display; or enter to continue: ")
                if len(changes) == 0:
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
                fig_name = raw_input("File name for figure: ")
            compare_nic_seasonal(cice_file, var_name, colour_bounds, save, fig_name)
        else:
            break


    
