from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from cartesian_grid_3d import *
from rotate_vector_roms import *

# Make a circumpolar Antarctic plot of the given (horizontal) ROMS variable.
# Input:
# file_path = path to ocean history/averages file
# var_name = name of variable in file_path to plot
# tstep = timestep in file_path to plot (1-indexed)
# depth_key = integer flag indicating whether to plot the surface level (0), 
#             the bottom level (1), a specific depth to interpolate to (2), 
#             the vertical average throughout the entire water column (3), 
#             or the vertical average between two specific depths (4)
# depth = if depth_key=2, the specific depth to interpolate to (negative, 
#         in metres)
# depth_bounds = if depth_key=4, the specific depths to average between 
#                (negative, in metres), stored as an array of size 2 with the 
#                shallow bound first.
# colour_bounds = optional bounds on colour scale, stored as an array of size
#                 2 with the lower bound first. If colour_bounds = None, then
#                 determine colour scale bounds automatically.
# save = optional boolean flag indicating that the plot should be saved to a 
#        file rather than displayed on the screen
# fig_name = if save=True, filename for figure
# grid_path = path to grid file; only needed if var_name is a vector component
def circumpolar_plot (file_path, var_name, tstep, depth_key, depth, depth_bounds, colour_bounds=None, save=False, fig_name=None, grid_path=None):

    # Grid parameters
    theta_s = 7.0
    theta_b = 2.0
    hc = 250
    N = 31
    deg2rad = pi/180

    # Read the variable and figure out if 2D or 3D (not including time)
    id = Dataset(file_path, 'r')
    if len(id.variables[var_name].shape) == 4:
        # 3D variable; will have to choose depth later
        data_full = id.variables[var_name][tstep-1,:,:-15,:]
        choose_depth = True
    elif len(id.variables[var_name].shape) == 3:
        # 2D variable
        data = id.variables[var_name][tstep-1,:-15,:]
        choose_depth = False
    if var_name == 'salt':
        units = 'psu'
    elif var_name == 'm':
        # Convert ice shelf melt rate from m/s to m/yr
        units = 'm/year'
        data = data*60.*60.*24.*365.25
    else:
        units = id.variables[var_name].units
    long_name = id.variables[var_name].long_name

    # Check for vector variables that need to be rotated
    if var_name in ['ubar', 'vbar', 'u', 'v', 'sustr', 'svstr', 'bustr', 'bvstr']:
        grid_id = Dataset(grid_path, 'r')
        angle = grid_id.variables['angle'][:-15,:]
        grid_id.close()
        if var_name in ['ubar', 'sustr', 'bustr']:
            # 2D u-variable
            u_data = data[:,:]
            v_data = id.variables[var_name.replace('u','v')][tstep-1,:-15,:]
            u_data_lonlat, v_data_lonlat = rotate_vector_roms(u_data, v_data, angle)
            data = u_data_lonlat
        elif var_name in ['vbar', 'svstr', 'bvstr']:
            # 2D v-variable
            v_data = data[:,:]
            u_data = id.variables[var_name.replace('v','u')][tstep-1,:-15,:]
            u_data_lonlat, v_data_lonlat = rotate_vector_roms(u_data, v_data, angle)
            data = v_data_lonlat
        elif var_name in ['u']:
            # 3D u-variable
            data_full_ugrid = data_full[:,:,:]
            data_full = ma.empty([data_full_ugrid.shape[0],data_full_ugrid.shape[1],data_full_ugrid.shape[2]+1])
            for k in range(N):
                u_data = data_full_ugrid[k,:,:]
                v_data = id.variables[var_name.replace('u','v')][tstep-1,k,:-15,:]
                u_data_lonlat, v_data_lonlat = rotate_vector_roms(u_data, v_data, angle)
                data_full[k,:,:] = u_data_lonlat
        elif var_name in ['v']:
            # 3D v-variable
            data_full_vgrid = data_full[:,:,:]
            data_full = ma.empty([data_full_vgrid.shape[0],data_full_vgrid.shape[1]+1,data_full_vgrid.shape[2]])
            for k in range(N):
                v_data = data_full_vgrid[k,:,:]
                u_data = id.variables[var_name.replace('v','u')][tstep-1,k,:-15,:]
                u_data_lonlat, v_data_lonlat = rotate_vector_roms(u_data, v_data, angle)
                data_full[k,:,:] = v_data_lonlat        

    id.close()
    id = Dataset(grid_path, 'r')
    # Read grid variables
    h = id.variables['h'][:-15,:]
    zice = id.variables['zice'][:-15,:]
    lon = id.variables['lon_rho'][:-15,:]
    lat = id.variables['lat_rho'][:-15,:]
    id.close()

    # Throw away the overlapping periodic boundary
    if choose_depth:
        data_full = data_full[:,:,:-1]
    else:
        data = data[:,:-1]
    lon = lon[:,:-1]
    lat = lat[:,:-1]
    h = h[:,:-1]
    zice = zice[:,:-1]

    # Convert to spherical coordinates
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)

    # Choose what to write on the title about depth
    if choose_depth:
        if depth_key == 0:
            depth_string = 'at surface'
        elif depth_key == 1:
            depth_string = 'at bottom'
        elif depth_key == 2:
            depth_string = 'at '+str(int(round(-depth)))+' m'
        elif depth_key == 3:
            depth_string = 'vertically averaged'
        elif depth_key == 4:
            depth_string = 'vertically averaged between '+str(int(round(-depth_bounds[0])))+' and '+str(int(round(-depth_bounds[1])))+' m'
    else:
        depth_string = ''

    if choose_depth:
        # For 3D variables, select data corresponding to depth choice
        if depth_key == 0:
            # Surface layer
            data = data_full[-1,:,:]
        elif depth_key == 1:
            # Bottom layer
            data = data_full[0,:,:]
        else:
            # We will need z-coordinates and possibly dz
            dx, dy, dz, z = cartesian_grid_3d(lon, lat, h, zice, theta_s, theta_b, hc, N)
            if depth_key == 2:
                # Interpolate to given depth
                data = interp_depth(data_full, z, depth)
            elif depth_key == 3:
                # Vertically average entire water column
                data = sum(data_full*dz, axis=0)/sum(dz, axis=0)
            elif depth_key == 4:
                # Vertically average between given depths
                data = average_btw_depths(data_full, z, dz, depth_bounds)

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
        if var_name in ['u', 'v', 'ubar', 'vbar', 'm', 'shflux', 'ssflux', 'sustr', 'svstr', 'bustr', 'bvstr', 'ssflux_restoring']:
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
    title(long_name+' ('+units+')\n'+depth_string, fontsize=30)
    axis('off')

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Linearly interpolate data to the specified depth at each horizontal point.
# Input:
# data_3d = array of data, dimension depth x lat x lon
# z_3d = array of depth values (negative, in metres), dimension depth x lat x lon
# z0 = depth to interpolate to (negative, in metres)
# Output:
# data = array of data interpolated to z0, dimension lat x lon
def interp_depth (data_3d, z_3d, z0):

    # Save horizontal dimensions
    num_lat = size(data_3d, 1)
    num_lon = size(data_3d, 2)
    # Set up output array; initialise this way to get the surface land mask
    data = data_3d[-1,:,:]*0.0

    # Loop over each horizontal point; can't find a cleaner way to do this
    for j in range(num_lat):
        for i in range(num_lon):
            # Extract the data and depth values of the current water column
            z_col = z_3d[:,j,i]
            data_col = data_3d[:,j,i]
            if data[j,i] is ma.masked:
                # This is a land point; leave as is                
                pass
            elif all(z_col < z0):
                # z0 is too shallow (i.e. in an ice shelf)
                data[j,i] = ma.masked
            elif all(z_col > z0):
                # z0 is too deep (i.e. in the seafloor)
                data[j,i] = ma.masked
            else:
                # Find the first index (starting from the bottom) shallower than z0
                k_above = nonzero(z_col > z0)[0][0]
                # The index before it will be the last index deeper than z0
                k_below = k_above - 1
                # Linearly interpolate data to z0
                coeff1 = (z_col[k_below] - z0)/(z_col[k_below] - z_col[k_above])
                coeff2 = 1 - coeff1
                data[j,i] = coeff1*data_col[k_above] + coeff2*data_col[k_below]

    return data


# Vertically average data between the specified depths at each horizontal point.
# Input:
# data_3d = array of data, dimension depth x lat x lon
# z_3d = array of depth values (negative, in metres), dimension depth x lat x lon
# dz_3d = array of vertical cell thicknesses (positive, in metres), dimension depth x lat x lon
# z_bounds = array containing the two depth values to average between (shallower depth first)
# Output:
# data = array of data averaged between z_bounds, dimension lat x lon
def average_btw_depths (data_3d, z_3d, dz_3d, z_bounds):

    # Save horizontal dimensions
    num_lat = size(data_3d, 1)
    num_lon = size(data_3d, 2)
    # Set up output array; initialise this way to get the surface land mask
    data = data_3d[-1,:,:]*0.0
    # Unpack z_bounds
    z_shallow = z_bounds[0]
    z_deep = z_bounds[1]

    # Loop over each horizontal point; can't find a cleaner way to do this
    for j in range(num_lat):
        for i in range(num_lon):
            # Extract the data, depth, and thickness values of the current water column
            z_col = z_3d[:,j,i]
            dz_col = dz_3d[:,j,i]
            data_col = data_3d[:,j,i]
            if data[j,i] is ma.masked:
                # This is a land point, leave as is
                pass
            elif all(z_col < z_deep):
                # both depth bounds are too shallow (i.e. in an ice shelf)
                data[j,i] = ma.masked
            elif all(z_col > z_shallow):
                # both depth bounds are too deep (i.e. in the seafloor)
                data[j,i] = ma.masked
            else:
                # Initialise integral and depth_range (vertical distance we're integrating over) to 0
                integral = 0.0
                depth_range = 0.0
                if any(z_col < z_deep):
                    # There exist ocean cells below z_deep
                    # Linearly interpolate to z_deep
                    k_above_deep = nonzero(z_col > z_deep)[0][0]
                    k_below_deep = k_above_deep - 1
                    coeff1 = (z_col[k_below_deep] - z_deep)/(z_col[k_below_deep] - z_col[k_above_deep])
                    coeff2 = 1 - coeff1
                    data_deep = coeff1*data_col[k_above_deep] + coeff2*data_col[k_below_deep]
                    # Now integrate between z_deep and z_col[k_above_deep]
                    dz_curr = z_col[k_above_deep] - z_deep
                    integral += 0.5*(data_deep + data_col[k_above_deep])*dz_curr # Trapezoidal rule
                    depth_range += dz_curr
                    # Save index of k_above_deep; we will start normal (non-interpolated) integration here
                    k_start = k_above_deep
                else:
                    # z_deep is deeper than the seafloor at this location
                    # Start normal (non-interpolated) integration at the seafloor
                    k_start = 0
                if any(z_col > z_shallow):
                    # There exist ocean cells above z_shallow
                    # Linearly interpolate to z_shallow
                    k_above_shallow = nonzero(z_col > z_shallow)[0][0]
                    k_below_shallow = k_above_shallow - 1
                    coeff1 = (z_col[k_below_shallow] - z_shallow)/(z_col[k_below_shallow] - z_col[k_above_shallow])
                    coeff2 = 1 - coeff1
                    data_shallow = coeff1*data_col[k_above_shallow] + coeff2*data_col[k_below_shallow]
                    # Now integrate between z_col[k_below_shallow] and z_shallow
                    dz_curr = z_shallow - z_col[k_below_shallow]
                    integral += 0.5*(data_col[k_below_shallow] + data_shallow)*dz_curr # Trapezoidal rule
                    depth_range += dz_curr
                    # Save index of k_above_shallow; we will stop normal (non-interpolated) integration one level below it
                    k_end = k_above_shallow
                else:
                    # z_shallow is shallower than the sea surface (or ice shelf draft) at this location
                    # Continue normal (non-interpolated) integration all the way to the surface/ice shelf
                    k_end = size(z_col)
                # Now integrate between k_start and k_end
                if k_start < k_end:
                    integral += sum(data_col[k_start:k_end]*dz_col[k_start:k_end])
                    depth_range += sum(dz_col[k_start:k_end])
                # Divide integral by depth_range to get the vertical average
                data[j,i] = integral/depth_range

    return data    


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to ocean history/averages file: ")
    var_name = raw_input("Variable name: ")

    # Figure out if we need to ask for depth information
    id = Dataset(file_path, 'r')
    if len(id.variables[var_name].shape) == 4:
        # 3D variable; ask for depth information
        depth_type = raw_input("Single depth (s) or vertical average (v)? ")
        if depth_type == 's':
            depth_input = raw_input("Surface layer (s), bottom layer (b), or specific depth (d)? ")
            if depth_input == 's':
                depth_key = 0
                depth = NaN
                depth_bounds = None
            elif depth_input == 'b':
                depth_key = 1
                depth = NaN
                depth_bounds = None
            elif depth_input == 'd':
                depth_key = 2
                depth = -1*float(raw_input("Enter depth (positive, in metres): "))
                depth_bounds = None
        elif depth_type == 'v':
            depth_input = raw_input("Vertical average throughout the entire water column (w) or between two specific depths (d)? ")
            if depth_input == 'w':
                depth_key = 3
                depth = NaN
                depth_bounds = None
            elif depth_input == 'd':
                depth_key = 4
                depth = NaN
                shallow_bound = -1*float(raw_input("Enter shallow depth bound (positive, in metres): "))
                deep_bound = -1*float(raw_input("Enter deep depth bound (positive, in metres): "))
                depth_bounds = [shallow_bound, deep_bound]
    elif len(id.variables[var_name].shape) == 3:
        # 2D variable; depth doesn't apply
        depth_key = 0
        depth = NaN
        depth_bounds = None
    id.close()

    if var_name in ['ubar', 'vbar', 'u', 'v', 'sustr', 'svstr', 'bustr', 'bvstr']:
        # Will need the grid file to get the angle
        grid_path = raw_input("Path to ROMS grid file: ")
    else:
        grid_path = raw_input("Path to ROMS grid file: ") #grid_path = None

    # Get index of time axis in ROMS history/averages file
    tstep = int(raw_input("Timestep number (starting at 1): "))

    # Get colour bounds if necessary
    colour_bounds = None
    get_bounds = raw_input("Set bounds on colour scale (y/n)? ")
    if get_bounds == 'y':
        lower_bound = float(raw_input("Lower bound: "))
        upper_bound = float(raw_input("Upper bound: "))
        colour_bounds = [lower_bound, upper_bound]

    # Get save/display choice
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None

    # Make the plot
    circumpolar_plot(file_path, var_name, tstep, depth_key, depth, depth_bounds, colour_bounds, save, fig_name, grid_path)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the input parameters; repeat until the user is finished
                changes = raw_input("Enter a parameter to change: (1) file path, (2) variable name, (3) depth, (4) timestep number, (5) colour bounds, (6) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters.
                    break
                else:
                    if int(changes) == 1:
                        # New file path
                        file_path = raw_input("Path to ocean history/averages file: ")
                    elif int(changes) == 2:
                        # New variable name
                        var_name = raw_input("Variable name: ")
                        # Figure out if we need to ask for depth information
                        id = Dataset(file_path, 'r')
                        if len(id.variables[var_name].shape) == 4:
                            # 3D variable; ask  for depth information
                            depth_type = raw_input("Single depth (s) or vertical average (v)? ")
                            if depth_type == 's':
                                depth_input = raw_input("Surface layer (s), bottom layer (b), or specific depth (d)? ")
                                if depth_input == 's':
                                    depth_key = 0
                                    depth = NaN
                                    depth_bounds = None
                                elif depth_input == 'b':
                                    depth_key = 1
                                    depth = NaN
                                    depth_bounds = None
                                elif depth_input == 'd':
                                    depth_key = 2
                                    depth = -1*float(raw_input("Enter depth (positive, in metres): "))
                                    depth_bounds = None
                            elif depth_type == 'v':
                                depth_input = raw_input("Vertical average throughout the entire water column (w) or between two specific depths (d)? ")
                                if depth_input == 'w':
                                    depth_key = 3
                                    depth = NaN
                                    depth_bounds = None
                                elif depth_input == 'd':
                                    depth_key = 4
                                    depth = NaN
                                    shallow_bound = -1*float(raw_input("Enter shallow depth bound (positive, in metres): "))
                                    deep_bound = -1*float(raw_input("Enter deep depth bound (positive, in metres): "))
                                    depth_bounds = [shallow_bound, deep_bound]
                        elif len(id.variables[var_name].shape) == 3:
                            # 2D variable; depth doesn't apply
                            depth_key = 0
                            depth = NaN
                            depth_bounds = None
                        id.close()
                        if var_name in ['ubar', 'vbar', 'u', 'v', 'sustr', 'svstr', 'bustr', 'bvstr'] and grid_path is None:
                            # Will need the grid file to get the angle
                            grid_path = raw_input("Path to ROMS grid file: ")
                    elif int(changes) == 3:
                        # New depth information
                        depth_type = raw_input("Single depth (s) or vertical average (v)? ")
                        if depth_type == 's':
                            depth_input = raw_input("Surface layer (s), bottom layer (b), or specific depth (d)? ")
                            if depth_input == 's':
                                depth_key = 0
                                depth = NaN
                                depth_bounds = None
                            elif depth_input == 'b':
                                depth_key = 1
                                depth = NaN
                                depth_bounds = None
                            elif depth_input == 'd':
                                depth_key = 2
                                depth = -1*float(raw_input("Enter depth (positive, in metres): "))
                                depth_bounds = None
                        elif depth_type == 'v':
                            depth_input = raw_input("Vertical average throughout the entire water column (w) or between two specific depths (d)? ")
                            if depth_input == 'w':
                                depth_key = 3
                                depth = NaN
                                depth_bounds = None
                            elif depth_input == 'd':
                                depth_key = 4
                                depth = NaN
                                shallow_bound = -1*float(raw_input("Enter shallow depth bound (positive, in metres): "))
                                deep_bound = -1*float(raw_input("Enter deep depth bound (positive, in metres): "))
                                depth_bounds = [shallow_bound, deep_bound]
                    elif int(changes) == 4:
                        # New timestep
                        tstep = int(raw_input("Timestep number (starting at 1): "))
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
            circumpolar_plot(file_path, var_name, tstep, depth_key, depth, depth_bounds, colour_bounds, save, fig_name, grid_path)

        else:
            break

            
                    
        
        

    
    
