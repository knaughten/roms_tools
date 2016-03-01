from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from calc_z import *

# Make a circumpolar Antarctic plot of the given (horizontal) variable from ROMS.
# Input:
# grid_path = path to ROMS grid file
# file_path = path to ocean history/averages file
# var_name = name of variable in file_path to plot
# tstep = timestep in file_path to plot (1-indexed)
# depth_key = integer flag indicating whether to plot the surface level (0), the bottom level 
#             (1), a specific depth to interpolate to (2), the vertical average throughout the 
#             entire water column (3), or the vertical average between two specific depths (4)
# depth = if depth_key=2, the specific depth to interpolate to (negative, in metres)
# depth_bounds = if depth_key=4, the specific depths to average between (negative, in metres),
#                stored as an array of size 2 with the shallow bound first.
# save = optional boolean flag indicating that the plot should be saved to a file rather than 
#        displayed on the screen
# fig_name = if save=True, filename for figure
def circumpolar_plot (grid_path, file_path, var_name, tstep, depth_key, depth, depth_bounds, save=False, fig_name=None):

    # Grid parameters
    theta_s = 0.9
    theta_b = 4.0
    hc = 40
    N = 31
    deg2rad = pi/180

    # Read the variable and figure out if 2D or 3D (not including time)
    id = Dataset(file_path, 'r')
    if len(id.variables[var_name].shape) == 4:
        # 3D variable; will have to choose depth later
        data_full = id.variables[var_name][tstep-1,:,:,:]
        choose_depth = True
    elif len(id.variables[var_name].shape) == 3:
        # 2D variable
        data = id.variables[var_name][tstep-1,:,:]
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

    # Figure out what grid the variable is on
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
    # h and zice are on the rho-grid; interpolate if necessary
    if grid_name == 'u':
        h = 0.5*(h[:,0:-1] + h[:,1:])
        zice = 0.5*(zice[:,0:-1] + zice[:,1:])
    elif grid_name == 'v':
        h = 0.5*(h[0:-1,:] + h[1:,:])
        zice = 0.5*(zice[0:-1,:] + zice[1:,:])
    # Read the correct lat and lon for this grid
    lon = id.variables[lon_name][:,:]
    lat = id.variables[lat_name][:,:]
    id.close()

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
            depth_string = 'at '+str(round(-depth))+' m'
        elif depth_key == 3:
            depth_string = 'vertically averaged'
        elif depth_key == 4:
            depth_string = 'vertically averaged between '+str(round(-depth_bounds[0]))+' and '+str(round(-depth_bounds[1]))+' m'
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
            # We will need to calculate z-coordinates
            z, sc_r, Cs_r = calc_z(h, zice, lon, lat, theta_s, theta_b, hc, N)
            if depth_key == 2:
                # Interpolate to given depth
                data = interp_depth(data_full, z, depth)
            elif depth_key in [3, 4]:
                # We need dz for averaging
                # We have z on the midpoint of each cell, now find it on the top and bottom edges
                # of each cell
                z_edges = zeros((size(z,0)+1, size(z,1), size(z,2)))
                z_edges[1:-1,:,:] = 0.5*(z[0:-1,:,:] + z[1:,:,:])
                # At the surface, z=zice; at bottom, extrapolate
                z_edges[-1,:,:] = zice[:,:]
                z_edges[0,:,:] = 2*z[0,:,:] - z_edges[1,:,:]
                # Now find dz
                dz = z_edges[1:,:,:] - z_edges[0:-1,:,:]
                if depth_key == 3:
                    # Vertically average entire water column
                    data = sum(data_full*dz, axis=0)/sum(dz, axis=0)
                elif depth_key == 4:                    
                    # Vertically average between given depths
                    data = average_btw_depths (data_full, z, dz, depth_bounds)

    # Center levels on 0 for certain variables, with a blue-to-red colourmap
    if var_name in ['u', 'v', 'ubar', 'vbar', 'm']:
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
    title(long_name+' ('+units+')\n'+depth_string, fontsize=30)
    axis('off')

    if save:
        savefig(fig_name)
    else:
        show()


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

    grid_path = raw_input("Path to ROMS grid file: ")
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

    # Get index of time axis in ROMS history/averages file
    tstep = int(raw_input("Timestep number (starting at 1): "))

    # Get save/display choice
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None

    # Make the plot
    circumpolar_plot(grid_path, file_path, var_name, tstep, depth_key, depth, depth_bounds, save, fig_name)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the input parameters; repeat until the user is finished
                changes = raw_input("Enter a parameter to change: (1) grid path, (2) file path, (3) variable name, (4) depth, (5) timestep number, (6) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters.
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
                    elif int(changes) == 4:
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
                    elif int(changes) == 5:
                        # New timestep
                        tstep = int(raw_input("Timestep number (starting at 1): "))
                    elif int(changes) == 6:
                        # Change from display to save, or vice versa
                        save = not save
            if save:
                # Get file name for figure
                fig_name = raw_input("File name for figure: ")

            # Make the plot
            circumpolar_plot(grid_path, file_path, var_name, tstep, depth_key, depth, depth_bounds, save, fig_name)

        else:
            break

            
                    
        
        

    
    
