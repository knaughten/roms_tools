from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from calc_z import *
from rotate_vector_roms import *

# Create a zonally averaged or zonally sliced plot (i.e. depth vs latitude)
# of the given variable.
# Input:
# file_path = path to ocean history or output file
# var_name = variable name in file_path; must be a variable of dimension
#            time x depth x latitude x longitude, and on the rho, u, or v grid
# tstep = timestep in file_path to plot (1-indexed)
# lon_key = integer flag indicating whether to plot a single longitude (0),
#           the zonal average over all longitudes (1), or the zonal average
#           between two specific longitudes (2)
# lon0 = if lon_key=0, the specific longitude to plot (between -180 and 180)
# lon_bounds = if lon_key=2, the specific longitudes to average between,
#              stored as an array of size 2 with the western bound first
# depth_min = deepest depth to plot (negative, metres)
# colour_bounds = optional bounds on colour scale, stored as an array of size
#                 2 with the lower bound first. If colour_bounds = None, then
#                 determine colour scale bounds automatically.
# save = optional boolean flag; if True, the figure will be saved with file name
#        fig_name, if False, the figure will display on the screen
# fig_name = optional string containing filename for figure, if save=True
# grid_path = path to grid file; only needed if var_name is a vector component
def zonal_plot (file_path, var_name, tstep, lon_key, lon0, lon_bounds, depth_min, colour_bounds=None, save=False, fig_name=None, grid_path=None):

    # Grid parameters
    theta_s = 0.9
    theta_b = 4.0
    hc = 40
    N = 31

    # Read the variable
    id = Dataset(file_path, 'r')
    data_3d = id.variables[var_name][tstep-1,:,:-15,:]
    # Also read sea surface height
    zeta = id.variables['zeta'][tstep-1,:-15,:]
    if var_name == 'salt':
        units = 'psu'
    else:
        units = id.variables[var_name].units
    long_name = id.variables[var_name].long_name

    # Rotate velocity if necessary
    if var_name in ['u', 'v']:
        grid_id = Dataset(grid_path, 'r')
        angle = grid_id.variables['angle'][:-15,:]
        grid_id.close()
        if var_name == 'u':
            data_3d_ugrid = data_3d[:,:,:]
            data_3d = ma.empty([data_3d_ugrid.shape[0], data_3d_ugrid.shape[1], data_3d_ugrid.shape[2]+1])
            for k in range(N):
                u_data = data_3d_ugrid[k,:,:]
                v_data = id.variables['v'][tstep-1,k,:-15,:]
                u_data_lonlat, v_data_lonlat = rotate_vector_roms(u_data, v_data, angle)
                data_3d[k,:,:] = u_data_lonlat
        elif var_name == 'v':
            data_3d_vgrid = data_3d[:,:,:]
            data_3d = ma.empty([data_3d_vgrid.shape[0], data_3d_vgrid.shape[1]+1, data_3d_vgrid.shape[2]])
            for k in range(N):
                v_data = data_3d_vgrid[k,:,:]
                u_data = id.variables['u'][tstep-1,k,:-15,:]
                u_data_lonlat, v_data_lonlat = rotate_vector_roms(u_data, v_data, angle)
                data_3d[k,:,:] = v_data_lonlat

    # Read grid variables
    h = id.variables['h'][:-15,:]
    zice = id.variables['zice'][:-15,:]
    lon_2d = id.variables['lon_rho'][:-15,:]
    lat_2d = id.variables['lat_rho'][:-15,:]
    id.close()

    # Throw away periodic boundary overlap
    data_3d = data_3d[:,:,:-2]
    zeta = zeta[:,:-2]
    h = h[:,:-2]
    zice = zice[:,:-2]
    lon_2d = lon_2d[:,:-2]
    lat_2d = lat_2d[:,:-2]        

    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    z_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N, zeta)

    # Choose what to write on the title about longitude
    if lon_key == 0:
        if lon0 < 0:
            lon_string = 'at ' + str(int(round(-lon0))) + r'$^{\circ}$W'
        else:
            lon_string = 'at ' + str(int(round(lon0))) + r'$^{\circ}$E'
    elif lon_key == 1:
        lon_string = 'zonally averaged'
    elif lon_key == 2:
        lon_string = 'zonally averaged between '
        if lon_bounds[0] < 0:
            lon_string += str(int(round(-lon_bounds[0]))) + r'$^{\circ}$W and '
        else:
            lon_string += str(int(round(lon_bounds[0]))) + r'$^{\circ}$E and '
        if lon_bounds[1] < 0:
            lon_string += str(int(round(-lon_bounds[1]))) + r'$^{\circ}$W'
        else:
            lon_string += str(int(round(lon_bounds[1]))) + r'$^{\circ}$E'

    # Edit longitude bounds to be from 0 to 360, to fit with ROMS convention
    if lon_key == 0:
        if lon0 < 0:
            lon0 += 360
    elif lon_key == 2:
        if lon_bounds[0] < 0:
            lon_bounds[0] += 360
        if lon_bounds[1] < 0:
            lon_bounds[1] += 360

    # Interpolate or average data
    if lon_key == 0:
        # Interpolate to lon0
        data, z, lat = interp_lon(data_3d, z_3d, lat_2d, lon_2d, lon0)
    elif lon_key == 1:
        # Zonally average over all longitudes
        # dlon is constant on this grid (0.25 degrees) so this is easy
        data = mean(data_3d, axis=2)
        z = mean(z_3d, axis=2)
        # Zonally average latitude, and copy into N depth levels
        lat = tile(mean(lat_2d, axis=1), (N,1))
    elif lon_key == 2:
        # Zonally average between lon_bounds
        data, z, lat = average_btw_lons(data_3d, z_3d, lat_2d, lon_2d, lon_bounds)

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
        if var_name in ['u', 'v']:
            # Center levels on 0 for certain variables, with a blue-to-red
            # colourmap
            max_val = amax(abs(data))
            lev = linspace(-max_val, max_val, num=40)
            colour_map = 'RdYlBu_r'
        else:
            lev = linspace(amin(data), amax(data), num=40)
            colour_map = 'jet'

    # Plot
    fig = figure(figsize=(18,6))
    contourf(lat, z, data, lev, cmap=colour_map, extend='both')
    colorbar()

    title(long_name + ' (' + units + ')\n' + lon_string)
    xlabel('Latitude')
    ylabel('Depth (m)')

    # Choose latitude bounds based on land mask
    data_sum = sum(data, axis=0)    
    # Find southernmost and northernmost unmasked j-indices
    edges = ma.flatnotmasked_edges(data_sum)
    j_min = edges[0]
    j_max = edges[1]
    if j_min == 0:
        # There are ocean points right to the southern boundary
        # Don't do anything special
        lat_min = min(lat[:,j_min])
    else:
        # There is land everywhere at the southern boundary
        # Show the last 2 degrees of this land mask
        lat_min = min(lat[:,j_min]) - 2
    if j_max == size(data_sum) - 1:
        # There are ocean points right to the northern boundary
        # Don't do anything special
        lat_max = max(lat[:,j_max])
    else:
        # There is land everywhere at the northern boundary
        # Show the first 2 degrees of this land mask
        lat_max = max(lat[:,j_max]) + 2
#    lat_max = -65
    xlim([lat_min, lat_max])
    ylim([depth_min, 0])

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()

    # Reset lon0 or lon_bounds to (-180, 180) range in case we
    # use them again for the next plot
    if lon_key == 0:
        if lon0 > 180:
            lon0 -= 360
    elif lon_key == 2:
        if lon_bounds[0] > 180:
            lon_bounds[0] -= 360
        if lon_bounds[1] > 180:
            lon_bounds[1] -= 360


# Linearly interpolate data, z, and latitude to the specified longitude.
# Input:
# data_3d = array of data, dimension depth x lat x lon
# z_3d = array of depth values (negative, in metres), dimension depth x lat x lon
# lat_2d = array of latitudevalues, dimension lat x lon
# lon_2d = array of longitude values, dimension lat x lon (between -180 and 180)
# lon0 = longitude to interpolate to (between -180 and 180)
# Output:
# data = array of data interpolated to lon0, dimension depth x lat
# z = array of depth values interpolated to lon0, dimension depth x lat
# lat = array of latitude values interpolated to lon0, dimension depth x lat
def interp_lon (data_3d, z_3d, lat_2d, lon_2d, lon0):

    # Save dimensions
    num_depth = size(data_3d, 0)
    num_lat = size(data_3d, 1)
    num_lon = size(data_3d, 2)
    # Set up output arrays
    data = ma.empty([num_depth, num_lat])
    z = ma.empty([num_depth, num_lat])
    lat = ma.empty([num_depth, num_lat])

    # Loop over latitudes; can't find a cleaner way to do this
    for j in range(num_lat):
        # Extract the longitude values of this slice
        lon_tmp = lon_2d[j,:]
        # Get indices and coefficients for interpolation
        ie, iw, coeffe, coeffw = interp_lon_helper(lon_tmp, lon0)        
        data[:,j] = coeffe*data_3d[:,j,ie] + coeffw*data_3d[:,j,iw]
        z[:,j] = coeffe*z_3d[:,j,ie] + coeffw*z_3d[:,j,iw]
        lat[:,j] = coeffe*lat_2d[j,ie] + coeffw*lat_2d[j,iw]

    return data, z, lat


# Zonally average data, z, and latitude between the specified longitude bounds.
# Input:
# data_3d = array of data, dimension depth x lat x lon
# z_3d = array of depth values (negative, in metres), dimension depth x lat x lon
# lat_2d = array of latitude values, dimension lat x lon
# lon_2d = array of longitude values, dimension lat x lon (between -180 and 180)
# lon_bounds = longitudes to average between, stored as an array of size 2 with
#              the western bound first
# Output:
# data = array of data averaged between lon_bounds, dimension depth x lat
# z = array of depth values averaged between lon_bounds, dimension depth x lat
# lat = array of latitude values averaged between lon_bounds, dimension 
#       depth x lat
def average_btw_lons (data_3d, z_3d, lat_2d, lon_2d, lon_bounds):

    dlon = 0.25 # Regular lon spacing: change this for new grids

    # Save horizontal dimensions
    num_depth = size(data_3d, 0)
    num_lat = size(data_3d, 1)
    num_lon = size(data_3d, 2)
    # Set up output arrays
    data = ma.empty([num_depth, num_lat])
    z = ma.empty([num_depth, num_lat])
    lat = ma.empty([num_depth, num_lat])
    # Unpack lon_bounds
    lon0_w = lon_bounds[0]
    lon0_e = lon_bounds[1]

    # Loop over latitudes; can't find a cleaner way to do this
    for j in range(num_lat):

        # Extract the longitude values of this slice
        lon_tmp = lon_2d[j,:]
        # Initialise integrals for data, z, and lon
        data_int = zeros([num_depth])
        z_int = zeros([num_depth])
        lat_int = 0.0
        lon_int = 0.0

        # Linearly interpolate data, z, and lat to lon0_w
        ie, iw, coeffe, coeffw = interp_lon_helper(lon_tmp, lon0_w)
        data_w = coeffe*data_3d[:,j,ie] + coeffw*data_3d[:,j,iw]
        z_w = coeffe*z_3d[:,j,ie] + coeffw*z_3d[:,j,iw]
        lat_w = coeffe*lat_2d[j,ie] + coeffw*lat_2d[j,iw]

        # Integrate data, z, lat, and lon between lon0_w and the nearest point to the east
        if all(data_3d[:,j,ie].mask):
            # We're on a land point; skip this integration
            pass
        else:
            dlon_w = lon_tmp[ie] - lon0_w
            # Take mod 360 if needed
            if dlon_w > 300:
                dlon_w -= 360
            elif dlon_w < -60:
                dlon_w += 360
            data_int += 0.5*(data_w + data_3d[:,j,ie])*dlon_w # Trapezoidal rule
            z_int += 0.5*(z_w + z_3d[:,j,ie])*dlon_w
            lat_int += 0.5*(lat_w + lat_2d[j,ie])*dlon_w
            lon_int += dlon_w
        # Save the index to start normal (non-interpolated) integration later
        istart = ie

        # Linearly interpolate data, z, and lat to lon0_e
        ie, iw, coeffe, coeffw = interp_lon_helper(lon_tmp, lon0_e)
        data_e = coeffe*data_3d[:,j,ie] + coeffw*data_3d[:,j,iw]
        z_e = coeffe*z_3d[:,j,ie] + coeffw*z_3d[:,j,iw]
        lat_e = coeffe*lat_2d[j,ie] + coeffw*lat_2d[j,iw]

        # Integrate data, z, lat, and lon between lon0_e and the nearest point to the west (iw)
        if all(data_3d[:,j,iw].mask):
            # We're on a land point; skip this integration
            pass
        else:
            dlon_e = lon0_e - lon_tmp[iw]
            # Take mod 360 if needed
            if dlon_e > 300:
                dlon_e -= 360
            elif dlon_e < -60:
                dlon_e += 360
            data_int += 0.5*(data_3d[:,j,iw] + data_e)*dlon_e # Trapezoidal rule
            z_int += 0.5*(z_3d[:,j,iw] + z_e)*dlon_e
            lat_int += 0.5*(lat_2d[j,iw] + lat_e)*dlon_e
            lon_int += dlon_e
        # Save the index which normal (non-interpolated) integration should stop before
        iend = ie

        # Normal integration between istart and iend
        if istart <= iend:
            # Normal case
            i_list = range(istart, iend)
        elif istart > iend:
            # lon_bounds crosses the periodic boundary
            i_list = range(istart, num_lon) + range(0, iend)
        # Have to loop over i values because of masking checks
        for i in i_list:
            if all(data_3d[:,j,i].mask):
                # We're on a land point; skip this integration
                pass
            else:
                data_int += data_3d[:,j,i]*dlon
                z_int += z_3d[:,j,i]*dlon
                lat_int += lat_2d[j,i]*dlon
                lon_int += dlon

        if lon_int == 0:
            # There weren't any ocean points between lon_bounds
            data[:,j] = ma.masked
            z[:,j] = ma.masked
            lat[:,j] = ma.masked
        else:
            # Divide integrals of data, z, and lat by lon_int to get averages
            data[:,j] = data_int/lon_int
            z[:,j] = z_int/lon_int
            lat[:,j] = lat_int/lon_int

    return data, z, lat


# Calculate indices and coefficients for linear interpolation of longitude.
# This takes care of all the mod 360 nonsense.
# Input:
# lon = 1D array of longitude values (straight out of ROMS i.e. between slightly < 0 and slightly > 360)
# lon0 = longitude to interpolate to (between 0 and 360)
# Output:
# ie, iw, coeffe, coeffw = integers (ie and iw) and coefficients (coeffe and coeffw) such that
#                          coeffe*lon[ie] + coeffw*lon[iw] = lon0, which will also hold for any
#                          variable on this longitude grid. ie is the index of the nearest point
#                          to the east of lon0; iw the nearest point to the west.
def interp_lon_helper (lon, lon0):

    if lon0 < amin(lon) or lon0 > amax(lon):
        # Special case: lon0 on periodic boundary
        # Be careful with mod 360 here

        # Find the periodic boundary
        dlon = lon[1:] - lon[0:-1]
        bdry = argmax(abs(dlon))
        if dlon[bdry] < -300:
            # Jumps from almost 360 to just over 0
            iw = bdry
            ie = bdry + 1
        else:
            # Periodic boundary lines up with the array boundary
            iw = size(lon) - 1
            ie = 0
        # Calculate difference between lon0 and lon[iw], mod 360 if necessary
        dlon_num = lon0 - lon[iw]
        if dlon_num < -300:
            dlon_num += 360
        # Calculate difference between lon[ie] and lon[iw], mod 360
        dlon_den = lon[ie] - lon[iw] + 360

    else:
        # General case

        # Add or subtract 360 from longitude values which wrap around
        # so that longitude increases monotonically from west to east
        i = arange(1, size(lon)+1)
        index1 = nonzero((i > 1200)*(lon < 100))
        lon[index1] = lon[index1] + 360
        index2 = nonzero((i < 200)*(lon > 300))
        lon[index2] = lon[index2] - 360

        # Take mod 360 of lon0 if necessary
        if all(lon < lon0):
            lon0 -= 360
        if all(lon > lon0):
            lon0 += 360
        
        # Find the first index eastward of lon0
        ie = nonzero(lon > lon0)[0][0]
        # The index before it will be the last index westward of lon0
        iw = ie - 1

        dlon_num = lon0 - lon[iw]
        dlon_den = lon[ie] - lon[iw]

    if dlon_num > 5 or dlon_den > 5:
        print 'interp_lon_helper: Problem at periodic boundary'
        return
    coeff1 = dlon_num/dlon_den
    coeff2 = 1 - coeff1

    return ie, iw, coeff1, coeff2


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to ocean history/averages file: ")
    var_name = raw_input("Variable name: ")
    if var_name in ['u', 'v']:
        # Will need the grid file to get the angle
        grid_path = raw_input("Path to ROMS grid file: ")
    else:
        grid_path = None
    tstep = int(raw_input("Timestep number (starting at 1): "))

    lon_type = raw_input("Single longitude (s) or zonal average (z)? ")
    if lon_type == 's':
        lon_key = 0
        lon0 = float(raw_input("Enter longitude (-180 to 180): "))
        lon_bounds = None
    elif lon_type == 'z':
        avg_type = raw_input("Zonal average over all longitudes (a) or between two specific longitudes (s)? ")
        if avg_type == 'a':
            lon_key = 1
            lon0 = NaN
            lon_bounds = None
        elif avg_type == 's':
            lon_key = 2
            lon0 = NaN
            w_bound = float(raw_input("Enter western bound on longitude (-180 to 180): "))
            e_bound = float(raw_input("Enter eastern bound on longitude (-180 to 180): "))
            lon_bounds = [w_bound, e_bound]

    depth_min = -1*float(raw_input("Deepest depth to plot (positive, metres): "))

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

    zonal_plot(file_path, var_name, tstep, lon_key, lon0, lon_bounds, depth_min, colour_bounds, save, fig_name, grid_path)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the input parameters; repeat until the user is finished
                changes = raw_input("Enter a parameter to change: (1) file path, (2) variable name, (3) timestep number, (4) longitude, (5) deepest depth, (6) colour bounds, (7) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New file path
                        file_path = raw_input("Path to ocean history/averages file: ")
                    elif int(changes) == 2:
                        # New variable name
                        var_name = raw_input("Variable name: ")
                        if var_name in ['u', 'v'] and grid_path is None:
                            # Will need the grid file to get the angle
                            grid_path = raw_input("Path to ROMS grid file: ")
                    elif int(changes) == 3:
                        # New timestep
                        tstep = int(raw_input("Timestep number (starting at 1): "))
                    elif int(changes) == 4:
                        # New longitude information
                        lon_type = raw_input("Single longitude (s) or zonal average (z)? ")
                        if lon_type == 's':
                            lon_key = 0
                            lon0 = float(raw_input("Enter longitude (-180 to 180): "))
                            lon_bounds = None
                        elif lon_type == 'z':
                            avg_type = raw_input("Zonal average over all longitudes (a) or between two specific longitudes (s)? ")
                            if avg_type == 'a':
                                lon_key = 1
                                lon0 = NaN
                                lon_bounds = None
                            elif avg_type == 's':
                                lon_key = 2
                                lon0 = NaN
                                w_bound = float(raw_input("Enter western bound on longitude (-180 to 180): "))
                                e_bound = float(raw_input("Enter eastern bound on longitude (-180 to 180): "))
                                lon_bounds = [w_bound, e_bound]
                    elif int(changes) == 5:
                        # New depth bound
                        depth_min = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
                    elif int(changes) == 6:
                        # Get colour bounds if necessary
                        colour_bounds = None
                        get_bounds = raw_input("Set bounds on colour scale (y/n)? ")
                        if get_bounds == 'y':
                            lower_bound = float(raw_input("Lower bound: "))
                            upper_bound = float(raw_input("Upper bound: "))
                            colour_bounds = [lower_bound, upper_bound]
                    elif int(changes) == 7:
                        # Change from display to save, or vice versa
                        save = not save
            if save:
                # Get file name for figure
                fig_name = raw_input("File name for figure: ")

            # Make the plot
            zonal_plot(file_path, var_name, tstep, lon_key, lon0, lon_bounds, depth_min, colour_bounds, save, fig_name, grid_path)

        else:
            break
                        

    

    

    
