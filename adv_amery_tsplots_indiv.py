from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from calc_z import *

# For each advection experiment, plot zonal slices of temperature and salinity
# through 71E (Amery Ice Shelf) at the end of the simulation.
def adv_amery_tsplots_indiv ():

    num_simulations = 6
    # Paths to simulation directories
    paths = ['/short/m68/kaa561/ROMS-CICE-MCT/tmproms/run/advection/c4_lowdif/', '/short/m68/kaa561/ROMS-CICE-MCT/tmproms/run/advection/c4_highdif/', '/short/m68/kaa561/ROMS-CICE-MCT/tmproms/run/advection/a4_lowdif/', '/short/m68/kaa561/ROMS-CICE-MCT/tmproms/run/advection/a4_highdif/', '/short/m68/kaa561/ROMS-CICE-MCT/tmproms/run/advection/u3_lowdif/', '/short/m68/kaa561/ROMS-CICE-MCT/tmproms/run/advection/u3limiters_lowdif/']
    # End of figure names for each simulation
    labels = ['_c4_lowdif.png', '_c4_highdif.png', '_a4_lowdif.png', '_a4_highdif.png', '_u3.png', '_u3_lim.png']
    # Name of ocean output file to read
    ocn_file = 'ocean_avg_0001.nc'
    # Timestep to plot (average over last day in 1992)
    tstep = 366
    # Longitude to plot
    lon0 = 71
    # Deepest depth to plot
    depth_min = -500
    # Bounds on colour scale for each variable
    temp_bounds = [-2, 3]
    salt_bounds = [33.8, 34.8]
    # Bounds on latitudes to plot
    lat_min = -72
    lat_max = -50

    # Grid parameters
    theta_s = 4.0
    theta_b = 0.9
    hc = 40
    N = 31

    # Build titles for each variable based on longitude
    if lon0 < 0:
        temp_title = r'Temperature ($^{\circ}$C) at ' + str(int(round(-lon0))) + r'$^{\circ}$W'
        salt_title = r'Salinity (psu) at ' + str(int(round(-lon0))) + r'$^{\circ}$W'
        # Edit longitude to be between0 and 360, following ROMS convention
        lon0 += 360
    else:
        temp_title = r'Temperature ($^{\circ}$C) at ' + str(int(round(lon0))) + r'$^{\circ}$E'
        salt_title = r'Salinity (psu) at ' + str(int(round(lon0))) + r'$^{\circ}$E'

    # Loop over simulations
    for sim in range(num_simulations):
        # Loop over variables
        for var_name in ['temp', 'salt']:
            # Read variable, sea surface height, and grid variables
            id = Dataset(paths[sim] + ocn_file, 'r')
            data_3d = id.variables[var_name][tstep-1,:,:-15,:]
            zeta = id.variables['zeta'][tstep-1,:-15,:]
            if sim == 0 and var_name == 'temp':
                # Grid variables are the same for all simulations so we
                # only need to read them once
                h = id.variables['h'][:-15,:]
                zice = id.variables['zice'][:-15,:]
                lon_2d = id.variables['lon_rho'][:-15,:]
                lat_2d = id.variables['lat_rho'][:-15,:]
            id.close()
            # Get a 3D array of z-coordinates
            z_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N, zeta)
            # Interpolate the variable, z, and latitude to lon0
            data, z, lat = interp_lon(data_3d, z_3d, lat_2d, lon_2d, lon0)
            # Set up colour levels for plotting
            if var_name == 'temp':
                lev = linspace(temp_bounds[0], temp_bounds[1], num=40)
            elif var_name == 'salt':
                lev = linspace(salt_bounds[0], salt_bounds[1], num=40)
            # Plot
            fig = figure(figsize=(12,6))
            contourf(lat, z, data, lev, cmap='jet', extend='both')
            colorbar()
            if var_name == 'temp':
                title(temp_title)
            elif var_name == 'salt':
                title(salt_title)
            xlabel('Latitude')
            ylabel('Depth (m)')
            xlim([lat_min, lat_max])
            ylim([depth_min, 0])
            # Save plot
            fig.savefig(var_name + labels[sim])


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

    adv_amery_tsplots_indiv()
