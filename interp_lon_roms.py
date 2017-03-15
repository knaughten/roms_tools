from numpy import *

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
def interp_lon_roms (data_3d, z_3d, lat_2d, lon_2d, lon0):

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
# ie, iw, coeffe, coeffw = integers (ie and iw) and coefficients 
#                          (coeffe and coeffw) such that 
#                          coeffe*lon[ie] + coeffw*lon[iw] = lon0, 
#                          which will also hold for any variable on this 
#                          longitude grid. ie is the index of the nearest 
#                          point to the east of lon0; iw the nearest point 
#                          to the west.
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
