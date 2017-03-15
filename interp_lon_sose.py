from numpy import *

# Linearly interpolate SOSE data to the specified longitude.
# Input:
# data_3d = arary of data, dimension depth x lat x lon
# lon = 1D array of longitude values (between 0 and 360)
# lon0 = longitude to interpolate to (between 0 and 360)
# Output:
# data = array of data interpolated to lon0, dimension depth x lat
def interp_lon_sose (data_3d, lon, lon0):

    if lon0 < lon[0] or lon0 > lon[-1]:
        # Special case: lon0 on periodic boundary
        # Be careful with mod 360 here
        iw = size(lon)-1
        ie = 0
        # Calculate difference between lon0 and lon[iw], mod 360 if necessary
        dlon_num = lon0 - lon[iw]
        if dlon_num < -300:
            dlon_num += 360
        # Calculate difference between lon[ie] and lon[iw], mod 360
        dlon_den = lon[ie] - lon[iw] + 360
    else:
        # General case
        # Find the first index eastwards of lon0
        ie = nonzero(lon > lon0)[0][0]
        # The index before it will be the last index westward of lon0
        iw = ie - 1
        dlon_num = lon0 - lon[iw]
        dlon_den = lon[ie] - lon[iw]
    # Coefficients for interpolation
    coeff1 = dlon_num/dlon_den
    coeff2 = 1 - coeff1

    # Interpolate
    data = coeff1*data_3d[:,:,ie] + coeff2*data_3d[:,:,iw]
    return data
