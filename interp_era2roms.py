from numpy import *
from scipy.interpolate import LinearNDInterpolator, RectBivariateSpline

# Given an array on the ERA-Interim grid, interpolate any missing values, and
# then interpolate to the ROMS grid.
# Input:
# A = array of size nxm containing values on the ERA-Interim grid (dimension
#     longitude x latitude)
# lon_era = array of length n containing ERA-Interim longitude values
# lat_era = array of length m containing ERA-Interim latitude values
# lon_roms = array of size pxq containing ROMS longitude values
# lat_roms = array of size pxq containing ROMS latitude values
# Output:
# B = array of size pxq containing values on the ROMS grid (dimension latitude x
#     longitude)
def interp_era2roms (A, lon_era, lat_era, lon_roms, lat_roms):

    # Save the sizes of ROMS axes
    num_lon = size(lon_roms, 1)
    num_lat = size(lon_roms, 0)

    # Missing values are something <<0, but these can change when the offset
    # and scale factor attributes are automatically applied. Either way,
    # missing values will be the minimum values in A.
    flag = amin(A)

    # Interpolate missing values
    # I got this bit of code from Stack Exchange
    # It seems to work, not exactly sure how
    valid_mask = A > flag
    coords = array(nonzero(valid_mask)).T
    values = A[valid_mask]
    fill_function = LinearNDInterpolator(coords, values)
    Afill = fill_function(list(ndindex(A.shape))).reshape(A.shape)
    # Fill any still-missing values with the mean
    Afill[isnan(Afill)] = nanmean(Afill)

    # Now interpolate from ERA-Interim grid to ROMS grid
    # First build a function to approximate A with 2D splines
    # Note that latitude has to be multiplied by -1 so both axes are strictly
    # ascending
    interp_function = RectBivariateSpline(lon_era, -lat_era, Afill)
    B = zeros(shape(lon_roms))
    # Call this function for each grid point individually - if you try to do
    # it all at once it throws a MemoryError
    for i in range(num_lon):
        for j in range(num_lat):
            B[j,i] = interp_function(lon_roms[j,i], -lat_roms[j,i])

    # Enforce periodic boundary
    B[:,0] = B[:,-2]
    B[:,-1] = B[:,1]

    return B
