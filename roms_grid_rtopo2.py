from netCDF4 import Dataset
from numpy import *
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage.filters import gaussian_filter

# Interpolate RTopo-2 data to an existing ROMS grid, and do some smoothing. 
# Then overwrite the bathymetry, ice shelf draft, and masks in the ROMS grid
# file with the new RTopo-2 values.

# Main routine
def run (rtopo_data, rtopo_aux, roms_grid, min_h, max_h, min_zice, h_smooth, zice_smooth):

    # Read RTopo-2 data
    id = Dataset(rtopo_data, 'r')
    lon_rtopo = id.variables['lon'][:]
    lat_rtopo = id.variables['lat'][:]
    # Make bathymetry positive for ROMS conventions
    h_rtopo = -1*transpose(id.variables['bathy'][:,:])
    zice_rtopo = transpose(id.variables['ice_bottom'][:,:])
    id.close()

    # Read RTopo-2 mask
    id = Dataset(rtopo_aux, 'r')
    mask = transpose(id.variables['amask'][:,:])
    id.close()
    # amask has 4 values: 0 = open ocean, 1 = grounded ice, 2 = ice shelves,
    # 3 = bare rock
    ocean_flag = 0.0
    iceshelf_flag = 2.0
    # Make a new mask which is 1 at all ocean points (open ocean or ice shelves)
    # and 0 at all land points (ice sheet or rock)
    mask_ocn_rtopo = zeros(shape(mask))    
    mask_ocn_rtopo[abs(mask-ocean_flag) < 0.1] = 1.0
    mask_ocn_rtopo[abs(mask-iceshelf_flag) < 0.1] = 1.0
    # Make another mask which is 1 at all ice shelf points and 0 elsewhere
    mask_zice_rtopo = zeros(shape(mask))
    mask_zice_rtopo[abs(mask-iceshelf_flag) < 0.1] = 1.0
    # At land points, set bathymetry to min_h
    h_rtopo[mask_ocn_rtopo < 0.1] = min_h
    # At non-ice shelf points, set ice shelf draft to 0
    zice_rtopo[mask_zice_rtopo < 0.1] = 0.0

    # Read ROMS latitude and longitude on all grids (rho, u, v, psi)
    id = Dataset(roms_grid, 'r')
    lon_rho = id.variables['lon_rho'][:,:]
    lat_rho = id.variables['lat_rho'][:,:]
    lon_u = id.variables['lon_u'][:,:]
    lat_u = id.variables['lat_u'][:,:]
    lon_v = id.variables['lon_v'][:,:]
    lat_v = id.variables['lat_v'][:,:]
    lon_psi = id.variables['lon_psi'][:,:]
    lat_psi = id.variables['lat_psi'][:,:]
    id.close()

    # ROMS longitude goes from 0 to 360; make it go from -180 to 180 instead
    # (taking mod 360 of values outside this range) so it matches RTopo-2
    index = lon_rho > 180
    lon_rho[index] = lon_rho[index] - 360
    index = lon_u > 180
    lon_u[index] = lon_u[index] - 360
    index = lon_v > 180
    lon_v[index] = lon_v[index] - 360
    index = lon_psi > 180
    lon_psi[index] = lon_psi[index] - 360

    # Interpolate bathymetry, ice shelf draft, and ice shelf mask to rho grid
    h = interp_rtopo2roms(h_rtopo, lon_rtopo, lat_rtopo, lon_rho, lat_rho)
    zice = interp_rtopo2roms(zice_rtopo, lon_rtopo, lat_rtopo, lon_rho, lat_rho)
    mask_zice = interp_rtopo2roms(mask_zice_rtopo, lon_rtopo, lat_rtopo, lon_rho, lat_rho)
    # Interpolate land-ocean mask to each ROMS grid
    mask_rho = interp_rtopo2roms(mask_ocn_rtopo, lon_rtopo, lat_rtopo, lon_rho, lat_rho)
    mask_u = interp_rtopo2roms(mask_ocn_rtopo, lon_rtopo, lat_rtopo, lon_u, lat_u)
    mask_v = interp_rtopo2roms(mask_ocn_rtopo, lon_rtopo, lat_rtopo, lon_v, lat_v)
    mask_psi = interp_rtopo2roms(mask_ocn_rtopo, lon_rtopo, lat_rtopo, lon_psi, lat_psi)

    # Interpolation of masks near the coast will lead to values between 0 and 1.
    # Divide these values into 0s and 1s with the cutoff point 0.1 (i.e. coastal
    # points are only considered ocean if they are at least 90% ocean in the
    # interpolation)
    mask_zice[mask_zice < 0.1] = 0.0
    mask_zice[mask_zice >= 0.1] = 1.0
    mask_rho[mask_rho < 0.1] = 0.0
    mask_rho[mask_rho >= 0.1] = 1.0
    mask_u[mask_u < 0.1] = 0.0
    mask_u[mask_u >= 0.1] = 1.0
    mask_v[mask_v < 0.1] = 0.0
    mask_v[mask_v >= 0.1] = 1.0
    mask_psi[mask_psi < 0.1] = 0.0
    mask_psi[mask_psi >= 0.1] = 1.0

    # RTopo-2 has a stranded bit of water beneath an ice shelf in East
    # Antarctica, cut off from the ocean. ROMS won't like this and we don't
    # really care about it, so fill it in with land.
    mask_zice[0:50,400:600] = 0.0
    mask_rho[0:50,400:600] = 0.0
    mask_u[0:50,400:600] = 0.0
    mask_v[0:50,400:600] = 0.0
    mask_psi[0:50,400:600] = 0.0

    # At land points, set bathymetry to min_h again, in case the interpolation
    # changed this near the coast
    h[mask_rho < 0.1] = min_h
    # Similarly with zice, as before
    zice[mask_zice < 0.1] = 0.0
    # Apply clipping depths to bathymetry and zice
    h[h < min_h] = min_h
    h[h > max_h] = max_h    
    zice[-zice < min_zice] = -min_zice

    # Smooth h and zice
    h = gaussian_filter(h, h_smooth, mode='nearest')
    zice = gaussian_filter(zice, zice_smooth, mode='nearest')

    # Overwrite h, zice, and masks in the ROMS grid file
    id = Dataset(roms_grid, 'a')
    id.variables['h'][:,:] = h
    id.variables['zice'][:,:] = zice
    id.variables['mask_zice'][:,:] = mask_zice
    id.variables['mask_rho'][:,:] = mask_rho
    id.variables['mask_u'][:,:] = mask_u
    id.variables['mask_v'][:,:] = mask_v
    id.variables['mask_psi'][:,:] = mask_psi
    id.close()   

# Given an array on the RTopo-2 grid, interpolate to the ROMS grid.
# Input:
# A = array of size nxm containing values on the RTopo-2 grid (dimension
#     longitude x latitude)
# lon_rtopo = array of length n containing RTopo-2 longitude values
# lat_rtopo = array of length m containing RTopo-2 latitude values
# lon_roms = array of size pxq containing ROMS longitude values (converted to
#            the range -180 to 180, as shown in the main routine)
# lat_roms = array of size pxq containing ROMS latitude values
# Output:
# B = array of size pxq containing values on the ROMS grid (dimension latitude x
#     longitude)
def interp_rtopo2roms (A, lon_rtopo, lat_rtopo, lon_roms, lat_roms):

    # Build a function for linear interpolation of A
    # Setting bounds_error=False and fill_value=None allows it to extrapolate
    # if necessary; this will only happen at the northernmost 0.5 degrees
    interp_function = RegularGridInterpolator((lon_rtopo, lat_rtopo), A, bounds_error=False, fill_value=None)
    # Call this function at the ROMS grid points
    B = interp_function((lon_roms, lat_roms))
    return B


# User parameters
if __name__ == "__main__":

    # Path to RTopo-2 file with variables lon, lat, bathy, ice_bottom
    rtopo_data = 'RTopo_2.0_30sec_Antarctica_data__2015-11-02.nc'
    # Path to RTopo-2 file with variable amask
    rtopo_aux = 'RTopo_2.0_30sec_Antarctica_aux__2015-11-02.nc'
    # Path to existing ROMS grid file to edit
    roms_grid = 'rtopo2_circumpolar_quarterdegree.nc'
    # Minimum and maximum seafloor depth
    min_h = 50.0
    max_h = 6000.0
    # Minimum ice shelf draft
    min_zice = 50.0
    # Standard deviations of Gaussian filters for smoothing bathymetry (h)
    # and ice shelf draft (zice); higher values mean more smoothing
    h_smooth = 0.5
    zice_smooth = 0.4

    run(rtopo_data, rtopo_aux, roms_grid, min_h, max_h, min_zice, h_smooth, zice_smooth)




