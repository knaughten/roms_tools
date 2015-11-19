# Build an initialisation file for ROMS using climatology temperature and
# salinity values from the World Ocean Atlas. Set initial velocities and sea
# surface height to zero.

# NB for raijin users: RegularGridInterpolator needs python/2.7.6 but the
# default is 2.7.3. Before running this script, switch them as follows:
# module unload python/2.7.3
# module unload python/2.7.3-matplotlib
# module load python/2.7.6
# module load python/2.7.6-matplotlib

from netCDF4 import Dataset
from numpy import *
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import KDTree
from calc_z import *


# Main routine
def run (grid_file, woa_file, output_file, Tcline, theta_s, theta_b, hc, N, nbdry_woa):

    # Read WOA data and grid
    print 'Reading World Ocean Atlas data'
    woa_id = Dataset(woa_file, 'r')
    lon_woa = woa_id.variables['longitude'][:]
    lat_woa = woa_id.variables['latitude'][:nbdry_woa]
    depth_woa = woa_id.variables['depth'][:]
    temp_woa = transpose(woa_id.variables['temp'][:,:nbdry_woa,:])
    salt_woa = transpose(woa_id.variables['salt'][:,:nbdry_woa,:])
    woa_id.close()

    # Read ROMS grid
    print 'Reading ROMS grid'
    grid_id = Dataset(grid_file, 'r')
    lon_roms = grid_id.variables['lon_rho'][:,:]
    lat_roms = grid_id.variables['lat_rho'][:,:]
    h = grid_id.variables['h'][:,:]
    zice = grid_id.variables['zice'][:,:]
    mask_rho = grid_id.variables['mask_rho'][:,:]
    mask_zice = grid_id.variables['mask_zice'][:,:]
    grid_id.close()
    num_lon = size(lon_roms, 1)
    num_lat = size(lon_roms, 0)
    # Mask h and zice with zeros
    h = h*mask_rho
    zice = zice*mask_zice

    # Get 3D array of ROMS z-coordinates, as well as 1D arrays of s-coordinates
    # and stretching curves
    z_roms_3d, sc_r, Cs_r = calc_z(h, zice, lon_roms, lat_roms, theta_s, theta_b, hc, N)
    # Copy the latitude and longitude values into 3D arrays of the same shape
    lon_roms_3d = tile(lon_roms, (N,1,1))
    lat_roms_3d = tile(lat_roms, (N,1,1))

    # Regridding happens here
    print 'Interpolating temperature'
    temp = interp_woa2roms(temp_woa, lon_woa, lat_woa, depth_woa, lon_roms_3d, lat_roms_3d, z_roms_3d, mask_rho, mask_zice, -0.5)
    print 'Interpolating salinity'
    salt = interp_woa2roms(salt_woa, lon_woa, lat_woa, depth_woa, lon_roms_3d, lat_roms_3d, z_roms_3d, mask_rho, mask_zice, 34.5)

    # Set initial velocities and sea surface height to zero
    u = zeros((N, num_lat, num_lon-1))
    v = zeros((N, num_lat-1, num_lon))
    ubar = zeros((num_lat, num_lon-1))
    vbar = zeros((num_lat-1, num_lon))
    zeta = zeros((num_lat, num_lon))

    print 'Writing to NetCDF file'
    out_id = Dataset(output_file, 'w')
    # Define dimensions
    out_id.createDimension('xi_u', num_lon-1)
    out_id.createDimension('xi_v', num_lon)
    out_id.createDimension('xi_rho', num_lon)
    out_id.createDimension('eta_u', num_lat)
    out_id.createDimension('eta_v', num_lat-1)
    out_id.createDimension('eta_rho', num_lat)
    out_id.createDimension('s_rho', N)
    out_id.createDimension('ocean_time', None)
    out_id.createDimension('one', 1);
    # Define variables and assign values
    out_id.createVariable('tstart', 'f8', ('one'))
    out_id.variables['tstart'].long_name = 'start processing day'
    out_id.variables['tstart'].units = 'day'
    out_id.variables['tstart'][:] = 0.0
    out_id.createVariable('tend', 'f8', ('one'))
    out_id.variables['tend'].long_name = 'end processing day'
    out_id.variables['tend'].units = 'day'
    out_id.variables['tend'][:] = 0.0
    out_id.createVariable('theta_s', 'f8', ('one'))
    out_id.variables['theta_s'].long_name = 'S-coordinate surface control parameter'
    out_id.variables['theta_s'][:] = theta_s
    out_id.createVariable('theta_b', 'f8', ('one'))
    out_id.variables['theta_b'].long_name = 'S-coordinate bottom control parameter'
    out_id.variables['theta_b'].units = 'nondimensional'
    out_id.variables['theta_b'][:] = theta_b
    out_id.createVariable('Tcline', 'f8', ('one'))
    out_id.variables['Tcline'].long_name = 'S-coordinate surface/bottom layer width'
    out_id.variables['Tcline'].units = 'meter'
    out_id.variables['Tcline'][:] = Tcline
    out_id.createVariable('hc', 'f8', ('one'))
    out_id.variables['hc'].long_name = 'S-coordinate parameter, critical depth'
    out_id.variables['hc'].units = 'meter'
    out_id.variables['hc'][:] = hc
    out_id.createVariable('Cs_r', 'f8', ('s_rho'))
    out_id.variables['Cs_r'].long_name = 'S-coordinate stretching curves at RHO-points'
    out_id.variables['Cs_r'].units = 'nondimensional'
    out_id.variables['Cs_r'].valid_min = -1.0
    out_id.variables['Cs_r'].valid_max = 0.0
    out_id.variables['Cs_r'][:] = Cs_r
    out_id.createVariable('ocean_time', 'f8', ('ocean_time'))
    out_id.variables['ocean_time'].long_name = 'time since initialization'
    out_id.variables['ocean_time'].units = 'seconds'
    out_id.variables['ocean_time'][0] = 0.0
    out_id.createVariable('u', 'f8', ('ocean_time', 's_rho', 'eta_u', 'xi_u'))
    out_id.variables['u'].long_name = 'u-momentum component'
    out_id.variables['u'].units = 'meter second-1'
    out_id.variables['u'][0,:,:,:] = u
    out_id.createVariable('v', 'f8', ('ocean_time', 's_rho', 'eta_v', 'xi_v'))
    out_id.variables['v'].long_name = 'v-momentum component'
    out_id.variables['v'].units = 'meter second-1'
    out_id.variables['v'][0,:,:,:] = v
    out_id.createVariable('ubar', 'f8', ('ocean_time', 'eta_u', 'xi_u'))
    out_id.variables['ubar'].long_name = 'vertically integrated u-momentum component'
    out_id.variables['ubar'].units = 'meter second-1'
    out_id.variables['ubar'][0,:,:] = ubar
    out_id.createVariable('vbar', 'f8', ('ocean_time', 'eta_v', 'xi_v'))
    out_id.variables['vbar'].long_name = 'vertically integrated v-momentum component'
    out_id.variables['vbar'].units = 'meter second-1'
    out_id.variables['vbar'][0,:,:] = vbar
    out_id.createVariable('zeta', 'f8', ('ocean_time', 'eta_rho', 'xi_rho'))
    out_id.variables['zeta'].long_name = 'free-surface'
    out_id.variables['zeta'].units = 'meter'
    out_id.variables['zeta'][0,:,:] = zeta
    out_id.createVariable('temp', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'))
    out_id.variables['temp'].long_name = 'potential temperature'
    out_id.variables['temp'].units = 'Celsius'
    out_id.variables['temp'][0,:,:,:] = temp
    out_id.createVariable('salt', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'))
    out_id.variables['salt'].long_name = 'salinity'
    out_id.variables['salt'].units = 'PSU'
    out_id.variables['salt'][0,:,:,:] = salt
    out_id.createVariable('sc_r', 'f8', ('s_rho'))
    out_id.variables['sc_r'].long_name = 'S-coordinate at rho-points'
    out_id.variables['sc_r'].units = 'nondimensional'
    out_id.variables['sc_r'].valid_min = -1.0
    out_id.variables['sc_r'].valid_max = 0.0
    out_id.variables['sc_r'][:] = sc_r
    out_id.close()


# Given an array on the World Ocean Atlas grid, interpolate onto the ROMS grid,
# fill the ice shelf cavities with nearest-neighbour values from the same depth
# level, and fill the land mask with constant values.
# Input:
# A = array of size nxmxo containing values on the World Ocean Atlas grid
#     (dimension longitude x latitude x depth)
# lon_woa = array of length n containing WOA longitude values
# lat_woa = array of length m containing WOA latitude values
# depth_woa = array of length o containing WOA depth values in positive
#             z-coordinates
# lon_roms_3d = array of size pxqxr containing ROMS longitude values
# lat_roms_3d = array of size pxqxr containing ROMS latitude values
# z_roms_3d = array of size pxqxr containing ROMS depth values in negative
#             z-coordinates (converted from grid file using calc_z.py, as
#             shown below in the main script)
# mask_rho = array of size qxr containing ROMS land mask (ocean/ice shelf 1,
#            land 0)
# mask_zice = array of size qxr containing ROMS ice shelf mask (ice shelf 1,
#             ocean/land 0)
# land_fill = scalar containing the value with which to fill the ROMS land mask
#             (really doesn't matter, don't use NaN though)
# Output:
# B = array of size pxqxr containing values on the ROMS grid (dimension depth x
#     latitude x longitude)
def interp_woa2roms (A, lon_woa, lat_woa, depth_woa, lon_roms_3d, lat_roms_3d, z_roms_3d, mask_rho, mask_zice, land_fill):

    # Calculate N based on size of ROMS grid
    N = size(lon_roms_3d, 0)

    # Build a function for linear interpolation of A; set flag to fill
    # out-of-bounds values with NaN
    interp_function = RegularGridInterpolator((lon_woa, lat_woa, depth_woa), A, bounds_error=False, fill_value=NaN)
    B = zeros(shape(lon_roms_3d))

    # Call this function once at each depth level - 3D vectorisation uses too
    # much memory!
    for k in range(N):

        print '...vertical level ', str(k+1), ' of ', str(N)
        tmp = interp_function((lon_roms_3d[k,:,:], lat_roms_3d[k,:,:], -z_roms_3d[k,:,:]))

        # Fill missing values in ice shelf cavities with nearest-neighbour
        # values at this depth level
        # I got this code from Stack Exchange and am not entirely sure
        # how it works
        j,i = mgrid[0:tmp.shape[0], 0:tmp.shape[1]]
        # Good indices are not NaN and not land points
        good_index = nonzero(~isnan(tmp)*(mask_rho > 0.5))
        # Bad indicies are NaN and are ice shelf cavity points
        bad_index = nonzero(isnan(tmp)*(mask_zice > 0.5))
        jigood = array((j[good_index], i[good_index])).T
        jibad = array((j[bad_index], i[bad_index])).T
        # Magic happens here
        tmp[bad_index] = tmp[good_index][KDTree(jigood).query(jibad)[1]]

        # Fill land mask with constant value
        tmp[mask_rho==0] = land_fill
        # Save this depth level
        B[k,:,:] = tmp

    return B


# User parameters
if __name__ == "__main__":

    # Path to ROMS grid file
    grid_file = '../apps/common/grid/rtopo2_circumpolar_quarterdegree.nc'
    # Path to World Ocean Atlas NetCDF file (converted from FESOM input using 
    # woa_netcdf.py)
    woa_file = '/short/y99/kaa561/FESOM/woa01_ts.nc'
    # Path to desired output file
    output_file = '../data/woa_ini.nc'

    # Grid parameters: check grid_file and *.in file to make sure these are correct
    Tcline = 40
    theta_s = 0.9
    theta_b = 4
    hc = 40
    N = 31

    # Northernmost index of WOA grid to read (1-based)
    nbdry_woa = 32

    run(grid_file, woa_file, output_file, Tcline, theta_s, theta_b, hc, N, nbdry_woa)

    

