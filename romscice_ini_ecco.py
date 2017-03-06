# Build a 1992 initialisation file for ROMS from the ECCO2 temperature and 
# salinity reanalysis of January 1992. Set initial velocities and sea surface
# height to zero. Under the ice shelves, extrapolate temperature and salinity
# from the ice shelf front.
# NB: Users will likely need to edit paths to ECCO2 data! Scroll down below
# the interp_ecco2roms function to find where these filenames are defined.

# NB for raijin users: RegularGridInterpolator needs python/2.7.6 but the
# default is 2.7.3. Before running this script, switch them as follows:
# module unload python/2.7.3
# module unload python/2.7.3-matplotlib
# module load python/2.7.6
# module load python/2.7.6-matplotlib

from netCDF4 import Dataset
from numpy import *
from scipy.interpolate import RegularGridInterpolator
from calc_z import *


# Main routine
def run (grid_file, theta_file, salt_file, output_file, Tcline, theta_s, theta_b, hc, N, nbdry_ecco):

    # Read ECCO2 data and grid
    print 'Reading ECCO2 data'
    theta_fid = Dataset(theta_file, 'r')
    lon_ecco_raw = theta_fid.variables['LONGITUDE_T'][:]
    lat_ecco = theta_fid.variables['LATITUDE_T'][0:nbdry_ecco]
    depth_ecco_raw = theta_fid.variables['DEPTH_T'][:]
    theta_raw = transpose(theta_fid.variables['THETA'][0,:,0:nbdry_ecco,:])
    theta_fid.close()
    salt_fid = Dataset(salt_file, 'r')
    salt_raw = transpose(salt_fid.variables['SALT'][0,:,0:nbdry_ecco,:])
    salt_fid.close()

    # The ECCO2 longitude axis doesn't wrap around; there is a gap between
    # almost-180W and almost-180E, and the ROMS grid has points in this gap.
    # So copy the last longitude value (mod 360) to the beginning, and the
    # first longitude value (mod 360) to the end.
    lon_ecco = zeros(size(lon_ecco_raw)+2)
    lon_ecco[0] = lon_ecco_raw[-1]-360
    lon_ecco[1:-1] = lon_ecco_raw
    lon_ecco[-1] = lon_ecco_raw[0]+360

    # The shallowest ECCO2 depth value is 5 m, but ROMS needs 0 m. So add the
    # index depth = 0 m to the beginning. Later we will just copy the 5 m values
    # for theta and salt into this index. Similarly, add the index depth = 6000 m
    # to the end.
    depth_ecco = zeros(size(depth_ecco_raw)+2)
    depth_ecco[0] = 0.0
    depth_ecco[1:-1] = depth_ecco_raw
    depth_ecco[-1] = 6000.0

    # Copy the theta and salt values to the new longitude and depth indices,
    # making sure to preserve the mask.
    theta = ma.array(zeros((size(lon_ecco), size(lat_ecco), size(depth_ecco))))
    theta[1:-1,:,1:-1] = ma.copy(theta_raw)
    theta[0,:,1:-1] = ma.copy(theta_raw[-1,:,:])
    theta[-1,:,1:-1] = ma.copy(theta_raw[0,:,:])
    theta[:,:,0] = ma.copy(theta[:,:,1])
    theta[:,:,-1] = ma.copy(theta[:,:,-2])
    salt = ma.array(zeros((size(lon_ecco), size(lat_ecco), size(depth_ecco))))
    salt[1:-1,:,1:-1] = ma.copy(salt_raw)
    salt[0,:,1:-1] = ma.copy(salt_raw[-1,:,:])
    salt[-1,:,1:-1] = ma.copy(salt_raw[0,:,:])
    salt[:,:,0] = ma.copy(salt[:,:,1])
    salt[:,:,-1] = ma.copy(salt[:,:,-2])

    # Read ROMS grid
    print 'Reading ROMS grid'
    grid_fid = Dataset(grid_file, 'r')
    lon_roms = grid_fid.variables['lon_rho'][:,:]
    lat_roms = grid_fid.variables['lat_rho'][:,:]
    h = grid_fid.variables['h'][:,:]
    zice = grid_fid.variables['zice'][:,:]
    mask_rho = grid_fid.variables['mask_rho'][:,:]
    mask_zice = grid_fid.variables['mask_zice'][:,:]
    grid_fid.close()
    num_lon = size(lon_roms, 1)
    num_lat = size(lon_roms, 0)
    # Mask h and zice with zeros
    h = h*mask_rho
    zice = zice*mask_zice

    # Get a 3D array of ROMS z-coordinates, as well as 1D arrays of s-coordinates
    # and stretching curves
    z_roms_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N)
    # Copy the latitude and longitude values into 3D arrays of the same shape
    lon_roms_3d = tile(lon_roms, (N,1,1))
    lat_roms_3d = tile(lat_roms, (N,1,1))

    # Regridding happens here...
    print 'Interpolating temperature'
    temp = interp_ecco2roms(theta, lon_ecco, lat_ecco, depth_ecco, lon_roms_3d, lat_roms_3d, z_roms_3d, mask_rho, mask_zice, -0.5)
    print 'Interpolating salinity'
    salt = interp_ecco2roms(salt, lon_ecco, lat_ecco, depth_ecco, lon_roms_3d, lat_roms_3d, z_roms_3d, mask_rho, mask_zice, 34.5)

    # Set initial velocities and sea surface height to zero
    u = zeros((N, num_lat, num_lon-1))
    v = zeros((N, num_lat-1, num_lon))
    ubar = zeros((num_lat, num_lon-1))
    vbar = zeros((num_lat-1, num_lon))
    zeta = zeros((num_lat, num_lon))

    print 'Writing to NetCDF file'
    out_fid = Dataset(output_file, 'w')
    # Define dimensions
    out_fid.createDimension('xi_u', num_lon-1)
    out_fid.createDimension('xi_v', num_lon)
    out_fid.createDimension('xi_rho', num_lon)
    out_fid.createDimension('eta_u', num_lat)
    out_fid.createDimension('eta_v', num_lat-1)
    out_fid.createDimension('eta_rho', num_lat)
    out_fid.createDimension('s_rho', N)
    out_fid.createDimension('ocean_time', None)
    out_fid.createDimension('one', 1);
    # Define variables and assign values
    out_fid.createVariable('tstart', 'f8', ('one'))
    out_fid.variables['tstart'].long_name = 'start processing day'
    out_fid.variables['tstart'].units = 'day'
    out_fid.variables['tstart'][:] = 0.0
    out_fid.createVariable('tend', 'f8', ('one'))
    out_fid.variables['tend'].long_name = 'end processing day'
    out_fid.variables['tend'].units = 'day'
    out_fid.variables['tend'][:] = 0.0
    out_fid.createVariable('theta_s', 'f8', ('one'))
    out_fid.variables['theta_s'].long_name = 'S-coordinate surface control parameter'
    out_fid.variables['theta_s'][:] = theta_s
    out_fid.createVariable('theta_b', 'f8', ('one'))
    out_fid.variables['theta_b'].long_name = 'S-coordinate bottom control parameter'
    out_fid.variables['theta_b'].units = 'nondimensional'
    out_fid.variables['theta_b'][:] = theta_b
    out_fid.createVariable('Tcline', 'f8', ('one'))
    out_fid.variables['Tcline'].long_name = 'S-coordinate surface/bottom layer width'
    out_fid.variables['Tcline'].units = 'meter'
    out_fid.variables['Tcline'][:] = Tcline
    out_fid.createVariable('hc', 'f8', ('one'))
    out_fid.variables['hc'].long_name = 'S-coordinate parameter, critical depth'
    out_fid.variables['hc'].units = 'meter'
    out_fid.variables['hc'][:] = hc
    out_fid.createVariable('Cs_r', 'f8', ('s_rho'))
    out_fid.variables['Cs_r'].long_name = 'S-coordinate stretching curves at RHO-points'
    out_fid.variables['Cs_r'].units = 'nondimensional'
    out_fid.variables['Cs_r'].valid_min = -1.0
    out_fid.variables['Cs_r'].valid_max = 0.0
    out_fid.variables['Cs_r'][:] = Cs_r
    out_fid.createVariable('ocean_time', 'f8', ('ocean_time'))
    out_fid.variables['ocean_time'].long_name = 'time since initialization'
    out_fid.variables['ocean_time'].units = 'seconds'
    out_fid.variables['ocean_time'][0] = 0.0
    out_fid.createVariable('u', 'f8', ('ocean_time', 's_rho', 'eta_u', 'xi_u'))
    out_fid.variables['u'].long_name = 'u-momentum component'
    out_fid.variables['u'].units = 'meter second-1'
    out_fid.variables['u'][0,:,:,:] = u
    out_fid.createVariable('v', 'f8', ('ocean_time', 's_rho', 'eta_v', 'xi_v'))
    out_fid.variables['v'].long_name = 'v-momentum component'
    out_fid.variables['v'].units = 'meter second-1'
    out_fid.variables['v'][0,:,:,:] = v
    out_fid.createVariable('ubar', 'f8', ('ocean_time', 'eta_u', 'xi_u'))
    out_fid.variables['ubar'].long_name = 'vertically integrated u-momentum component'
    out_fid.variables['ubar'].units = 'meter second-1'
    out_fid.variables['ubar'][0,:,:] = ubar
    out_fid.createVariable('vbar', 'f8', ('ocean_time', 'eta_v', 'xi_v'))
    out_fid.variables['vbar'].long_name = 'vertically integrated v-momentum component'
    out_fid.variables['vbar'].units = 'meter second-1'
    out_fid.variables['vbar'][0,:,:] = vbar
    out_fid.createVariable('zeta', 'f8', ('ocean_time', 'eta_rho', 'xi_rho'))
    out_fid.variables['zeta'].long_name = 'free-surface'
    out_fid.variables['zeta'].units = 'meter'
    out_fid.variables['zeta'][0,:,:] = zeta
    out_fid.createVariable('temp', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'))
    out_fid.variables['temp'].long_name = 'potential temperature'
    out_fid.variables['temp'].units = 'Celsius'
    out_fid.variables['temp'][0,:,:,:] = temp
    out_fid.createVariable('salt', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'))
    out_fid.variables['salt'].long_name = 'salinity'
    out_fid.variables['salt'].units = 'PSU'
    out_fid.variables['salt'][0,:,:,:] = salt
    out_fid.createVariable('sc_r', 'f8', ('s_rho'))
    out_fid.variables['sc_r'].long_name = 'S-coordinate at rho-points'
    out_fid.variables['sc_r'].units = 'nondimensional'
    out_fid.variables['sc_r'].valid_min = -1.0
    out_fid.variables['sc_r'].valid_max = 0.0
    out_fid.variables['sc_r'][:] = sc_r
    out_fid.close()

    print 'Finished'


# Given an array on the ECCO2 grid, interoplate onto the ROMS grid, extrapolate
# under ice shelf cavities, and fill the land mask with constant values.
# Input:
# A = array of size nxmxo containing values on the ECCO2 grid (dimension
#     longitude x latitude x depth)
# lon_ecco = array of length n containing ECCO2 longitude values (must have
#            first and last indices repeated  mod 360 to cover the gap 
#            180W=180E, as shown below in the main script)
# lat_ecco = array of length m containing ECCO2 latitude values
# depth_ecco = array of length o containing ECCO2 depth values in positive 
#              z-coordinates (must have index z=0 appended to the beginning, 
#              as shown below in the main script)
# lon_roms_3d = array of size pxqxr containing ROMS longitude values
# lat_roms_3d = array of size pxqxr containing ROMS latitude values
# z_roms_3d = array of size pxqxr containing ROMS depth values in negative
#             z-coordinates (converted from grid file using calc_z.py, as
#             shown below in the main script)
# mask_rho = array of size qxr containing ROMS land mask (ocean/ice shelf 1,
#            land 0)
# mask_zice = array of size qxr containing ROMS ice shelf mask (ice shelf 1,
#             ocean/land 0)
# fill = scalar containing the value with which to fill the ROMS land mask
#        (doesn't really matter, don't use NaN)
# Output:
# B = array of size pxqxr containing values on the ROMS grid (dimension depth x
#     latitude x longitude)
def interp_ecco2roms (A, lon_ecco, lat_ecco, depth_ecco, lon_roms_3d, lat_roms_3d, z_roms_3d, mask_rho, mask_zice, fill):

    # Radius of the Earth in m
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0

    # Calculate N based on size of ROMS grid
    N = size(lon_roms_3d, 0)

    # Unmask A and fill with NaNs
    A_unmask = A.data
    A_unmask[A.mask]=NaN

    # Build a function for linear interpolation of A; set flag to fill
    # out-of-bounds values with NaN
    interp_function = RegularGridInterpolator((lon_ecco, lat_ecco, depth_ecco), A_unmask, bounds_error=False, fill_value=NaN)
    B = zeros(shape(lon_roms_3d))

    # Interpolate each z-level individually - 3D vectorisation uses too
    # much memory!
    for k in range(N):
        print '...vertical level ', str(k+1), ' of ', str(N)
        # Pass positive values for ROMS depth
        tmp = interp_function((lon_roms_3d[k,:,:], lat_roms_3d[k,:,:], -z_roms_3d[k,:,:]))
        # Fill NaNs with constant value
        index = isnan(tmp)
        tmp[index] = fill
        # Fill land mask with constant value
        tmp[mask_rho==0] = fill
        # Save this depth level
        B[k,:,:] = tmp

    print '...selecting the ice shelf front'
    # Find the ice shelf front
    front = zeros(shape(lon_roms_3d))
    # Find northernmost latitude with ice shelves
    nbdry_zice = where(sum(mask_zice, axis=1) > 0)[-1][-1] + 2
    # Loop over all cells; can't find a cleaner way to do this
    for j in range(nbdry_zice):
        for i in range(size(mask_zice,1)):
            # First make sure this is an ocean point
            if mask_rho[j,i] == 1:
                # Find bounds on window of radius 1
                j_min = max(j-1,0)
                j_max = min(j+2, size(mask_zice,0))
                i_min = max(i-1,0)
                i_max = min(i+2, size(mask_zice,1))
                # Points at the ice shelf front are not ice shelf points, but
                # have at least one neighbour that is
                if mask_zice[j,i] == 0 and any(mask_zice[j_min:j_max,i_min:i_max] == 1):
                    front[:,j,i] = 1

    print '...extrapolating under ice shelves'
    mask_zice_3d = tile(mask_zice, (N,1,1))
    for j in range(nbdry_zice, -1, -1):
        print '...latitude index ' + str(nbdry_zice-j+1) + ' of ' + str(nbdry_zice+1)
        # Find coordinates of ice shelf front points
        lon_front = lon_roms_3d[front==1]
        lat_front = lat_roms_3d[front==1]
        z_front = z_roms_3d[front==1]
        # Also find the values of the given variable here
        B_front = B[front==1]
        for i in range(size(mask_zice,1)):
            if mask_zice[j,i] == 1:
                for k in range(N):
                    # Calculate Cartesian distance from this point to each point
                    # at the ice shelf front
                    dlon = lon_front - lon_roms_3d[k,j,i]
                    index = dlon > 300
                    dlon[index] -= 360
                    index = dlon < -300
                    dlon[index] += 360
                    dlon = abs(dlon)
                    dlat = abs(lat_front - lat_roms_3d[k,j,i])
                    dz = abs(z_front - z_roms_3d[k,j,i])
                    dx = r*cos(lat_roms_3d[k,j,i]*deg2rad)*dlon*deg2rad
                    dy = r*dlat*deg2rad
                    dist = sqrt(dx**2 + dy**2 + dz**2)
                    # Find the index with the minimum distance: this is the
                    # nearest-neighbour
                    nearest = argmin(dist)
                    B[k,j,i] = B_front[nearest]
        # Update the ice shelf front points to include the new points we've just
        # extrapolated to
        front_tmp = front[:,j,:]
        mask_zice_tmp = mask_zice_3d[:,j,:]
        front_tmp[mask_zice_tmp==1] = 1
        front[:,j,:] = front_tmp
        
    return B


# User parameters
if __name__ == "__main__":

    # File paths to edit here:
    # Path to ROMS grid file
    grid_file = '/short/m68/kaa561/metroms_iceshelf/apps/common/grid/circ30S_quarterdegree.nc'
    # Path to ECCO2 files for temperature and salinity in January 1992
    theta_file = '/short/m68/kaa561/metroms_iceshelf/data/THETA.1440x720x50.199201.nc'
    salt_file = '/short/m68/kaa561/metroms_iceshelf/data/SALT.1440x720x50.199201.nc'
    # Path to desired output file
    output_file = '/short/m68/kaa561/metroms_iceshelf/data/ecco2_ini.nc'

    # Grid parameters; check grid_file and *.in file to make sure these are correct
    Tcline = 40
    theta_s = 4.0
    theta_b = 0.9
    hc = 40
    N = 31

    # Northernmost index of ECCO2 grid to read (1-based)
    nbdry_ecco = 241

    run(grid_file, theta_file, salt_file, output_file, Tcline, theta_s, theta_b, hc, N, nbdry_ecco)

    

        
    
    
        
        
    
    
    
