from netCDF4 import Dataset
from numpy import *
from scipy.interpolate import RegularGridInterpolator
from calc_z import *

# Build a ROMS lateral boundary condition file from ECCO2 temperature, salinity,
# and velocity output, containing boundary conditions at the northern boundary.
# Output 12 sets of monthly-averaged data, one for each month of the given year.
# NB: Users will likely need to edit paths to ECCO2 data! Scroll down below
# the interp_ecco2roms function to find where these filenames are defined.

# NB for raijin users: RegularGridInterpolator needs python/2.7.6 but the
# default is 2.7.3. Before running this script, switch them as follows:
# module unload python/2.7.3
# module unload python/2.7.3-matplotlib
# module load python/2.7.6
# module load python/2.7.6-matplotlib

# This script can be easily called by a batch script; see convert_ecco.job
# for an example.

def convert_file (year):

    # Make sure input argument is an integer (sometimes the batch script likes
    # to pass it as a string)
    year = int(year)

    # Paths of ROMS grid file, input ECCO2 files (without the tail yyyymm.nc),
    # and output ROMS-CICE boundary condition file; other users will need to
    # change these
    grid_file = '../ROMS-CICE-MCT/apps/common/grid/circ30S_quarterdegree.nc'
    theta_base = '../ROMS-CICE-MCT/data/ECCO2/THETA.1440x720x50.' + str(year)
    salt_base = '../ROMS-CICE-MCT/data/ECCO2/SALT.1440x720x50.' + str(year)
    uvel_base = '../ROMS-CICE-MCT/data/ECCO2/UVEL.1440x720x50.' + str(year)
    vvel_base = '../ROMS-CICE-MCT/data/ECCO2/VVEL.1440x720x50.' + str(year)
    output_file = '../ROMS-CICE-MCT/data/ecco2_cube92_lbc_' + str(year) + '.nc'

    # Grid parameters; check grid_file and *.in to make sure these are correct
    Tcline = 40
    theta_s = 0.9
    theta_b = 4
    hc = 40
    N = 31
    # Radius of the Earth (needed for spherical coordinates transformations)
    r = 6.371e6
    # Constant to convert from degrees to radians
    deg2rad = pi/180.0
    # Northernmost index of ECCO2 grid to read (1-based)
    nbdry_ecco = 241

    # Read ECCO2 grid
    print 'Reading ECCO2 grid'
    ecco_fid = Dataset(theta_base + '01.nc', 'r')
    lon_ecco_raw = ecco_fid.variables['LONGITUDE_T'][:]
    lat_ecco = ecco_fid.variables['LATITUDE_T'][0:nbdry_ecco]
    depth_ecco_raw = ecco_fid.variables['DEPTH_T'][:]
    ecco_fid.close()

    # The ECCO2 longitude axis doesn't wrap around; there is a gap between
    # almost-180W and almost-180E, and the ROMS grid has points in this gap.
    # So copy the last longitude value (mod 360) to the beginning, and the
    # first longitude value (mod 360) to the end.
    lon_ecco = zeros(size(lon_ecco_raw)+2)
    lon_ecco[0] = lon_ecco_raw[-1]-360
    lon_ecco[1:-1] = lon_ecco_raw
    lon_ecco[-1] = lon_ecco_raw[0]+360

    # The shallowest ECCO2 depth value is 5 m, but ROMS needs 0 m. So add the
    # index depth = 0 m to the beginning. Later we will just copy the 5 m
    # values for theta and salt into this index. Similarly, the deepest ECCO2
    # depth value is not deep enough for ROMS, so make a 6000 m index at the end.
    depth_ecco = zeros(size(depth_ecco_raw)+2)
    depth_ecco[0] = 0.0
    depth_ecco[1:-1] = depth_ecco_raw
    depth_ecco[-1] = 6000.0

    # Read ROMS grid
    print 'Reading ROMS grid'
    grid_fid = Dataset(grid_file, 'r')
    lon_rho = grid_fid.variables['lon_rho'][:,:]
    lat_rho = grid_fid.variables['lat_rho'][:,:]
    lon_u = grid_fid.variables['lon_u'][:,:]
    lat_u = grid_fid.variables['lat_u'][:,:]
    lon_v = grid_fid.variables['lon_v'][:,:]
    lat_v = grid_fid.variables['lat_v'][:,:]
    h = grid_fid.variables['h'][:,:]    
    zice = grid_fid.variables['zice'][:,:]
    mask_rho = grid_fid.variables['mask_rho'][:,:]
    mask_zice = grid_fid.variables['mask_zice'][:,:]
    grid_fid.close()

    # Save the lengths of the longitude axis for each grid
    num_lon_rho = size(lon_rho, 1)
    num_lon_u = size(lon_u, 1)
    num_lon_v = size(lon_v, 1)
    # Mask h and zice with zeros
    h = h*mask_rho
    zice = zice*mask_zice
    # Interpolate h and zice to u and v grids
    h_u = 0.5*(h[:,0:-1] + h[:,1:])
    h_v = 0.5*(h[0:-1,:] + h[1:,:])
    zice_u = 0.5*(zice[:,0:-1] + zice[:,1:])
    zice_v = 0.5*(zice[0:-1,:] + zice[1:,:])

    # Find thickness of northern boundary cells in latitude-direction
    dlat = lat_v[-1,:] - lat_v[-2,:]
    # dy = r*dlat where dlat is converted to radians
    # Also copy this into an array of dimension depth x longitude
    dy = tile(r*dlat*deg2rad, (N,1))
    # dlon is constant at this latitude
    dlon = 0.25
    lat = lat_v[-1,:]
    # dx = r*cos(lat)*dlon where lat and dlon are converted to radians
    # Also copy this into an array of dimension depth x longitude
    dx = tile(r*cos(lat*deg2rad)*dlon*deg2rad, (N,1))

    # Get a 3D array of ROMS z-coordinates on the rho-grid, as well as 1D
    # arrays of s-coordinates and stretching curves
    z_rho, sc_r, Cs_r = calc_z(h, zice, lon_rho, lat_rho, theta_s, theta_b, hc, N)
    # Only keep the depth x longitude slice of z_rho at the northern boundary
    z_rho = z_rho[:,-1,:]
    # Repeat for the u and v grids
    z_u, sc_r, Cs_r = calc_z(h_u, zice_u, lon_u, lat_u, theta_s, theta_b, hc, N)
    z_u = z_u[:,-1,:]

    # Now calculate dz on the u and v grids
    # We have z at the midpoint of each cell, now find it on the top and
    # bottom edges
    z_edges_u = zeros((size(z_u, 0)+1, size(z_u, 1)))
    z_edges_u[1:-1,:] = 0.5*(z_u[0:-1,:] + z_u[1:,:])
    # At surface, z = 0; at bottom, set z to be the same as the midpoint of
    # the deepest cell
    z_edges_u[0,:] = z_u[0,:]
    # Now find dz
    dz_u = z_edges_u[1:,:] - z_edges_u[0:-1,:]
    # Repeat for the v grid
    z_v, sc_r, Cs_r = calc_z(h_v, zice_v, lon_v, lat_v, theta_s, theta_b, hc, N)
    z_v = z_v[:,-1,:]
    z_edges_v = zeros((size(z_v, 0)+1, size(z_v, 1)))
    z_edges_v[1:-1,:] = 0.5*(z_v[0:-1,:] + z_v[1:,:])
    z_edges_v[0,:] = z_v[0,:]
    dz_v = z_edges_v[1:,:] - z_edges_v[0:-1,:]

    # Copy longitude and latitude into arrays of dimension depth x longitude
    lon_rho = tile(lon_rho[-1,:], (N,1))
    lat_rho = tile(lat_rho[-1,:], (N,1))
    lon_u = tile(lon_u[-1,:], (N,1))
    lat_u = tile(lat_u[-1,:], (N,1))
    lon_v = tile(lon_v[-1,:], (N,1))
    lat_v = tile(lat_v[-1,:], (N,1))

    # Set up output file
    print 'Setting up ', output_file
    out_fid = Dataset(output_file, 'w')
    out_fid.createDimension('xi_u', num_lon_u)
    out_fid.createDimension('xi_v', num_lon_v)
    out_fid.createDimension('xi_rho', num_lon_rho)
    out_fid.createDimension('s_rho', N)
    out_fid.createDimension('ocean_time', None)
    out_fid.createDimension('one', 1);
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
    out_fid.createVariable('sc_r', 'f8', ('s_rho'))
    out_fid.variables['sc_r'].long_name = 'S-coordinate at rho-points'
    out_fid.variables['sc_r'].units = 'nondimensional'
    out_fid.variables['sc_r'].valid_min = -1
    out_fid.variables['sc_r'].valid_max = 0
    out_fid.variables['sc_r'][:] = sc_r
    out_fid.createVariable('Cs_r', 'f8', ('s_rho'))
    out_fid.variables['Cs_r'].long_name = 'S-coordinate stretching curves at RHO-points'
    out_fid.variables['Cs_r'].units = 'nondimensional'
    out_fid.variables['Cs_r'].valid_min = -1
    out_fid.variables['Cs_r'].valid_max = 0
    out_fid.variables['Cs_r'][:] = Cs_r
    out_fid.createVariable('ocean_time', 'f8', ('ocean_time'))
    out_fid.variables['ocean_time'].long_name = 'time since initialization'
    out_fid.variables['ocean_time'].units = 'days'
    out_fid.createVariable('temp_north', 'f8', ('ocean_time', 's_rho', 'xi_rho'))
    out_fid.variables['temp_north'].long_name = 'northern boundary potential temperature'
    out_fid.variables['temp_north'].units = 'Celsius'
    out_fid.createVariable('salt_north', 'f8', ('ocean_time', 's_rho', 'xi_rho'))
    out_fid.variables['salt_north'].long_name = 'northern boundary salinity'
    out_fid.variables['salt_north'].units = 'PSU'
    out_fid.createVariable('u_north', 'f8', ('ocean_time', 's_rho', 'xi_u'))
    out_fid.variables['u_north'].long_name = 'northern boundary u-momentum component'
    out_fid.variables['u_north'].units = 'meter second-1'
    out_fid.createVariable('v_north', 'f8', ('ocean_time', 's_rho', 'xi_v'))
    out_fid.variables['v_north'].long_name = 'northern boundary v-momentum component'
    out_fid.variables['v_north'].units = 'meter second-1'
    out_fid.createVariable('ubar_north', 'f8', ('ocean_time', 'xi_u'))
    out_fid.variables['ubar_north'].long_name = 'northern boundary vertically integrated u-momentum component'
    out_fid.variables['ubar_north'].units = 'meter second-1'
    out_fid.createVariable('vbar_north', 'f8', ('ocean_time', 'xi_v'))
    out_fid.variables['vbar_north'].long_name = 'northern boundary vertically integrated v-momentum component'
    out_fid.variables['vbar_north'].units = 'meter second-1'
    out_fid.createVariable('zeta_north', 'f8', ('ocean_time', 'xi_rho'))
    out_fid.variables['zeta_north'].long_name = 'northern boundary sea surface height'
    out_fid.variables['zeta_north'].units = 'meter'
    out_fid.close()

    # Loop through each month of this year
    for month in range(12):

        print 'Processing month ', str(month+1), ' of 12'
        # Construct the rest of the file paths
        if month+1 < 10:
            tail = '0' + str(month+1) + '.nc'
        else:
            tail = str(month+1) + '.nc'

        # Read temperature, salinity, velocity data
        theta_fid = Dataset(theta_base + tail, 'r')
        theta_raw = transpose(theta_fid.variables['THETA'][0,:,0:nbdry_ecco,:])
        theta_fid.close()
        salt_fid = Dataset(salt_base + tail, 'r')
        salt_raw = transpose(salt_fid.variables['SALT'][0,:,0:nbdry_ecco,:])
        salt_fid.close()
        uvel_fid = Dataset(uvel_base + tail, 'r')
        uvel_raw = transpose(uvel_fid.variables['UVEL'][0,:,0:nbdry_ecco,:])
        uvel_fid.close()
        vvel_fid = Dataset(vvel_base + tail, 'r')
        vvel_raw = transpose(vvel_fid.variables['VVEL'][0,:,0:nbdry_ecco,:])
        vvel_fid.close()

        # Copy the data to the new longitude and depth indices, making sure
        # to preserve the mask.
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
        uvel = ma.array(zeros((size(lon_ecco), size(lat_ecco), size(depth_ecco))))
        uvel[1:-1,:,1:-1] = ma.copy(uvel_raw)
        uvel[0,:,1:-1] = ma.copy(uvel_raw[-1,:,:])
        uvel[-1,:,1:-1] = ma.copy(uvel_raw[0,:,:])
        uvel[:,:,0] = ma.copy(uvel[:,:,1])
        uvel[:,:,-1] = ma.copy(uvel[:,:,-2])
        vvel = ma.array(zeros((size(lon_ecco), size(lat_ecco), size(depth_ecco))))
        vvel[1:-1,:,1:-1] = ma.copy(vvel_raw)
        vvel[0,:,1:-1] = ma.copy(vvel_raw[-1,:,:])
        vvel[-1,:,1:-1] = ma.copy(vvel_raw[0,:,:])
        vvel[:,:,0] = ma.copy(vvel[:,:,1])
        vvel[:,:,-1] = ma.copy(vvel[:,:,-2])

        # Regridding happens here...
        print 'Interpolating temperature'
        temp_interp = interp_ecco2roms(theta, lon_ecco, lat_ecco, depth_ecco, lon_rho, lat_rho, z_rho, -0.5)
        print 'Interpolating salinity'
        salt_interp = interp_ecco2roms(salt, lon_ecco, lat_ecco, depth_ecco, lon_rho, lat_rho, z_rho, 34.5)
        print 'Interpolating u'
        u_interp = interp_ecco2roms(uvel, lon_ecco, lat_ecco, depth_ecco, lon_u, lat_u, z_u, 0)
        print 'Interpolating v'
        v_interp = interp_ecco2roms(vvel, lon_ecco, lat_ecco, depth_ecco, lon_v, lat_v, z_v, 0)

        # Find the volume-average of v
        v_int = sum(v_interp*dx*dy*dz_v)
        volume = sum(dx*dy*dz_v)
        v_avg = v_int/volume
        # Subtract this average from the data to make sure that volume is
        # conserved at the northern boundary
        v_interp = v_interp - v_avg

        # Calculate vertical averages of u and v to get ubar and vbar
        # Be sure to treat land mask carefully so we don't divide by 0
        ubar_interp = sum(u_interp*dz_u, axis=0)
        wct_u = h_u[-1,:]
        index = wct_u == 0
        ubar_interp[~index] = ubar_interp[~index]/wct_u[~index]
        ubar_interp[index] = 0.0
        vbar_interp = sum(v_interp*dz_v, axis=0)
        wct_v = h_v[-1,:]
        index = wct_v == 0
        vbar_interp[~index] = vbar_interp[~index]/wct_v[~index]
        vbar_interp[index] = 0.0

        # Calculate time values centered in the middle of each month,
        # relative to 1995
        time = 365.25*(year-1995) + 365.25/12*(month+0.5)

        # Save data to NetCDF file
        out_fid = Dataset(output_file, 'a')
        out_fid.variables['ocean_time'][month] = time
        out_fid.variables['temp_north'][month,:,:] = temp_interp
        out_fid.variables['salt_north'][month,:,:] = salt_interp
        out_fid.variables['u_north'][month,:,:] = u_interp
        out_fid.variables['v_north'][month,:,:] = v_interp
        out_fid.variables['ubar_north'][month,:] = ubar_interp
        out_fid.variables['vbar_north'][month,:] = vbar_interp
        out_fid.variables['zeta_north'][month,:] = 0.0
        out_fid.close()


# Given an array on the ECCO2 grid, fill the land mask with constant values,
# and then interpolate to the ROMS grid for the northern boundary.
# Input:
# A = array of size nxmxo containing values on the ECCO2 grid (dimension
#     longitude x latitude x depth)
# lon_ecco = array of length n containing ECCO2 longitude values (must have
#            first and last indices repeated mod 360 to cover the gap
#            180W=180E, as shown in the main function)
# lat_ecco = array of length m containing ECCO2 latitude values
# depth_ecco = array of length o containing ECCO2 depth values in positive
#              z-coordinates (must have index z=0 appended to the beginning,
#              as shown in the main script)
# lon_roms = array of size pxr containing ROMS longitude values
# lat_roms = array of size pxr containing ROMS latitude values
# z_roms = array of size pxr containing ROMS depth values in negative
#          z-coordinates (converted from grid file using calc_z.py, as shown
#          in the main script)
# fill_value = scalar containing the value with which to fill the ECCO2 land
#              mask: something close to the mean value of A so that the splines
#              don't freak out (suggest -0.5 for temperature, 34.5 for salinity,
#              0 for u and v)
# Output:
# B = array of size pxr containing values on the ROMS grid (dimension depth x
#     longitude)
def interp_ecco2roms (A, lon_ecco, lat_ecco, depth_ecco, lon_roms, lat_roms, z_roms, fill_value):

    # Fill the ECCO2 land mask with a constant value close to the mean of A.
    # This is less than ideal because it will skew the interpolation slightly
    # along the coast, but it doesn't matter because the northern boundary is
    # far away from the coast.
    Afill = A
    Afill[A.mask] = fill_value

    # Build a function for linear interpolation of Afill
    interp_function = RegularGridInterpolator((lon_ecco, lat_ecco, depth_ecco), Afill)
    # Call this function at the ROMS grid points (pass positive values for
    # ROMS depth)
    B = interp_function((lon_roms, lat_roms, -z_roms))

    return B


        
    
    
        
        
    
    
    
