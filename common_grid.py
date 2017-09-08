from numpy import *
from netCDF4 import Dataset, num2date
from scipy.interpolate import griddata
from rotate_vector_roms import *
from rotate_vector_cice import *
from monthly_avg_roms import *
from monthly_avg_cice import *

# Interpolate ROMS output to a regular quarter-degree grid for easy comparison
# with FESOM. Write monthly averages of surface temperature, salinity, surface
# heat and salt flux, sea ice concentration and thickness, ocean and sea ice
# velocity, the surface stress vector and its curl.
# Input:
# roms_file = path to ROMS output file containing 5-day averages for the entire
#             simulation, starting on 1 January
# cice_file = path to CICE output file containing 5-day averages for the same
#             timesteps as ROMS_file
# out_file = path to desired output file
def common_grid (roms_file, cice_file, out_file):

    # Resolution of common grid (degrees, same for lat and lon)
    res = 0.25
    # Northern boundary to interpolate to
    nbdry = -50
    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    N = 31

    print 'Calculating grids'

    # Make the latitude and longitude arrays for the common grid
    lon_common = arange(-180, 180+res, res)
    lat_common = arange(-90, nbdry+res, res)
    # Get a 2D version of each to calculate dx and dy in metres
    lon_2d, lat_2d = meshgrid(lon_common, lat_common)
    # dx = r*cos(lat)*dlon where lat and dlon (i.e. res) are in radians    
    dx = r*cos(lat_2d*deg2rad)*res*deg2rad
    # dy = r*dlat where dlat (i.e. res) is in radians
    # This is constant so reshape to an array of the right dimensions
    dy = zeros(shape(dx)) + r*res*deg2rad

    # Read the ROMS grid
    id = Dataset(roms_file, 'r')
    # We only need lat and lon on the rho grid
    lon_rho = id.variables['lon_rho'][:,:]
    lat_rho = id.variables['lat_rho'][:,:]
    # Get shape of u and v grids
    u_shape = id.variables['lon_u'].shape
    v_shape = id.variables['lon_v'].shape
    # Read land mask 
    mask_roms = id.variables['mask_rho'][:,:]
    # Mask out ice shelves too
    zice = id.variables['zice'][:,:]
    mask_roms[zice != 0] = 0.0
    # Read angle (for rotation of vector components)
    angle_roms = id.variables['angle'][:,:]
    # Get time as an array of Date objects
    time_id = id.variables['ocean_time']
    time = num2date(time_id[:], units=time_id.units, calendar=time_id.calendar.lower())
    id.close()

    # Read the CICE grid
    id = Dataset(cice_file, 'r')
    # We only need lat and lon on the velocity grid
    lon_cice = id.variables['ULON'][:,:]
    lat_cice = id.variables['ULAT'][:,:]
    # Read angle (for rotation of vector components)
    angle_cice = id.variables['ANGLE'][:,:]
    id.close()

    # Make sure longitude is between -180 and 180
    index = lon_rho > 180
    lon_rho[index] = lon_rho[index] - 360
    index = lon_cice > 180
    lon_cice[index] = lon_cice[index] - 360

    print 'Counting months'
    # Assume we start at the beginning of a year
    # Figure out how many complete years have happened since then
    num_full_years = time[-1].year - time[0].year
    if time[-1].month == 12 and time[-1].day in range(29,31+1):
        # We happen to end at the very end of a year
        num_full_years += 1
    else:
        # Count the complete months that have happened this year
        num_extra_months = time[-1].month - 1
        # Don't bother with the hassle of considering cases where we end at
        # the very end of a month. Just ignore the month.
    num_months = 12*num_full_years + num_extra_months

    print 'Interpolating land mask to new grid'
    mask_common = interp_roms2common(lon_common, lat_common, lon_rho, lat_rho, mask_roms)
    mask_common[isnan(mask_common)] = 0
    # Cut it off at 1
    mask_common[mask_common < 0.5] = 0
    mask_common[mask_common >= 0.5] = 1

#    print 'Setting up ' + out_file
#    id = Dataset(out_file, 'w')
#    id.createDimension('longitude', size(lon_common))
#    id.createDimension('latitude', size(lat_common))
#    id.createDimension('time', None)
#    id.createVariable('longitude', 'f8', ('longitude'))
#    id.variables['longitude'].units = 'degrees'
#    id.variables['longitude'][:] = lon_common
#    id.createVariable('latitude', 'f8', ('latitude'))
#    id.variables['latitude'].units = 'degrees'
#    id.variables['latitude'][:] = lat_common
#    id.createVariable('time', 'f8', ('time'))
#    id.variables['time'].units = 'months'
#    id.createVariable('mask', 'f8', ('latitude', 'longitude'))
#    id.variables['mask'].units = '1'
#    id.variables['mask'][:,:] = mask_common
#    id.createVariable('sst', 'f8', ('time', 'latitude', 'longitude'))
#    id.variables['sst'].long_name = 'sea surface temperature'
#    id.variables['sst'].units = 'C'
#    id.createVariable('sss', 'f8', ('time', 'latitude', 'longitude'))
#    id.variables['sss'].long_name = 'sea surface salinity'
#    id.variables['sss'].units = 'psu'
#    id.createVariable('shflux', 'f8', ('time', 'latitude', 'longitude'))
#    id.variables['shflux'].long_name = 'surface heat flux into ocean'
#    id.variables['shflux'].units = 'W/m^2'
#    id.createVariable('ssflux', 'f8', ('time', 'latitude', 'longitude'))
#    id.variables['ssflux'].long_name = 'surface virtual salinity flux into ocean'
#    id.variables['ssflux'].units = 'psu m/s'
#    id.createVariable('aice', 'f8', ('time', 'latitude', 'longitude'))
#    id.variables['aice'].long_name = 'sea ice concentration'
#    id.variables['aice'].units = '1'
#    id.createVariable('hice', 'f8', ('time', 'latitude', 'longitude'))
#    id.variables['hice'].long_name = 'sea ice thickness'
#    id.variables['hice'].units = 'm'
#    id.createVariable('uocn', 'f8', ('time', 'latitude', 'longitude'))
#    id.variables['uocn'].long_name = 'ocean surface velocity eastward'
#    id.variables['uocn'].units = 'm/s'
#    id.createVariable('vocn', 'f8', ('time', 'latitude', 'longitude'))
#    id.variables['vocn'].long_name = 'ocean surface velocity northward'
#    id.variables['vocn'].units = 'm/s'
#    id.createVariable('uice', 'f8', ('time', 'latitude', 'longitude'))
#    id.variables['uice'].long_name = 'sea ice velocity eastward'
#    id.variables['uice'].units = 'm/s'
#    id.createVariable('vice', 'f8', ('time', 'latitude', 'longitude'))
#    id.variables['vice'].long_name = 'sea ice velocity northward'
#    id.variables['vice'].units = 'm/s'
#    id.createVariable('sustr', 'f8', ('time', 'latitude', 'longitude'))
#    id.variables['sustr'].long_name = 'zonal surface stress'
#    id.variables['sustr'].units = 'N/m^2'
#    id.createVariable('svstr', 'f8', ('time', 'latitude', 'longitude'))
#    id.variables['svstr'].long_name = 'meridional surface stress'
#    id.variables['svstr'].units = 'N/m^2'
#    id.createVariable('curl_str', 'f8', ('time', 'latitude', 'longitude'))
#    id.variables['curl_str'].long_name = 'curl of surface stress'
#    id.variables['curl_str'].units = 'N/m^3'
#    id.close()

    # Loop over months
    for month in range(18,num_months):
        print 'Processing month ' + str(month+1) + ' of ' + str(num_months)
        id = Dataset(out_file, 'a')
        # Write time value for this month
        id.variables['time'][month] = month+1

        print '...sea surface temperature'
        # Get monthly average of 3D variable
        temp_roms = monthly_avg_roms(roms_file, 'temp', [N, size(lon_rho,0), size(lon_rho,1)], month%12, instance=month/12+1)
        # Select surface layer
        sst_roms = temp_roms[-1,:,:]
        # Interpolate to common grid
        sst_common = interp_roms2common(lon_common, lat_common, lon_rho, lat_rho, sst_roms)
        # Apply land mask
        sst = ma.masked_where(mask_common==0, sst_common)
        # Write to file
        id.variables['sst'][month,:,:] = sst

        print '...sea surface salinity'
        # Get monthly average of 3D variable
        salt_roms = monthly_avg_roms(roms_file, 'salt', [N, size(lon_rho,0), size(lon_rho,1)], month%12, instance=month/12+1)
        # Select surface layer
        sss_roms = salt_roms[-1,:,:]
        # Interpolate to common grid
        sss_common = interp_roms2common(lon_common, lat_common, lon_rho, lat_rho, sss_roms)
        # Apply land mask
        sss = ma.masked_where(mask_common==0, sss_common)
        # Write to file
        id.variables['sss'][month,:,:] = sss

        print '...surface heat flux'
        # Get monthly average
        shflux_roms = monthly_avg_roms(roms_file, 'shflux', shape(lon_rho), month%12, instance=month/12+1)
        # Interpolate to common grid
        shflux_common = interp_roms2common(lon_common, lat_common, lon_rho, lat_rho, shflux_roms)
        # Apply land mask
        shflux = ma.masked_where(mask_common==0, shflux_common)
        # Write to file
        id.variables['shflux'][month,:,:] = shflux

        print '...surface salt flux'
        # Get monthly average
        ssflux_roms = monthly_avg_roms(roms_file, 'ssflux', shape(lon_rho), month%12, instance=month/12+1)
        # Interpolate to common grid
        ssflux_common = interp_roms2common(lon_common, lat_common, lon_rho, lat_rho, ssflux_roms)
        # Apply land mask
        ssflux = ma.masked_where(mask_common==0, ssflux_common)
        # Write to file
        id.variables['ssflux'][month,:,:] = ssflux

        print '...sea ice concentration'
        # Get monthly average (use CICE file)
        aice_cice = monthly_avg_cice(cice_file, 'aice', shape(lon_cice), month%12, instance=month/12+1)
        # Interpolate to common grid (note CICE grid not ROMS)
        aice_common = interp_roms2common(lon_common, lat_common, lon_cice, lat_cice, aice_cice)
        # Apply land mask
        aice = ma.masked_where(mask_common==0, aice_common)
        # Write to file
        id.variables['aice'][month,:,:] = aice

        print '...sea ice thickness'
        # Get monthly average (use CICE file)
        hice_cice = monthly_avg_cice(cice_file, 'hi', shape(lon_cice), month%12, instance=month/12+1)
        # Interpolate to common grid (note CICE grid not ROMS)
        hice_common = interp_roms2common(lon_common, lat_common, lon_cice, lat_cice, hice_cice)
        # Apply land mask
        hice = ma.masked_where(mask_common==0, hice_common)
        # Write to file
        id.variables['hice'][month,:,:] = hice

        print '...surface ocean velocity vector'
        # Surface ocean velocity
        # Get monthly averages of both 3D vector components
        uocn_3d_tmp = monthly_avg_roms(roms_file, 'u', [N, u_shape[0], u_shape[1]], month%12, instance=month/12+1)
        vocn_3d_tmp = monthly_avg_roms(roms_file, 'v', [N, v_shape[0], v_shape[1]], month%12, instance=month/12+1)
        # Select surface layer
        uocn_tmp = uocn_3d_tmp[-1,:,:]
        vocn_tmp = vocn_3d_tmp[-1,:,:]
        # Rotate to lon-lat space (note they are on the rho grid now)
        uocn_roms, vocn_roms = rotate_vector_roms(uocn_tmp, vocn_tmp, angle_roms)
        # Interpolate to common grid
        uocn_common = interp_roms2common(lon_common, lat_common, lon_rho, lat_rho, uocn_roms)
        vocn_common = interp_roms2common(lon_common, lat_common, lon_rho, lat_rho, vocn_roms)
        # Apply land mask
        uocn = ma.masked_where(mask_common==0, uocn_common)
        vocn = ma.masked_where(mask_common==0, vocn_common)
        # Write to file
        id.variables['uocn'][month,:,:] = uocn
        id.variables['vocn'][month,:,:] = vocn

        print '...sea ice velocity vector'
        # Sea ice velocity (CICE variable not ROMS)
        # Get monthly averages of both vector components
        uice_tmp = monthly_avg_cice(cice_file, 'uvel', shape(lon_cice), month%12, instance=month/12+1)
        vice_tmp = monthly_avg_cice(cice_file, 'vvel', shape(lon_cice), month%12, instance=month/12+1)
        # Rotate to lon-lat space
        uice_cice, vice_cice = rotate_vector_cice(uice_tmp, vice_tmp, angle_cice)
        # Interpolate to common grid (note CICE grid not ROMS)
        uice_common = interp_roms2common(lon_common, lat_common, lon_cice, lat_cice, uice_cice)
        vice_common = interp_roms2common(lon_common, lat_common, lon_cice, lat_cice, vice_cice)
        # Apply land mask
        uice = ma.masked_where(mask_common==0, uice_common)
        vice = ma.masked_where(mask_common==0, vice_common)
        # Write to file
        id.variables['uice'][month,:,:] = uice
        id.variables['vice'][month,:,:] = vice

        print '...surface stress vector'
        # Surface stresses
        # Get monthly averages of both vector components
        sustr_tmp = monthly_avg_roms(roms_file, 'sustr', u_shape, month%12, instance=month/12+1)
        svstr_tmp = monthly_avg_roms(roms_file, 'svstr', v_shape, month%12, instance=month/12+1)
        # Rotate to lon-lat space (note they are on the rho grid now)
        sustr_roms, svstr_roms = rotate_vector_roms(sustr_tmp, svstr_tmp, angle_roms)
        # Interpolate to common grid
        sustr_common = interp_roms2common(lon_common, lat_common, lon_rho, lat_rho, sustr_roms)
        svstr_common = interp_roms2common(lon_common, lat_common, lon_rho, lat_rho, svstr_roms)
        # Apply land mask
        sustr = ma.masked_where(mask_common==0, sustr_common)
        svstr = ma.masked_where(mask_common==0, svstr_common)
        # Write to file
        id.variables['sustr'][month,:,:] = sustr
        id.variables['svstr'][month,:,:] = svstr

        print '...curl of surface stress vector'
        # Curl of surface stress = d/dx (svstr) - d/dy (sustr)
        # First calculate the two derivatives
        dsvstr_dx = ma.empty(shape(svstr_common))
        # Forward difference approximation
        dsvstr_dx[:,:-1] = (svstr_common[:,1:] - svstr_common[:,:-1])/dx[:,:-1]
        # Backward difference for the last row
        dsvstr_dx[:,-1] = (svstr_common[:,-1] - svstr_common[:,-2])/dx[:,-1]
        dsustr_dy = ma.empty(shape(sustr_common))
        dsustr_dy[:-1,:] = (sustr_common[1:,:] - sustr_common[:-1,:])/dy[:-1,:]
        dsustr_dy[-1,:] = (sustr_common[-1,:] - sustr_common[-2,:])/dy[-1,:]
        curl_str = dsvstr_dx - dsustr_dy
        curl_str = ma.masked_where(mask_common==0, curl_str)
        # Write to file
        id.variables['curl_str'][month,:,:] = curl_str

        id.close()

    print 'Finished'
    

# Interpolate from the ROMS grid to the common grid.
# This works for the CICE grid too.
# Input:
# lon_1d = 1D array of longitude on the common grid, -180 to 180 (size n)
# lat_1d = 1D array of latitude on the common grid (size m)
# lon_roms = 2D array of longitude on the ROMS grid, -180 to 180 (size pxq)
# lat_roms = 2D array of latitude on the ROMS grid (size pxq)
# data_roms = 2D array of any variable on the ROMS grid (size pxq)
# Output:
# data_common = 2D array of data_roms interpolated to the common grid (size mxn)
def interp_roms2common (lon_1d, lat_1d, lon_roms, lat_roms, data_roms):

    # Get a 2D field of common latitude and longitude
    lon_2d, lat_2d = meshgrid(lon_1d, lat_1d)

    # Make an array of all the ROMS coordinates, flattened
    points = empty([size(lon_roms), 2])
    points[:,0] = ravel(lon_roms)
    points[:,1] = ravel(lat_roms)
    # Also flatten the ROMS data
    values = ravel(data_roms)
    # Now make an array of all the common grid coordinates, flattened
    xi = empty([size(lon_2d), 2])
    xi[:,0] = ravel(lon_2d)
    xi[:,1] = ravel(lat_2d)
    # Now call griddata
    result = griddata(points, values, xi)
    # Un-flatten the result
    data_common = reshape(result, shape(lon_2d))

    return data_common


# Command-line interface
if __name__ == "__main__":

    roms_file = raw_input("Path to ROMS averages file for the entire simulation, starting 1 Jan: ")
    cice_file = raw_input("Path to CICE history file containing the same time indices: ")
    out_file = raw_input("Path to desired output file: ")
    common_grid(roms_file, cice_file, out_file)
    
