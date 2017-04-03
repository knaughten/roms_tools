from numpy import *
from netCDF4 import Dataset, num2date
from scipy.interpolate import griddata
from rotate_vector_roms import *
from rotate_vector_cice import *
from monthly_avg_roms import *
from monthly_avg_cice import *

def common_grid (roms_file, cice_file, out_file):

    # Resolution of common grid (degrees, same for lat and lon)
    res = 0.25
    # Northern boundary to interpolate to
    nbdry = -50
    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0

    print 'Calculating grids'

    # Make the latitude and longitude arrays for the common grid
    lon_common = arange(-180, 180+res, res)
    lat_common = arange(-90, nbdry+res, res)
    # Get a 2D version of each to calculate dx and dy in metres
    lat_2d, lon_2d = meshgrid(lat_common, lon_common)
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
    # Make the interpolation strict, i.e. cells with any land in their
    # interpolation neighbourhood are considered land. This way there is no 
    # worrying about interpolating variables on the boundary.
    mask_common[mask_common < 1] = 0

    print 'Setting up ' + out_file
    id = Dataset(out_file, 'w')
    id.createDimension('longitude', size(lon_common))
    id.createDimension('latitude', size(lat_common))
    id.createDimension('time', None)
    id.createVariable('longitude', 'f8', ('longitude'))
    id.variables['longitude'].units = 'degrees'
    id.variables['longitude'][:] = lon_common
    id.createVariable('latitude', 'f8', ('latitude'))
    id.variables['latitude'].units = 'degrees'
    id.variables['latitude'][:] = lat_common
    id.createVariable('time', 'f8', ('time'))
    id.variables['time'].units = 'months'
    id.createVariable('mask', 'f8', ('latitude', 'longitude'))
    id.variables['mask'].units = '1'
    id.variables['mask'][:,:] = mask_common
    id.createVariable('shflux', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['shflux'].long_name = 'surface heat flux into ocean'
    id.variables['shflux'].units = 'W/m^2'
    id.createVariable('ssflux', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['ssflux'].long_name = 'surface virtual salinity flux into ocean'
    id.variables['ssflux'].units = 'psu m/s'
    id.createVariable('sustr', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['sustr'].long_name = 'zonal surface stress'
    id.variables['sustr'].units = 'N/m^2'
    id.createVariable('svstr', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['svstr'].long_name = 'meridional surface stress'
    id.variables['svstr'].units = 'N/m^2'
    id.createVariable('curl_str', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['curl_str'].long_name = 'curl of surface stress'
    id.variables['curl_str'].units = 'N/m^3'
    id.createVariable('uice', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['uice'].long_name = 'sea ice velocity eastward'
    id.variables['uice'].units = 'm/s'
    id.createVariable('vice', 'f8', ('time', 'latitude', 'longitude'))
    id.variables['vice'].long_name = 'sea ice velocity northward'
    id.variables['vice'].units = 'm/s'

    # Loop over months
    for month in range(num_months):
        print 'Processing month ' + str(month+1) + ' of ' + str(num_months)
        # Write time value for this month
        id.variables['time'][month] = month+1

        print '...surface heat flux'
        # Get monthly average
        shflux_roms = monthly_avg_roms(roms_file, 'shflux', shape(lon_rho), month%12+1, instance=month/12+1)
        # Interpolate to common grid
        shflux_common = interp_roms2common(lon_common, lat_common, lon_rho, lat_rho, shflux_roms)
        # Apply land mask
        shflux = ma.masked_where(mask==0, shflux_common)
        # Write to file
        id.variables['shflux'][month,:,:] = shflux

        print '...surface salt flux'
        # Get monthly average
        ssflux_roms = monthly_avg_roms(roms_file, 'ssflux', shape(lon_rho), month%12+1, instance=month/12+1)
        # Interpolate to common grid
        ssflux_common = interp_roms2common(lon_common, lat_common, lon_rho, lat_rho, ssflux_roms)
        # Apply land mask
        ssflux = ma.masked_where(mask==0, ssflux_common)
        # Write to file
        id.variables['ssflux'][month,:,:] = ssflux

        print '...surface stress vector'
        # Surface stresses
        # Get monthly averages of both vector components
        sustr_tmp = monthly_avg_roms(roms_file, 'sustr', u_shape, month%12+1, instance=month/12+1)
        svstr_tmp = monthly_avg_roms(roms_file, 'svstr', v_shape, month%12+1, instance=month/12+1)
        # Rotate to lon-lat space (note they are on the rho grid now)
        sustr_roms, svstr_roms = rotate_vector_roms(sustr_tmp, svstr_tmp, angle_roms)
        # Interpolate to common grid
        sustr_common = interp_roms2common(lon_common, lat_common, lon_rho, lat_rho, sustr_roms)
        svstr_common = interp_roms2common(lon_common, lat_common, lon_rho, lat_rho, svstr_roms)
        # Apply land mask
        sustr = ma.masked_where(mask==0, sustr_common)
        svstr = ma.masked_where(mask==0, svstr_common)
        # Write to file
        id.variables['sustr'][month,:,:] = sustr
        id.variables['svstr'][month,:,:] = svstr

        print '...curl of surface stress vector'
        # Curl of surface stress = d/dx (svstr) - d/dy (sustr)
        # First calculate the two derivatives
        dsvstr_dx = ma.empty(shape(svstr_common))
        # Forward difference approximation
        dsvstr_dx[:,:-1] = (svstr_common[:,1:] - svstr_common[:,:-1])/(2*dx[:,:-1])
        # Backward difference for the last row
        dsvstr_dx[:,-1] = (svstr_common[:,-1] - svstr_common[:,-2])/(2*dx[:,-1])
        dsustr_dy = ma.empty(shape(sustr_common))
        dsustr_dy[:-1,:] = (sustr_common[1:,:] - sustr_common[:-1,:])/(2*dy[:-1,:])
        dsustr_dy[-1,:] = (sustr_common[-1,:] - sustr_common[-2,:])/(2*dy[-1,:])
        curl_str = dsvstr_dx - dsustr_dy
        # Write to file
        id.variables['curl_str'][month,:,:] = curl_str

        print '...sea ice velocity vector'
        # Sea ice velocity (CICE variable not ROMS)
        # Get monthly averages of both vector components
        uice_tmp = monthly_avg_cice(cice_file, 'uvel', shape(lon_cice), month%12+1, instance=month/12+1)
        vice_tmp = monthly_avg_cice(cice_file, 'vvel', shape(lon_cice), month%12+1, instance=month/12+1)
        # Rotate to lon-lat space
        uice_cice, vice_cice = rotate_vector_cice(uice_tmp, vice_tmp, angle_cice)
        # Interpolate to common grid (note CICE grid not ROMS)
        uice_common = interp_roms2common(lon_common, lat_common, lon_cice, lat_cice, uice_cice)
        vice_common = interp_roms2common(lon_common, lat_common, lon_cice, lat_cice, vice_cice)
        # Apply land mask
        uice = ma.masked_where(mask==0, uice_common)
        vice = ma.masked_where(mask==0, vice_common)
        # Write to file
        id.variables['uice'][month,:,:] = uice
        id.variables['vice'][month,:,:] = vice

    print 'Finished'
    id.close()        
    

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
    lat_2d, lon_2d = meshgrid(lat_1d, lon_1d)

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
    
