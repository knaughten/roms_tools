from netCDF4 import Dataset
from numpy import *
from scipy.interpolate import RectBivariateSpline

def romscice_gpcp (year):

    grid_file = '/short/m68/kaa561/metroms_iceshelf/apps/common/grid/circ30S_quarterdegree.nc'
    gpcp_head = '/short/m68/kaa561/gpcp/gpcp_cdr_v23rB1_y' + str(year) + '_m'
    output_file = '/short/m68/kaa561/metroms_iceshelf/data/GPCP/precip_' + str(year) + '_monthly.nc'
    conv_factor = 5e-4 # mm/day to m/12h
    days_per_month = array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

    if year % 4 == 0:
        days_per_month[1] = 29

    # Read ROMS latitude and longitude
    id = Dataset(grid_file, 'r')
    lon_roms = id.variables['lon_rho'][:,:]
    lat_roms = id.variables['lat_rho'][:,:]
    id.close()
    num_lon = size(lon_roms, 1)
    num_lat = size(lon_roms, 0)

    # Read GPCP latitude and longitude
    id = Dataset(gpcp_head + '01.nc', 'r')
    lon_gpcp_raw = id.variables['longitude'][:]
    lat_gpcp = id.variables['latitude'][:]
    id.close()

    # Wrap the GPCP longitude axis around 360 so there is no gap
    lon_gpcp = zeros(size(lon_gpcp_raw)+2)
    lon_gpcp[0] = lon_gpcp_raw[-1]-360
    lon_gpcp[1:-1] = lon_gpcp_raw
    lon_gpcp[-1] = lon_gpcp_raw[0]+360

    # Find number of days between 1 Jan 1992 and 1 Jan this year
    time_start = (year-1992)*365.0 + ceil((year-1992)/4.0)

    # Set up output file
    out_id = Dataset(output_file, 'w')
    # Define dimensions
    out_id.createDimension('xi_rho', num_lon)
    out_id.createDimension('eta_rho', num_lat)
    out_id.createDimension('time', None)
    out_id.createVariable('lon_rho', 'f8', ('eta_rho', 'xi_rho'))
    out_id.variables['lon_rho'].long_name = 'longitude of rho-points'
    out_id.variables['lon_rho'].units = 'degree_east'
    out_id.variables['lon_rho'][:,:] = lon_roms
    out_id.createVariable('lat_rho', 'f8', ('eta_rho', 'xi_rho'))
    out_id.variables['lat_rho'].long_name = 'latitude of rho-points'
    out_id.variables['lat_rho'].units = 'degree_north'
    out_id.variables['lat_rho'][:,:] = lat_roms
    out_id.createVariable('time', 'f8', ('time'))
    out_id.variables['time'].units = 'days since 1992-01-01 00:00:0.0'
    out_id.createVariable('rain', 'f8', ('time', 'eta_rho', 'xi_rho'))
    out_id.variables['rain'].long_name = 'rain fall rate'
    out_id.variables['rain'].units = 'm_per_12hr'

    # Loop over months
    for month in range(12):
        # Save output time value: 15th of this month at midnight
        out_id.variables['time'][month] = time_start + sum(days_per_month[0:month]) + 14
        # Construct path to GPCP file for this month
        if month+1 < 10:
            month_str = '0' + str(month+1)
        else:
            month_str = str(month+1)
        in_id = Dataset(gpcp_head + month_str + '.nc', 'r')
        # Read GPCP precip data and convert to m/12h
        precip_raw = in_id.variables['precip'][:,:]*conv_factor
        in_id.close()
        # Wrap around periodic boundary
        precip = zeros([size(precip_raw,0), size(precip_raw,1)+2])
        precip[:,0] = precip_raw[:,-1]
        precip[:,1:-1] = precip_raw
        precip[:,-1] = precip_raw[:,0]
        # Interpolate to ROMS grid
        rain = interp_gpcp2roms(precip, lon_gpcp, lat_gpcp, lon_roms, lat_roms)
        # Make sure there are no negative values
        rain[rain < 0] = 0.0
        # Save
        out_id.variables['rain'][month,:,:] = rain
    out_id.close()


def interp_gpcp2roms (A, lon_gpcp, lat_gpcp, lon_roms, lat_roms):

    # First build a function to approximate A with 2D splines
    interp_function = RectBivariateSpline(lat_gpcp, lon_gpcp, A)
    B = zeros(shape(lon_roms))
    # Call this function for each grid point individually - if you try to do
    # it all at once it throws a MemoryError
    for i in range(size(lon_roms,1)):
        for j in range(size(lon_roms,0)):
            B[j,i] = interp_function(lat_roms[j,i], lon_roms[j,i])

    # Enforce periodic boundary
    B[:,0] = B[:,-2]
    B[:,-1] = B[:,1]

    return B


# Command-line interface
if __name__ == "__main__":

    first_year = 1992
    last_year = 2005
    for year in range(first_year, last_year+1):
        print 'Processing ' + str(year)
        romscice_gpcp(year)
    
    
