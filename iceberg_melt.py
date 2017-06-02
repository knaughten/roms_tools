from netCDF4 import Dataset
from numpy import *
from scipy.interpolate import griddata

# Read Martin and Adcroft's monthly climatology of freshwater fluxes
# from iceberg melt, interpolate to the ROMS grid, and save as a
# forcing file.
# Input:
# out_file = path to desired output file
def iceberg_melt (out_file):

    # Naming conventions for iceberg files
    iceberg_head = '/short/m68/kaa561/metroms_iceshelf/data/originals/MartinAdcroft2010_iceberg_meltfluxes/icebergs.1861-1960.'
    iceberg_tail = '.melt.nc'
    # Iceberg grid file
    iceberg_grid = '/short/m68/kaa561/metroms_iceshelf/data/originals/MartinAdcroft2010_iceberg_meltfluxes/icebergs.static.nc'
    # ROMS grid file
    roms_grid ='/short/m68/kaa561/metroms_iceshelf/apps/common/grid/circ30S_quarterdegree.nc'

    # Read ROMS grid
    id = Dataset(roms_grid, 'r')
    lon_roms = id.variables['lon_rho'][:,:]
    lat_roms = id.variables['lat_rho'][:,:]
    id.close()
    num_lon = size(lon_roms, 1)
    num_lat = size(lon_roms, 0)

    # Read the iceberg grid
    id = Dataset(iceberg_grid, 'r')
    lon_iceberg = id.variables['lon'][:,:]
    lat_iceberg = id.variables['lat'][:,:]
    id.close()

    # Make sure longitudes are between 0 and 360
    index = lon_iceberg < 0
    lon_iceberg[index] = lon_iceberg[index] + 360

    # Set up output file
    out_id = Dataset(out_file, 'w')
    # Define dimensions
    out_id.createDimension('xi_rho', num_lon)
    out_id.createDimension('eta_rho', num_lat)
    out_id.createDimension('time', None)
    # Define variables
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
    out_id.variables['time'].cycle_length = 365.25
    out_id.createVariable('icebergs', 'f8', ('time', 'eta_rho', 'xi_rho'))
    out_id.variables['icebergs'].long_name = 'freshwater flux from iceberg melt'
    out_id.variables['icebergs'].units = 'kg/m^2/s'

    # Loop over months
    for month in range(12):
        print 'Processing month ' + str(month+1)
        # Reconstruct the filename of this month's iceberg data
        if month+1 < 10:
            month_str = '0' + str(month+1)
        else:
            month_str = str(month+1)
        iceberg_file = iceberg_head + month_str + iceberg_tail
        # Read iceberg freshwater flux in kg/m^2/s
        id = Dataset(iceberg_file, 'r')
        melt_iceberg = id.variables['melt'][0,:,:]
        id.close()
        # Interpolate to ROMS grid
        melt_roms = interp_iceberg2roms(melt_iceberg, lon_iceberg, lat_iceberg, lon_roms, lat_roms)
        # Calculate time values centered in the middle of each month
        out_id.variables['time'][month] = 365.25/12*(month+0.5)
        # Save data
        out_id.variables['icebergs'][month,:,:] = melt_roms
    out_id.close()


# Given a field A on the iceberg grid, linearly interpolate to the ROMS grid.
# Input:
# A = 2D array (m x n) containing data on the iceberg grid
# lon_iceberg = 2D array (m x n) contaning longitude values on the
#               iceberg grid, from 0 to 360
# lat_iceberg = 2D array (m x n) containing latitude values on the
#               iceberg grid
# lon_roms = 2D array (p x q) containing longitude values on the ROMS grid,
#            from 0 to 360
# lat_roms = 2D array (p x q) containing latitude values on the ROMS grid
# Output:
# A_interp = 2D array (p x q) containing A interpolated to the ROMS grid
def interp_iceberg2roms (A, lon_iceberg, lat_iceberg, lon_roms, lat_roms):

    # Set up an nx2 array containing the coordinates of each point in the
    # iceberg grid
    points = empty([size(lon_iceberg), 2])
    points[:,0] = ravel(lon_iceberg)
    points[:,1] = ravel(lat_iceberg)
    # Also flatten the data
    values = ravel(A)
    # Now set up an mx2 array containing the coordinates of each point
    # we want to interpolate to, in the ROMS grid
    xi = empty([size(lon_roms), 2])
    xi[:,0] = ravel(lon_roms)
    xi[:,1] = ravel(lat_roms)
    # Now call griddata; fill out-of-bounds values (such as under ice shelves)
    # with zeros
    result = griddata(points, values, xi, method='linear', fill_value=0.0)
    # Un-flatten the result
    A_interp = reshape(result, shape(lon_roms))
    # Enforce periodic boundary
    A_interp[:,0] = A_interp[:,-2]
    A_interp[:,-1] = A_interp[:,1]

    return A_interp


if __name__ == "__main__":

    file = raw_input('Path to desired output file: ')
    iceberg_melt(file)

    

        

        

    
