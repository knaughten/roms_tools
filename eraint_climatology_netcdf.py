from eraint_field import *
from numpy import *
from netCDF4 import Dataset

# Calculate the ERA-Interim monthly climatology from 1992-2005 inclusive, for 11
# atmospheric variables (all the variables which ROMS and/or FESOM depend on).
# Save to a NetCDF file.
def eraint_climatology_netcdf ():

    # Date range for climatology
    start_year = 1992
    end_year = 2005
    # Names of ERA-Interim variables
    var_names = ['sp', 't2m', 'd2m', 'tcc', 'u10', 'v10', 'tp', 'sf', 'e', 'ssrd', 'strd']
    # Variable names to use for NetCDF file
    var_names_output = ['Pair', 'Tair', 'Hair', 'cloud', 'Uwind', 'Vwind', 'precip', 'snow', 'evap', 'swrad', 'lwrad']
    # Units of final variables (note there are some conversions in eraint_field)
    var_units = ['kPa', 'degC', '1', '%', 'm/s', 'm/s', '10^6 kg/m^2/s', '10^6 kg/m^2/s', '10^6 kg/m^2/s', 'W/m^2', 'W/m^2']
    # Path to output NetCDF file
    output_file = '/short/y99/kaa561/CMIP5_forcing/atmos/ERA-Interim.nc'

    # Loop over variables
    for i in range(len(var_names)):
        var = var_names[i]
        print 'Processing variable ' + var
        # Read data and grid
        data, lon, lat = eraint_field(var, start_year, end_year)
        if i == 0:
            # Create NetCDF file on the first iteration
            print 'Setting up ' + output_file
            id = Dataset(output_file, 'w')
            # Define dimensions
            id.createDimension('longitude', size(lon))
            id.createDimension('latitude', size(lat))
            id.createDimension('time', 12)
            # Define dimension variables and fill with data
            id.createVariable('longitude', 'f8', ('longitude'))
            id.variables['longitude'].units = 'degrees'
            id.variables['longitude'][:] = lon
            id.createVariable('latitude', 'f8', ('latitude'))
            id.variables['latitude'].units = 'degrees'
            id.variables['latitude'][:] = lat
            id.createVariable('time', 'f8', ('time'))
            id.variables['time'].units = 'month'
            id.variables['time'][:] = arange(1, 12+1)
        # Define a new variable in the NetCDF file and fill with data
        id.createVariable(var_names_output[i], 'f8', ('time', 'latitude', 'longitude'))
        id.variables[var_names_output[i]].units = var_units[i]
        id.variables[var_names_output[i]][:,:,:] = data
    id.close()


# Command-line interface
if __name__ == "__main__":

    eraint_climatology_netcdf()
