from ecco2_field import *
from numpy import *
from netCDF4 import Dataset

# Calculate the ECCO2 monthly climatology from 1992-2005 inclusive, interpolated
# to the northern boundary of the ROMS domain (currently 30S, this is set in
# ecco2_field.py), for 4 ocean variables. Save to a NetCDF file.
def ecco2_climatology_netcdf ():

    # Date range for climatology
    start_year = 1992
    end_year = 2005
    # Names of ECCo2 variables
    var_names = ['THETA', 'SALT', 'UVEL', 'VVEL']
    # Variable names to use for NetCDF file
    var_names_output = ['temp', 'salt', 'u', 'v']
    # Units of final variables (note there are some conversions in ecco2_field)
    var_units = ['degC', 'psu', 'm/s', 'm/s']
    # Path to output NetCDF file
    output_file = '/short/y99/kaa561/CMIP5_forcing/ocean/ECCO2.nc'

    # Loop over variables
    for i in range(len(var_names)):
        var = var_names[i]
        print 'Processing variable ' + var
        # Read data and grid
        data, lon, depth = ecco2_field(var, start_year, end_year)
        if i == 0:
            # Save mask of the first variable (temperature)
            orig_mask = data.mask
            # Create NetCDF file on the first iteration
            print 'Setting up ' + output_file
            id = Dataset(output_file, 'w')
            # Define dimensions
            id.createDimension('longitude', size(lon))
            id.createDimension('depth', size(depth))
            id.createDimension('time', 12)
            # Define dimension variables and fill with data
            id.createVariable('longitude', 'f8', ('longitude'))
            id.variables['longitude'].units = 'degrees'
            id.variables['longitude'][:] = lon
            id.createVariable('depth', 'f8', ('depth'))
            id.variables['depth'].units = 'metres'
            id.variables['depth'][:] = depth
            id.createVariable('time', 'i4', ('time'))
            id.variables['time'].units = 'month'
            id.variables['time'][:] = arange(1, 12+1)            
        if all(~data.mask):
            # There is no mask on this data
            # This happens for u and v which are filled with zeros in the mask
            # Apply the mask from temperature, saved earlier
            data = ma.masked_where(orig_mask, data)
        # Define a new variable in the NetCDF file and fill with data
        id.createVariable(var_names_output[i], 'f8', ('time', 'depth', 'longitude'))
        id.variables[var_names_output[i]].units = var_units[i]
        id.variables[var_names_output[i]][:,:,:] = data
    id.close()


# Command-line interface
if __name__ == "__main__":

    ecco2_climatology_netcdf()

    

        
