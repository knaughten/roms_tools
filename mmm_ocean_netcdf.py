from cmip5_paths import *
from numpy import *
from netCDF4 import Dataset

# Calculate the multi-model mean of ocean climatology files created using
# cmip5_ocean_climatology_netcdf.py.
def mmm_ocean_netcdf ():

    # Directory containing CMIP5 NetCDF climatology ocean files, interpolated
    # to the ECCO2 grid at the northern boundary of ROMS (created using
    # cmip5_ocean_climatology_netcdf.py)
    directory = '/short/y99/kaa561/CMIP5_forcing/ocean/'
    # Path to ECCO2 file (created using cmip5_ocean_climatology_netcdf.py)
    ecco2_file = directory + 'ECCO2.nc'
    # Path to multi-model mean output file
    output_file = directory + 'MMM.nc'
    # Variable names in NetCDF files
    var_names = ['temp', 'salt', 'u', 'v']
    # Corresponding units
    var_units = ['degC', 'psu', 'm/s', 'm/s']

    # Read ECCO2 grid
    id = Dataset(ecco2_file, 'r')
    lon = id.variables['longitude'][:]
    depth = id.variables['depth'][:]
    id.close()

    # Set up output file
    print 'Setting up ' + output_file
    out_id = Dataset(output_file, 'w')
    # Define dimensions
    out_id.createDimension('longitude', size(lon))
    out_id.createDimension('depth', size(depth))
    out_id.createDimension('time', 12)
    # Define dimension variables and fill with axes
    out_id.createVariable('longitude', 'f8', ('longitude'))
    out_id.variables['longitude'].units = 'degrees'
    out_id.variables['longitude'][:] = lon
    out_id.createVariable('depth', 'f8', ('depth'))
    out_id.variables['depth'].units = 'metres'
    out_id.variables['depth'][:] = depth
    out_id.createVariable('time', 'f8', ('time'))
    out_id.variables['time'].units = 'month'
    out_id.variables['time'][:] = arange(1, 12+1)

    # Get a list of CMIP5 Model objects
    models = build_model_list()
    # Build a corresponding list of model names
    model_names = []
    for model in models:
        model_names.append(model.name)

    # Loop over variables
    for i in range(len(var_names)):
        var = var_names[i]
        print 'Variable ' + var

        num_models = 0  # Number of models in the multi-model mean
        multi_model_mean = None
        # Loop over models
        for model_name in model_names:
            # Read model data
            id = Dataset(directory + model_name + '.nc', 'r')
            model_data = id.variables[var][:,:,:]
            id.close()

            if all(model_data.mask):
                # Everything is masked; data is missing for this variable
                pass
            else:
                # Add to multi-model mean and increment num_models
                if multi_model_mean is None:
                    multi_model_mean = model_data[:,:,:]
                else:
                    multi_model_mean[:,:,:] += model_data[:,:,:]
                num_models += 1

        # Divide multi_model_mean (currently just a sum) by num_models
        multi_model_mean /= num_models
        # Define variable and fill with data
        out_id.createVariable(var_names[i], 'f8', ('time', 'depth', 'longitude'))
        out_id.variables[var_names[i]].units = var_units[i]
        out_id.variables[var_names[i]][:,:,:] = multi_model_mean
    out_id.close()


# Command-line interface
if __name__ == "__main__":

    mmm_ocean_netcdf()
    


    
