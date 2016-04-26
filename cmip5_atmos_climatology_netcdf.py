from cmip5_paths import *
from cmip5_field import *
from eraint_field import *
from numpy import *
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator

# NB for raijin users: RegularGridInterpolator needs python/2.7.6 but the
# default is 2.7.3. Before running this script, switch them as follows:
# module unload python/2.7.3
# module unload python/2.7.3-matplotlib
# module load python/2.7.6
# module load python/2.7.6-matplotlib

# For the given CMIP5 model, calculates the monthly climatology from 1992-2005
# inclusive for 11 atmospheric variables. Interpolates to the ERA-Interim grid
# and saves to a NetCDF file. Note that to run this script, you must have
# previously run eraint_climatology_netcdf.py.
# Input:
# model_name = name of model (this must match the list in cmip5_paths.py)
def cmip5_atmos_climatology_netcdf (model_name):

    # Experiment name
    expt = 'historical'
    # Years over which to calculate climatology
    start_year = 1992
    end_year = 2005
    # CMIP5 variable names
    var_names = ['ps', 'tas', 'huss', 'clt', 'uas', 'vas', 'pr', 'prsn', 'evspsbl', 'rsds', 'rlds']
    # Variable names to use in NetCDF file
    var_names_output = ['Pair', 'Tair', 'Hair', 'cloud', 'Uwind', 'Vwind', 'precip', 'snow', 'evap', 'swrad', 'lwrad']
    # Units of final variables (note some conversions in cmip5_field)
    var_units = ['kPa', 'degC', '1', '%', 'm/s', 'm/s', '10^6 kg/m^2/s', '10^6 kg/m^2/s', '10^6 kg/m^2/s', 'W/m^2', 'W/m^2']
    # Path to output NetCDF file
    output_file = '/short/y99/kaa561/CMIP5_forcing/atmos/' + model_name + '.nc'
    # Path to corresponding ERA-Interim file (created using
    # eraint_climatology_netcdf.py)
    eraint_file = '/short/y99/kaa561/CMIP5_forcing/atmos/ERA-Interim.nc'

    # Read ERA-Interim grid
    id = Dataset(eraint_file, 'r')
    era_lon = id.variables['longitude'][:]
    era_lat = id.variables['latitude'][:]
    id.close()

    # Build a list of Model objects for CMIP5
    all_models = build_model_list()
    # Build a corresponding list of all CMIP5 model names
    all_model_names = []
    for model in all_models:
        all_model_names.append(model.name)
    # Find index of model_name in all_model_names, and select the Model object
    # at the same index of all_models
    # Now model_name has a corresponding Model object
    model = all_models[all_model_names.index(model_name)]

    # Loop over variables
    for i in range(len(var_names)):
        var = var_names[i]
        print 'Processing variable ' + var

        # Read monthly climatology for this variable
        model_data, model_lon, model_lat, tmp = cmip5_field(model, expt, var, start_year, end_year)
        # Set up array for climatology interpolated to ERA-Interim grid
        model_data_interp = ma.empty([12, size(era_lat), size(era_lon)])
        if model_data is not None:
            # Interpolate one month at a time
            for t in range(size(model_data,0)):
                model_data_interp[t,:,:] = interp_model2era(model_data[t,:,:], model_lon, model_lat, era_lon, era_lat)
        else:
            # No data (missing variable in CMIP5 archive)
            model_data_interp[:,:,:] = ma.masked
        if i == 0:
            # Set up NetCDF file on the first iteration
            print 'Setting up ' + output_file
            id = Dataset(output_file, 'w')
            # Define dimensions
            id.createDimension('longitude', size(era_lon))
            id.createDimension('latitude', size(era_lat))
            id.createDimension('time', 12)
            # Define dimension variables and fill with data
            id.createVariable('longitude', 'f8', ('longitude'))
            id.variables['longitude'].units = 'degrees'
            id.variables['longitude'][:] = era_lon
            id.createVariable('latitude', 'f8', ('latitude'))
            id.variables['latitude'].units = 'degrees'
            id.variables['latitude'][:] = era_lat
            id.createVariable('time', 'f8', ('time'))
            id.variables['time'].units = 'month'
            id.variables['time'][:] = arange(1, 12+1)
        # Define new variable and fill with data
        id.createVariable(var_names_output[i], 'f8', ('time', 'latitude', 'longitude'))
        id.variables[var_names_output[i]].units = var_units[i]
        id.variables[var_names_output[i]][:,:,:] = model_data_interp
    id.close()


# Interpolate a given CMIP5 field to the ERA-Interim grid.
# Input:
# model_data = 2D array (size mxn) of model data for the given variable
# model_lon = 1D array (length n) of longitude for this model
# model_lat = 1D array (length m) of latitude for this model
# era_lon = 1D array (length q) of longitude on the ERA-Interim grid
# era_lat = 1D array (length p) of latitude on the ERA-Interim grid
# Output:
# data_interp = 2D array (size pxq) of model data interpolated to the
#               ERA-Interim grid 
def interp_model2era(model_data, model_lon, model_lat, era_lon, era_lat):

    # Make sure the model's longitude goes from 0 to 360, not -180 to 180
    index = model_lon < 0
    model_lon[index] = model_lon[index] + 360

    # CMIP5 model axes don't wrap around; there is a gap between almost-180W
    # and almost-180E (these are 0 to 360 but you get the idea) and depending
    # on the grid, we may need to interpolate in this gap.
    # So copy the last longitude value (mod 360) to the beginning, and the
    # first longitude value (mod 360) to the end.
    model_lon_wrap = zeros(size(model_lon)+2)
    model_lon_wrap[0] = model_lon[-1] - 360
    model_lon_wrap[1:-1] = model_lon
    model_lon_wrap[-1] = model_lon[0] + 360
    model_lon = model_lon_wrap
    # Copy the westernmost and easternmost data points to match
    model_data_wrap = ma.array(zeros((size(model_lat), size(model_lon))))
    model_data_wrap[:,1:-1] = model_data
    model_data_wrap[:,0] = model_data[:,-1]
    model_data_wrap[:,-1] = model_data[:,0]
    model_data = model_data_wrap

    if amin(model_lat) > amin(era_lat):
        # Add a point at 90S
        model_lat_new = zeros(size(model_lat)+1)
        model_data_new = ma.array(zeros((size(model_lat_new), size(model_lon))))
        if model_lat[0] > model_lat[1]:
            model_lat_new[0:-1] = model_lat
            model_lat_new[-1] = -90.0
            model_data_new[0:-1,:] = model_data
            model_data_new[-1,:] = model_data[-1,:]
        elif model_lat[0] < model_lat[1]:
            model_lat_new[1:] = model_lat
            model_lat_new[0] = -90.0
            model_data_new[1:,:] = model_data
            model_data_new[0,:] = model_data[0,:]
        model_lat = model_lat_new
        model_data = model_data_new
    if amax(model_lat) < amax(era_lat):
        # Add a point at 90N
        model_lat_new = zeros(size(model_lat)+1)
        model_data_new = ma.array(zeros((size(model_lat_new), size(model_lon))))
        if model_lat[0] > model_lat[1]:
            model_lat_new[1:] = model_lat
            model_lat_new[0] = 90.0
            model_data_new[1:,:] = model_data
            model_data_new[0,:] = model_data[0,:]
        elif model_lat[0] < model_lat[1]:
            model_lat_new[0:-1] = model_lat
            model_lat_new[-1] = 90.0
            model_data_new[0:-1,:] = model_data
            model_data_new[-1,:] = model_data[-1,:]
        model_lat = model_lat_new
        model_data = model_data_new

    # Get 2D mesh of ERA-Interim lat and lon
    era_lon_2d, era_lat_2d = meshgrid(era_lon, era_lat)
    # Build an interpolation function for model_data
    interp_function = RegularGridInterpolator((model_lat, model_lon), model_data)
    # Call it for the ERA-Interim grid
    data_interp = interp_function((era_lat_2d, era_lon_2d))

    return data_interp


# Command-line interface
if __name__ == "__main__":

    # Process one model at a time
    models = build_model_list()
    for model in models:
        print model.name
        cmip5_atmos_climatology_netcdf(model.name)
            
