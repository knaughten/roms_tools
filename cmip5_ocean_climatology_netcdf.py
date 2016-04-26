from cmip5_paths import *
from cmip5_field import *
from ecco2_field import *
from numpy import *
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator, griddata
from scipy.spatial import KDTree

# NB for raijin users: RegularGridInterpolator needs python/2.7.6 but the
# default is 2.7.3. Before running this script, switch them as follows:
# module unload python/2.7.3
# module unload python/2.7.3-matplotlib
# module load python/2.7.6
# module load python/2.7.6-matplotlib

# For the given CMIP5 model, calculates the monthly climatology from 1992-2005
# inclusive for 4 ocean variables. Interpolates to the ECCO2 grid at the
# northern boundary of ROMS (currently 30S, see ecco2_climatology_netcdf.py)
# and saves to a NetCDF file. Note that to run this script, you must have
# previously run ecco2_climatology_netcdf.py.
# Input:
# model_name = name of model (this must match the list in cmip5_paths.py)
def cmip5_ocean_climatology_netcdf (model_name):

    # Experiment name
    expt = 'historical'
    # Years over which to calculate climatology
    start_year = 1992
    end_year = 2005
    # Northern boundary of ROMS
    nbdry = -30.0
    # CMIP5 variable names
    var_names = ['thetao', 'so', 'uo', 'vo']
    # Variable names to use in NetCDF file
    var_names_output = ['temp', 'salt', 'u', 'v']
    # Units of final variables (note some conversions in cmip5_field)
    var_units = ['degC', 'psu', 'm/s', 'm/s']
    # Path to output NetCDF file
    output_file = '/short/y99/kaa561/CMIP5_forcing/ocean/' + model_name + '.nc'
    # Path to corresponding ECCO2 file (created using
    # ecco2_climatology_netcdf.py)
    ecco2_file = '/short/y99/kaa561/CMIP5_forcing/ocean/ECCO2.nc'

    # Read ECCO2 grid
    id = Dataset(ecco2_file, 'r')
    ecco_lon = id.variables['longitude'][:]
    ecco_depth = id.variables['depth'][:]
    # Save the land mask
    mask = id.variables['temp'][0,:,:].mask
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
        model_data, model_lon, model_lat, model_depth = cmip5_field(model, expt, var, start_year, end_year)
        # Set up array for climatology interpolated to ECCO2 grid
        model_data_interp = ma.empty([12, size(ecco_depth), size(ecco_lon)])
        if model_data is not None:
            # Interpolate one month at a time
            for t in range(size(model_data,0)):
                tmp = interp_model2ecco(model_data[t,:,:,:], model_lon, model_lat, model_depth, ecco_lon, nbdry, ecco_depth)
                # Mask with ECCO2 land mask
                model_data_interp[t,:,:] = ma.masked_where(mask, tmp)
        else:
            # No data (missing variable in CMIP5 archive)
            model_data_interp[:,:,:] = ma.masked
        if i == 0:
            # Set up NetCDF file on the first iteration
            print 'Setting up ' + output_file
            id = Dataset(output_file, 'w')
            # Define dimensions
            id.createDimension('longitude', size(ecco_lon))
            id.createDimension('depth', size(ecco_depth))
            id.createDimension('time', 12)
            # Define dimension variables and fill with data
            id.createVariable('longitude', 'f8', ('longitude'))
            id.variables['longitude'].units = 'degrees'
            id.variables['longitude'][:] = ecco_lon
            id.createVariable('depth', 'f8', ('depth'))
            id.variables['depth'].units = 'metres'
            id.variables['depth'][:] = ecco_depth
            id.createVariable('time', 'i4', ('time'))
            id.variables['time'].units = 'month'
            id.variables['time'][:] = arange(1, 12+1)
        # Define new variable and fill with data
        id.createVariable(var_names_output[i], 'f8', ('time', 'depth', 'longitude'))
        id.variables[var_names_output[i]].units = var_units[i]
        id.variables[var_names_output[i]][:,:,:] = model_data_interp
    id.close()


# Interpolate a given CMIP5 field to the ECCO2 grid at the northern boundary
# of ROMS.
# model_data = 3D array (size mxnxo) of model data for the given variable
# model_lon = 1D array (length o) or 2D array (size nxo) of longitude for this
#             model
# model_lat = 1D array (length n) or 2D array (size nxo) of latitude for this
#             model
# model_depth = 1D array (length m) or 3D array (size mxnxo) of depth for this
#               model, in z-coordinates
# ecco_lon = 1D array (length q) of longitude on the ECCO2 grid
# nbdry = latitude of the northern boundary of ROMS
# ecco_depth = 1D array (length p) of depth on the ECCO2 grid
# Output:
# data_interp = 2D array (size pxq) of model data interpolated to the ECCO2
#               grid at the northern boundary of ROMS
def interp_model2ecco (model_data, model_lon, model_lat, model_depth, ecco_lon, nbdry, ecco_depth):

    # Make sure the model's longitude goes from 0 to 360, not -180 to 180
    index = model_lon < 0
    model_lon[index] = model_lon[index] + 360

    # Check for model depth axis which starts at the bottom
    if len(model_depth.shape) == 1:
        if model_depth[0] > model_depth[1]:
            # Flip it, and the data, so it starts at the surface
            model_depth = flipud(model_depth)
            model_data = flipud(model_data)    

    # Interpolate to northern boundary

    if len(model_lon.shape) == 1 and len(model_lat.shape) == 1:
        # Case 1: longitude and latitude are 1D axis variables

        # Find jS and jN, the indices immediately south and immediately north
        # of nbdry
        if model_lat[0] > model_lat[1]:
            # Latitude is decreasing
            # Find the first index south of nbdry
            jS = nonzero(model_lat < nbdry)[0][0]
            # The index before it is the last index north of nbdry
            jN = jS - 1
        elif model_lat[0] < model_lat[1]:
            # Latitude is increasing
            # Find the first index north of nbdry
            jN = nonzero(model_lat > nbdry)[0][0]
            # The index before it is the last index south of nbdry
            jS = jN - 1
        # Linearly interpolate model_data to nbdry
        model_dataS = model_data[:,jS,:]
        model_dataN = model_data[:,jN,:]
        model_data_bdry = (model_dataN - model_dataS)/(model_lat[jN] - model_lat[jS])*(nbdry - model_lat[jS]) + model_dataS

    elif len(model_lon.shape) == 2 and len(model_lat.shape) == 2:
        # Case 2: longitude and latitude are 2D variables, i.e. rotated grid

        # Set up empty arrays for data and longitude interpolated to nbdry
        model_data_bdry = ma.empty((size(model_data,0), size(model_data,2)))
        model_lon_bdry = ma.empty((size(model_data,2)))
        if len(model_depth.shape) == 3:
            # We will also have to interpolate depth to nbdry
            model_depth_bdry = ma.empty(shape(model_data_bdry))
        # Loop over longitude axis and interpolate one index at a time
        for i in range(size(model_data,2)):
            # Select current latitude values
            model_lat_tmp = model_lat[:,i]
            # Check for non-monotonic latitude values
            dlat = model_lat_tmp[1:] - model_lat_tmp[0:-1]
            if any(dlat < 0) and any(dlat > 0):
                # Latitude doesn't strictly increase or strictly decrease with j
                # A few models have rotated grids where latitude is monotonic
                # in the middle but changes in the other direction near the ends
                # We will just mask out the ends, one index at a time, until
                # the leftover values in the middle are monotonic
                # This works because the non-monotonic sections are near the
                # poles, and far away from nbdry (currently 30S)
                model_lat_tmp = ma.masked_array(model_lat_tmp)
                posn = 1
                while True:
                    model_lat_tmp[posn-1] = ma.masked
                    model_lat_tmp[-posn] = ma.masked
                    dlat = model_lat_tmp[1:] - model_lat_tmp[0:-1]
                    if all(dlat < 0) or all(dlat > 0):
                        break
                    posn += 1
            else:
                posn = 0
            if model_lat_tmp[posn] > model_lat_tmp[posn+1]:
                # Latitude is decreasing
                # Find the first index south of nbdry
                jS = nonzero(model_lat_tmp < nbdry)[0][0]
                # The index before it is the last index north of nbdry
                jN = jS - 1
            elif model_lat_tmp[posn] < model_lat_tmp[posn+1]:
                # Latitude is increasing
                # Find the first index north of nbdry
                jN = nonzero(model_lat_tmp > nbdry)[0][0]
                # The index before it is the last index south of nbdry
                jS = jN - 1
            # Linearly interpolate model_data to nbdry
            model_dataS = model_data[:,jS,i]
            model_dataN = model_data[:,jN,i]
            model_data_bdry[:,i] = (model_dataN - model_dataS)/(model_lat[jN,i] - model_lat[jS,i])*(nbdry - model_lat[jS,i]) + model_dataS
            # Linearly interpolate model_lon to nbdry
            model_lonS = model_lon[jS,i]
            model_lonN = model_lon[jN,i]
            if model_lonS < 1 and model_lonN > 300:
                # The jump in longitude from almost-360 to almost-0 happens
                # between lonS and lonN. Take mod 360 to fix this.
                model_lonS += 360
            model_lon_bdry[i] = (model_lonN - model_lonS)/(model_lat[jN,i] - model_lat[jS,i])*(nbdry - model_lat[jS,i]) + model_lonS
            if len(model_depth.shape) == 3:
                # Linearly interpolate depth to nbdry
                model_depthS = model_depth[:,jS,i]
                model_depthN = model_depth[:,jN,i]
                model_depth_bdry[:,i] = (model_depthN - model_depthS)/(model_lat[jN,i] - model_lat[jS,i])*(nbdry - model_lat[jS,i]) + model_depthS
        # Overwrite model_lon with the interpolated values; now it's 1D
        model_lon = model_lon_bdry
        if len(model_depth.shape) == 3:
            # Overwrite model_depth with the interpolated values; now it's 2D
            model_depth = model_depth_bdry

    if model_lon[-2] == model_lon[0] and model_lon[-1] == model_lon[1]:
        # A few models have grids which overlap by 2 indices at the periodic
        # boundary; delete these last 2 indices
        model_lon = model_lon[0:-2]
        model_data_bdry = model_data_bdry[:,0:-2]

    # Most grids will have longitude jump from almost-360 to almost-0
    # at some point in the interior. Find the index where this happens, and
    # rearrange the model data and axes so the jump is moved to the periodic
    # boundary. This way model_lon will be monotonic, which is needed for
    # interpolation.
    dlon = model_lon[1:] - model_lon[0:-1]
    index = argmax(abs(dlon))
    if dlon[index] < -300:
        # Add 1 to index because dlon is one index shorter than model_lon
        jump = index+1
        # Rearrange longitude values and data
        model_lon_new = zeros(size(model_lon))
        model_lon_new[0:-jump] = model_lon[jump:]
        model_lon_new[-jump:] = model_lon[0:jump]
        model_lon = model_lon_new
        model_data_new = ma.empty(shape(model_data_bdry))
        model_data_new[:,0:-jump] = model_data_bdry[:,jump:]
        model_data_new[:,-jump:] = model_data_bdry[:,0:jump]
        model_data_bdry = model_data_new
        if len(model_depth.shape) == 2:
            # Rearrange depth values
            model_depth_new = ma.empty(shape(model_depth))
            model_depth_new[:,0:-jump] = model_depth[:,jump:]
            model_depth_new[:,-jump:] = model_depth[:,0:jump]
            model_depth = model_depth_new        

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
    model_data_wrap = ma.empty((size(model_data_bdry,0), size(model_lon)))
    model_data_wrap[:,1:-1] = model_data_bdry
    model_data_wrap[:,0] = model_data_bdry[:,-1]
    model_data_wrap[:,-1] = model_data_bdry[:,0]
    model_data_bdry = model_data_wrap
    if len(model_depth.shape) == 2:
        # Copy the westernmost and easternmost depth values to match
        model_depth_wrap = ma.empty((size(model_depth,0), size(model_lon)))
        model_depth_wrap[:,1:-1] = model_depth
        model_depth_wrap[:,0] = model_depth[:,-1]
        model_depth_wrap[:,-1] = model_depth[:,0]
        model_depth = model_depth_wrap

    if amin(model_depth) >= amin(ecco_depth):
        # The model grid doesn't extend far enough to the surface to be able
        # to interpolate to ECCO2's surface points (5 metres depth).
        # Add another point at 0m, and just copy the data values closest
        # to the surface.
        if len(model_depth.shape) == 1:
            model_depth_new = zeros(size(model_depth)+1)
            model_depth_new[1:] = model_depth
        elif len(model_depth.shape) == 2:
            model_depth_new = zeros(size(model_depth,0)+1, size(model_depth,1))
            model_depth_new[1:,:] = model_depth
        model_depth = model_depth_new
        model_data_new = ma.empty((size(model_data_bdry,0)+1, size(model_lon)))
        model_data_new[1:,:] = model_data_bdry
        model_data_new[0,:] = model_data_new[1,:]
        model_data_bdry = model_data_new
    if amax(model_depth) <= amax(ecco_depth):
        # The model grid doesn't extend deep enough to be able to interpolate
        # to ECCO2's bottom points (almost 6000m depth).
        # Add another point at 6000m, and just copy the data values closest
        # to the bottom.
        if len(model_depth.shape) == 1:        
            model_depth_new = zeros(size(model_depth)+1)
            model_depth_new[0:-1] = model_depth
            model_depth_new[-1] = 6000.0
        elif len(model_depth.shape) == 2:
            model_depth_new = zeros((size(model_depth,0)+1, size(model_depth,1)))
            model_depth_new[0:-1,:] = model_depth
            model_depth_new[-1,:] = 6000.0
        model_depth = model_depth_new
        model_data_new = ma.empty((size(model_data_bdry,0)+1, size(model_lon)))
        model_data_new[0:-1,:] = model_data_bdry
        model_data_new[-1,:] = model_data_new[-2,:]
        model_data_bdry = model_data_new

    # Fill in masked values with nearest neighbours so they don't screw up
    # the interpolation
    # I got this code from Stack Exchange, not really sure how it works
    j,i = mgrid[0:model_data_bdry.shape[0], 0:model_data_bdry.shape[1]]
    jigood = array((j[~model_data_bdry.mask], i[~model_data_bdry.mask])).T
    jibad = array((j[model_data_bdry.mask], i[model_data_bdry.mask])).T
    model_data_bdry[model_data_bdry.mask] = model_data_bdry[~model_data_bdry.mask][KDTree(jigood).query(jibad)[1]]

    # Get a 2D mesh of ECCO2 longitude and depth
    ecco_lon_2d, ecco_depth_2d = meshgrid(ecco_lon, ecco_depth)

    if len(model_depth.shape) == 1:

        # Build an interpolation function for model_data
        interp_function = RegularGridInterpolator((model_depth, model_lon), model_data_bdry)
        # Call it for the ECCO2 grid
        data_interp = interp_function((ecco_depth_2d, ecco_lon_2d))

    elif len(model_depth.shape) == 2:

        # Since the grid isn't regular in lon-depth space, we have to use
        # griddata instead of RegularGridInterpolator
        # First copy the model's longitude values into a 2D array of dimension
        # depth x longitude
        model_lon_2d = tile(model_lon, (size(model_data_bdry,0), 1))
        # Now set up an nx2 array containing the coordinates of each point
        # in the flattened array of model data
        points = empty([size(model_data_bdry),2])
        points[:,0] = ravel(model_lon_2d)
        points[:,1] = ravel(model_depth)
        # Also flatten the model data
        values = ravel(model_data_bdry)
        # Now set up an mx2 array containing the coordinates of each point
        # we want to interpolate to, in the ECCO2 grid
        xi = empty([size(ecco_lon_2d),2])
        xi[:,0] = ravel(ecco_lon_2d)
        xi[:,1] = ravel(ecco_depth_2d)
        # Now call griddata
        result = griddata(points, values, xi)
        # Un-flatten the result
        data_interp = reshape(result, shape(ecco_lon_2d))

    return data_interp


# Command-line interface
if __name__ == "__main__":

    # Process one model at a time
    models = build_model_list()
    for model in models:
        print model.name
        cmip5_ocean_climatology_netcdf(model.name)
