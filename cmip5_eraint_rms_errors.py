from cmip5_paths import *
from cmip5_field import *
from eraint_field import *
from ecco2_field import *
from numpy import *
from scipy.interpolate import RegularGridInterpolator
from matplotlib.pyplot import *

# NB for raijin users: RegularGridInterpolator needs python/2.7.6 but the
# default is 2.7.3. Before running this script, switch them as follows:
# module unload python/2.7.3
# module unload python/2.7.3-matplotlib
# module load python/2.7.6
# module load python/2.7.6-matplotlib

# Calculate root-mean-square errors (as in Gleckler et al., 2008) for each of
# 39 CMIP5 models and the multi-model mean, with respect to 11 ERA-Interim
# variables (all variables which ROMS and/or FESOM require as forcing). 
# The domain is the Southern Ocean (all longitudes, and latitudes from the
# northern boundary of ROMS to the southernmost ocean point not in an ice shelf
# cavity) and the monthly climatology averaged over 1992-2005 inclusive.
# Also calculate the relative errors as in Gleckler et al. and make a "portrait
# plot" of coloured tiles in a model vs variable matrix.
# Save both rms errors and relative errors into text files.
def cmip5_eraint_rms_errors():

    # Experiment name
    expt = 'historical'
    # Years to average over
    start_year = 1992
    end_year = 2005
    # CMIP5 atmospheric variables to read
    var_names_cmip5 = ['ps', 'tas', 'huss', 'clt', 'uas', 'vas', 'pr', 'prsn', 'evspsbl', 'rsds', 'rlds']
    # Corresponding ERA-Interim variables (conversions necessary, see
    # eraint_field.py)
    var_names_era = ['sp', 't2m', 'd2m', 'tcc', 'u10', 'v10', 'tp', 'sf', 'e', 'ssrd', 'strd']
    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Path to ROMS grid
    roms_grid = '/short/m68/kaa561/ROMS-CICE-MCT/apps/common/grid/circ30S_quarterdegree_rp5.nc'

    # Read the ROMS grid and find the minimum and maximum latitudes we care
    # about (ocean points which are not in ice shelf cavities)
    id = Dataset(roms_grid, 'r')
    lat_rho = id.variables['lat_rho'][:,:]
    mask_rho = id.variables['mask_rho'][:,:]
    mask_zice = id.variables['mask_zice'][:,:]
    id.close()
    lat_roms = ma.masked_where(mask_rho-mask_zice==0, lat_rho)
    latS = amin(lat_roms)
    latN = amax(lat_roms)

    # Build a list of Model objects containing 39 CMIP5 models
    models = build_model_list()

    # Set up output arrays of dimension variable x model, with one extra
    # space in the model dimension for the multi-model-mean
    rms_errors = ma.empty([len(var_names_cmip5), len(models)+1])
    relative_errors = ma.empty(shape(rms_errors))

    # Loop over variables
    for i in range(len(var_names_cmip5)):
        var_cmip5 = var_names_cmip5[i]
        print 'Variable ' + var_cmip5
        var_era = var_names_era[i]

        # Read ERA-Interim data
        print 'Processing ERA-Interim'
        era_data, era_lon, era_lat = eraint_field(var_era, start_year, end_year)
        # Trim to the domain we care about (will extend slightly south of latS
        # and slightly north of latN)
        j_min = nonzero(era_lat < latN)[0][0] - 1
        j_max = nonzero(era_lat < latS)[0][0] + 1
        era_lat = era_lat[j_min:j_max]
        era_data = era_data[:,j_min:j_max,:]

        # Make a list of month indices (1 to 12)
        era_months = (arange(size(era_data, 0)+1) % 12 + 1).tolist()
        # Count the number of times each month shows up
        era_num_months = []
        for month in range(1, 12+1):
            era_num_months.append(era_months.count(month))
        # Calculate monthly climatology
        era_climatology = zeros([12, size(era_data,1), size(era_data,2)])
        for t in range(size(era_data,0)):
            curr_month = era_months[t]-1
            # Accumulate averages for each month
            era_climatology[curr_month,:,:] += era_data[t,:,:] / float(era_num_months[curr_month])

        # ERA-Interim has constant spacing for longitude and latitude; find
        # these values (should be 0.75 degrees)
        era_dlon = mean(era_lon[1:] - era_lon[0:-1])
        era_dlat = mean(era_lat[1:] - era_lat[0:-1])
        # Calculate dx and dy in Cartesian space
        # dx = r*cos(lat)*dlon where lat and dlon are converted to radians
        era_dx = abs(r*cos(era_lat*deg2rad)*era_dlon*deg2rad)
        # dy = r*dlat where dlat is converted to radians
        era_dy = abs(r*era_dlat*deg2rad)
        # Define dt (number of seconds in each month)
        era_dt = array([31.0, 28.25, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0, 31.0])*24*60*60
        # Copy dx and dt to 3D arrays of dimension time x latitude x longitude
        # This is a bit complicated
        # We don't have to do this for dy because it's constant
        era_dx = transpose(tile(era_dx, (size(era_climatology,2), 1)))
        era_dx = tile(era_dx, (size(era_climatology,0), 1, 1))
        era_dt = tile(era_dt, (size(era_climatology,2), size(era_climatology,1), 1))
        era_dt = transpose(era_dt)        

        # Initialise the number of models used in the multi-model mean
        num_models = 0
        multi_model_mean = None
        # Loop over models
        for j in range(len(models)):
            model = models[j]
            print 'Processing model ' + model.name

            # Read model data including the grid and month indices
            model_data, model_lon, model_lat, tmp, model_months = cmip5_field(model, expt, var_cmip5, start_year, end_year)
            # Make sure there is actually data (some variables are missing)
            if model_data is not None:

                # Count the number of times each month shows up
                model_num_months = []
                for month in range(1, 12+1):
                    model_num_months.append(model_months.count(month))
                # Calculate monthly climatology  
                model_climatology = zeros([12, size(model_data,1), size(model_data,2)])
                for t in range(size(model_data,0)):
                    curr_month = model_months[t]-1
                    # Accumulate averages for each month
                    model_climatology[curr_month,:,:] += model_data[t,:,:]/float(model_num_months[curr_month])

                # Interpolate to ERA-Interim grid, one month at a time
                model_climatology_interp = ma.empty(shape(era_climatology))
                for t in range(size(model_climatology, 0)):
                    model_climatology_interp[t,:,:] = interp_model2era(model_climatology[t,:,:], model_lon, model_lat, era_lon, era_lat)

                # Add to multi-model mean and increment num_models
                if multi_model_mean is None:
                    multi_model_mean = model_climatology_interp[:,:,:]
                else:
                    multi_model_mean[:,:,:] += model_climatology_interp[:,:,:]
                num_models +=1

                # Root-mean-squared error is the square root of the sum of the
                # squares of the residuals at each point (x, y, t), weighted by
                # dx*dy*dt
                rms_errors[i,j] = sqrt(sum((model_climatology_interp - era_climatology)**2*era_dx*era_dy*era_dt)/sum(era_dx*era_dy*era_dt))
            else:
                rms_errors[i,j] = ma.masked

        # Divide multi_model_mean (currently just a sum) by num_models to get
        # the mean
        multi_model_mean /= num_models
        # Calculate rms error for multi-model mean as well
        rms_errors[i,-1] = sqrt(sum((multi_model_mean - era_climatology)**2*era_dx*era_dy*era_dt)/sum(era_dx*era_dy*era_dt))

        # Find median rms error for this variable across models (not including
        # multi-model mean)
        median_error = median(rms_errors[i,:-1])
        for j in range(len(models)+1):
            # Calculate relative error; for example 0.3 means the given model
            # has a 30% higher rms error than the median
            relative_errors[i,j] = (rms_errors[i,j] - median_error)/rms_errors[i,j]

    # Make a list of model names
    labels = []
    for model in models:
        labels.append(model.name)
    labels.append('MMM')

    # Make portrait plot
    fig = figure(figsize=(8,20))
    ax = fig.add_subplot(111)
    img = ax.pcolormesh(transpose(relative_errors), vmin=-0.8, vmax=0.8, cmap='RdBu_r')

    # Configure plot
    xticks(arange(0.5, len(var_names_cmip5)+0.5), var_names_cmip5, rotation=-45)
    yticks(arange(0.5, len(labels)+0.5), labels)
    xlim([0, len(var_names_cmip5)])
    ylim([len(labels), 0])
    title('Relative error', fontsize=16)
    box = ax.get_position()
    ax.set_position([box.x0+box.width*0.1, box.y0, box.width*0.9, box.height])
    cbaxes = fig.add_axes([0.3,0.03,0.5,0.025])
    colorbar(img, orientation='horizontal', ticks=arange(-0.75, 1, 0.25), cax=cbaxes)
    savefig('relative_errors_eraint.png')

    # Write rms errors to text file, tabulated nicely in a matrix of models x
    # variables
    f = open('rms_errors_eraint.txt', 'w')
    f.write('{0:16}'.format('Model'))
    for var in var_names_cmip5:
        f.write('{0:16}'.format(var))
    f.write('\n')
    for j in range(len(models)):
        f.write('{0:16}'.format(models[j].name))
        for i in range(len(var_names_cmip5)):
            f.write('{0:16}'.format(str(round(rms_errors[i,j],4))))
        f.write('\n')
    f.write('{0:16}'.format('MMM'))
    for i in range(len(var_names_cmip5)):
        f.write('{0:16}'.format(str(round(rms_errors[i,-1],4))))
    f.write('\n')
    f.close()

    # Similarly for relative errors
    f = open('relative_errors_eraint.txt', 'w')
    f.write('{0:16}'.format('Model'))
    for var in var_names_cmip5:
        f.write('{0:16}'.format(var))
    f.write('\n')
    for j in range(len(models)):
        f.write('{0:16}'.format(models[j].name))
        for i in range(len(var_names_cmip5)):
            f.write('{0:16}'.format(str(round(relative_errors[i,j],4))))
        f.write('\n')
    f.write('{0:16}'.format('MMM'))
    for i in range(len(var_names_cmip5)):
        f.write('{0:16}'.format(str(round(relative_errors[i,-1],4))))
    f.write('\n')
    f.close()


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
    model_lon_wrap[0] = model_lon[-1]-360
    model_lon_wrap[1:-1] = model_lon
    model_lon_wrap[-1] = model_lon[0]+360
    # Copy the westernmost and easternmost data points to match
    model_data_wrap = ma.array(zeros((size(model_lat), size(model_lon_wrap))))
    model_data_wrap[:,1:-1] = model_data
    model_data_wrap[:,0] = model_data_wrap[:,-1]
    model_data_wrap[:,-1] = model_data_wrap[:,0]

    # Get 2D mesh of ERA-Interim lat and lon
    era_lon_2d, era_lat_2d = meshgrid(era_lon, era_lat)
    # Build an interpolation function for model_data
    interp_function = RegularGridInterpolator((model_lat, model_lon_wrap), model_data_wrap)
    # Call it for the ERA-Interim grid
    data_interp = interp_function((era_lat_2d, era_lon_2d))

    return data_interp


# Command-line interface
if __name__ == '__main__':

    cmip5_eraint_rms_errors()
