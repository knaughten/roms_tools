from cmip5_paths import *
from cmip5_field import *
from eraint_field import *
from numpy import *
from netCDF4 import Dataset
from scipy.interpolate import interp1d

# Calculate skill scores (sum of squares of residuals, divided by number of
# points) for each CMIP5 model compared to ERA-Interim reanalyses, for each
# of 11 atmospheric variables, zonally averaged over the Southern Ocean (with
# latitude bounds corresponding to the ROMS circumpolar grid) and time-averaged
# between 1995 and 2005. For each variable, rank the models by their skill 
# scores and output the results in a plain text file.
def cmip5_eraint_skill ():

    # Set parameters
    # Experiment name
    expt = 'historical'
    # Years to average over
    start_year = 1995
    end_year = 2005
    # Path to ROMS grid file
    roms_grid = '/short/m68/kaa561/ROMS-CICE-MCT/apps/common/grid/circ38S_quarterdegree.nc'

    # Variable names for CMIP5
    var_names_cmip5 = ['ps', 'tas', 'huss', 'clt', 'uas', 'vas', 'pr', 'prsn', 'evspsbl', 'rsds', 'rlds']
    # Corresponding variable names for ERA-Interim
    var_names_era = ['sp', 't2m', 'd2m', 'tcc', 'u10', 'v10', 'tp', 'sf', 'e', 'ssrd', 'strd']

    # Read ROMS latitude and apply land mask
    id = Dataset(roms_grid, 'r')
    lat_rho = id.variables['lat_rho'][:,:]
    mask_rho = id.variables['mask_rho'][:,:]
    id.close()
    lat_roms = ma.masked_where(mask_rho==0, lat_rho)
    # Find latitude bounds of ROMS ocean points
    latS = amin(lat_roms)
    latN = amax(lat_roms)

    # Build an array of Model objects, one for each of 39 CMIP5 models
    models = build_model_list()

    # Loop over variables
    for i in range(len(var_names_cmip5)):

        # Select variable names
        var_cmip5 = var_names_cmip5[i]
        var_era = var_names_era[i]
        print 'Variable ' + var_cmip5

        print 'Processing ERA-Interim'
        era_data, era_lat = eraint_field(var_era, start_year, end_year)
        # Zonally average - this is easy on a regular grid
        era_data_zonalavg = mean(era_data, axis=2)
        # Time avearge - also easy because equally spaced time indices
        era_data_timeavg = mean(era_data_zonalavg, axis=0)

        # Trim latitude bounds so we are entirely within the ROMS domain
        j_min = nonzero(era_lat < latN)[0][0]
        j_max = nonzero(era_lat < latS)[0][0]
        era_lat = era_lat[j_min:j_max]
        era_data_timeavg = era_data_timeavg[j_min:j_max]

        # Loop through Model objects
        labels = []
        error_scores = []
        for model in models:

            print 'Processing ' + model.name
            # Get the model output for this variable, and the model's latitude
            # axis (they are all on different grids)
            model_data, model_lat = cmip5_field(model, expt, var_cmip5, start_year, end_year)

            if model_data is not None:
                # Zonally average - note all CMIP5 models have regular grids
                model_data_zonalavg = mean(model_data, axis=2)
                # Time average
                model_data_timeavg = mean(model_data_zonalavg, axis=0)

                # Interpolate the resulting 1D array to the ERA-Interim
                # latitude points
                interp_function = interp1d(model_lat, model_data_timeavg)
                model_data_interp = interp_function(era_lat)

                # Sum squares of residuals and divide by number of latitude
                # points to get error estimate
                error = sum((model_data_interp - era_data_timeavg)**2)/size(era_lat)
                # Append to list
                labels.append(model.name)
                error_scores.append(error)

        # Get indices of sorted error estimates
        sort_index = argsort(error_scores)

        # Write model names and their error estimates to a file, in increasing
        # order of error estimate (eg models most similar to ERA-Interim on top)
        error_file = var_cmip5 + '_errors.txt'
        print 'Writing ' + error_file
        f = open(error_file, 'w')
        for i in sort_index:
            f.write(labels[i] + '   ' + str(error_scores[i]) + '\n')
        f.close()


# Command-line interface
if __name__ == "__main__":

    cmip5_eraint_skill()
