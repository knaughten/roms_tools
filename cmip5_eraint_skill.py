from cmip5_paths import *
from cmip5_field_old import *
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
# Input:
# season = string specifying which season to average over ('djf', 'mam', 'jja',
#          'son') or whether to average over all seasons ('annual')
def cmip5_eraint_skill (season):

    # Set parameters
    # Experiment name
    expt = 'historical'
    # Years to average over
    start_year = 1992
    end_year = 2005
    # Path to ROMS grid file
    roms_grid = '/short/m68/kaa561/ROMS-CICE-MCT/apps/common/grid/circ30S_quarterdegree_rp5.nc'

    # Variable names for CMIP5
    var_names_cmip5 = ['ps', 'tas', 'huss', 'clt', 'uas', 'vas', 'pr', 'prsn', 'evspsbl', 'rsds', 'rlds']
    # Corresponding variable names for ERA-Interim
    var_names_era = ['sp', 't2m', 'd2m', 'tcc', 'u10', 'v10', 'tp', 'sf', 'e', 'ssrd', 'strd']

    # Read ROMS latitude and mask out land and ice shelf points
    id = Dataset(roms_grid, 'r')
    lat_rho = id.variables['lat_rho'][:,:]
    mask_rho = id.variables['mask_rho'][:,:]
    mask_zice = id.variables['mask_zice'][:,:]
    id.close()
    lat_roms = ma.masked_where(mask_rho-mask_zice==0, lat_rho)
    # Find latitude bounds of ROMS open ocean points
    latS = amin(lat_roms)
    latN = amax(lat_roms)

    # Work out which months we want to average over
    if season == 'djf':
        months = [12, 1, 2]
    elif season == 'mam':
        months = [3, 4, 5]
    elif season == 'jja':
        months = [6, 7, 8]
    elif season == 'son':
        months = [9, 10, 11]
    elif season == 'annual':
        months = range(1,12+1)
    print 'Time-averaging over season ' + season.upper()

    # Build an array of Model objects, one for each of 39 CMIP5 models
    models = build_model_list()

    # Loop over variables
    for i in range(len(var_names_cmip5)):

        # Select variable names
        var_cmip5 = var_names_cmip5[i]
        var_era = var_names_era[i]
        print 'Variable ' + var_cmip5

        print 'Processing ERA-Interim'
        era_data, era_lon, era_lat = eraint_field(var_era, start_year, end_year)
        # Zonally average - this is easy on a regular grid
        era_data_zonalavg = mean(era_data, axis=2)
        # Mask out months that we don't care about
        for t in range(size(era_data_zonalavg, 0)):
            curr_month = (t+1) % 12
            if curr_month not in months:
                era_data_zonalavg[t,:] = ma.masked
        # Time average (this will only capture the correct months)
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
            model_data, model_lat, model_months = cmip5_field_old(model, expt, var_cmip5, start_year, end_year)

            if model_data is not None:
                # Zonally average - note all CMIP5 models have regular grids
                model_data_zonalavg = mean(model_data, axis=2)
                # Mask out months that we don't care about
                for t in range(size(model_data_zonalavg, 0)):
                    if model_months[t] not in months:
                        model_data_zonalavg[t,:] = ma.masked
                # Time average (this will only capture the correct months)
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
        error_file = 'cmip5_skill_scores/' + var_cmip5 + '_' + season + '.txt'
        print 'Writing ' + error_file
        f = open(error_file, 'w')
        for i in sort_index:
            f.write(labels[i] + '   ' + str(error_scores[i]) + '\n')
        f.close()


# Command-line interface
if __name__ == "__main__":

    # Run through all seasons
    cmip5_eraint_skill('djf')
    cmip5_eraint_skill('mam')
    cmip5_eraint_skill('jja')
    cmip5_eraint_skill('son')
    cmip5_eraint_skill('annual')
    
