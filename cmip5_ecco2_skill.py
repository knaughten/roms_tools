from cmip5_paths import *
from cmip5_field import *
from ecco2_field import *
from numpy import *
from netCDF4 import Dataset
from scipy.interpolate import interp1d

# Calculate skill scores (sum of squares of residuals, divided by number of
# points) for each CMIP5 model compared to ECCO2 reanalyses, for each of 4
# ocean variables, zonally averaged over the northern boundary of the ROMS
# circumpolar grid (currently 30S) and time-averaged between 1995 and 2005.
# For each variable, rank the models by their skill scores and output the
# results in a plain text file.
# Input:
# season = string specifying which season to average over ('djf', 'mam', 'jja',
#          'son') or whether to average over all seasons ('annual')
def cmip5_ecco2_skill (season):

    # Set parameters
    # Experiment name
    expt = 'historical'
    # Years to average over
    start_year = 1992
    end_year = 2005

    # Variable names for CMIP5
    var_names_cmip5 = ['thetao', 'so', 'uo', 'vo']
    # Corresponding variable names for ECCO2
    var_names_ecco2 = ['THETA', 'SALT', 'UVEL', 'VVEL']    

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
        var_ecco2 = var_names_ecco2[i]
        print 'Variable ' + var_cmip5

        print 'Processing ECCO2'
        ecco2_data, ecco2_depth = ecco2_field(var_ecco2, start_year, end_year)
        # Zonally average - this is easy on a regular grid
        ecco2_data_zonalavg = mean(ecco2_data, axis=2)
        # Mask out months that we don't care about
        for t in range(size(ecco2_data_zonalavg, 0)):
            curr_month = (t+1) % 12
            if curr_month not in months:
                ecco2_data_zonalavg[t,:] = ma.masked
        # Time average (this will only capture the correct months)
        ecco2_data_timeavg = mean(ecco2_data_zonalavg, axis=0)

        # Loop through Model objects
        labels = []
        error_scores = []
        for model in models:

            print 'Processing ' + model.name
            # Get the model output for this variable, and the model's depth
            # axis (they are all on different grids)
            model_data, model_depth, model_months = cmip5_field(model, expt, var_cmip5, start_year, end_year)

            if model_data is not None:
                # Zonally average - note all CMIP5 models have regular grids
                model_data_zonalavg = mean(model_data, axis=2)
                # Mask out months that we don't care about
                for t in range(size(model_data_zonalavg, 0)):
                    if model_months[t] not in months:
                        model_data_zonalavg[t,:] = ma.masked
                # Time average (this will only capture the correct months)
                model_data_timeavg = mean(model_data_zonalavg, axis=0)

                if model_depth[1] < model_depth[0]:
                    # Model depth array starts at seafloor; flip it and the data
                    model_depth = flipud(model_depth)
                    model_data_timeavg = flipud(model_data_timeavg)

                # Select the ECCO2 depth points which fall within this model's
                # depth range
                index = nonzero((ecco2_depth>=amin(model_depth))*(ecco2_depth<=amax(model_depth)))                

                # Interpolate the resulting 1D array to the ECCO2 depth points
                interp_function = interp1d(model_depth, model_data_timeavg)
                model_data_interp = interp_function(ecco2_depth[index])

                # Sum squares of residuals and divide by the number of depth
                # points to get error estimate
                error = sum((model_data_interp - ecco2_data_timeavg[index])**2)/size(ecco2_depth[index])
                # Append to list
                labels.append(model.name)
                error_scores.append(error)

        # Get indices of sorted error estimates
        sort_index = argsort(error_scores)

        # Write model names and their error estimates to a file, in increasing
        # order of error estimate (eg models most similar to ECCO2 on top)
        error_file = 'cmip5_skill_scores/' + var_cmip5 + '_' + season + '.txt'
        print 'Writing ' + error_file
        f = open(error_file, 'w')
        for i in sort_index:
            f.write(labels[i] + '   ' + str(error_scores[i]) + '\n')
        f.close()


# Command-line interface
if __name__ == "__main__":

    # Run through all seasons
    #cmip5_ecco2_skill('djf')
    #cmip5_ecco2_skill('mam')
    cmip5_ecco2_skill('jja')
    cmip5_ecco2_skill('son')
    cmip5_ecco2_skill('annual')
