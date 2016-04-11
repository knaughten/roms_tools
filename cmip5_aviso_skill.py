from cmip5_paths import *
from cmip5_field_ss import *
from aviso_field import *
from numpy import *
from netCDF4 import Dataset

# Calculate skill scores (square of residual) for each CMIP5 model compared to
# AVISO reanalyses of sea surface height, zonally averaged over the northern
# boundary of the ROMS circumpolar grid (currently 30S) and time-averaged
# between 1995 and 2005. Rank the models by their skill scores and output the
# results in a plain text file.
def cmip5_aviso_skill ():

    # Set parameters
    # Experiment name
    expt = 'historical'
    # Years to average over
    start_year = 1995
    end_year = 2005    

    # Build an array of Model objects, one for each of 39 CMIP5 models
    models = build_model_list()

    print 'Processing AVISO data'
    aviso_data = aviso_field(start_year, end_year)
    # Zonally average - this is easy on a regular grid
    aviso_data_zonalavg = mean(aviso_data, axis=1)
    # Time average - also easy because equally spaced time indices
    aviso_data_timeavg = mean(aviso_data_zonalavg, axis=0)

    # Loop through Model objects
    labels = []
    error_scores = []
    for model in models:
        print 'Processing ' + model.name
        # Get the model output for this variable
        model_data = cmip5_field_ss(model, expt, 'zos', start_year, end_year)

        # Check if the output actually exists (if not, cmip5_field will
        # return None)
        if model_data is not None:
            # Zonally average - note all CMIP5 models have regular grids
            model_data_zonalavg = mean(model_data, axis=1)
            # Time average
            model_data_timeavg = mean(model_data_zonalavg, axis=0)

            # Calculate square of residual
            error = (model_data_timeavg - aviso_data_timeavg)**2
            labels.append(model.name)
            error_scores.append(error)

    # Get indices of sorted error estimates
    sort_index = argsort(error_scores)

    # Write model names and their error estimates to a file, in increasing
    # order of error estimate (eg models most similar to AVISO on top)
    error_file = 'zos_errors.txt'
    print 'Writing ' + error_file
    f = open(error_file, 'w')
    for i in sort_index:
        f.write(labels[i] + '   ' + str(error_scores[i]) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    cmip5_aviso_skill()
