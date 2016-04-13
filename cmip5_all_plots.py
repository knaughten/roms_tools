from cmip5_paths import *
from cmip5_plot import *

# Call cmip5_plot for all variables, all seasons, and all models.
# Input:
# season = string specifying which season to average over ('djf', 'mam', 'jja',
#          'son') or whether to average over all seasons ('annual')
def cmip5_all_plots (season):

    # Variable names for CMIP5
    var_names_cmip5 = ['ps', 'tas', 'huss', 'clt', 'uas', 'vas', 'pr', 'prsn', 'evspsbl', 'rsds', 'rlds', 'thetao', 'so', 'uo', 'vo']

    # Build an array of Model objects, one for each of 39 CMIP5 models
    models = build_model_list()

    # Loop over variables
    for i in range(len(var_names_cmip5)):
        var = var_names_cmip5[i]
        print 'Plotting ' + var
        cmip5_plot(var, season, models, True, 'cmip5_figures/' + var + '_' + season + '.png')


# Command-line interface
if __name__ == "__main__":

    # Run through all seasons
    #cmip5_all_plots('djf')
    #cmip5_all_plots('mam')
    #cmip5_all_plots('jja')
    cmip5_all_plots('son')
    cmip5_all_plots('annual')

    
                     
    
