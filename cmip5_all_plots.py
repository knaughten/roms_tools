from cmip5_paths import *
from cmip5_plot import *

# Call cmip5_plot.py for all models (including the multi-model mean), all
# variables (atmosphere and ocean), and all seasons.
# Note that in order to run this script, you must first run:
# eraint_climatology_netcdf.py, ecco2_climatology_netcdf.py,
# cmip5_atmos_climatology_netcdf.py, cmip5_ocean_climatology_netcdf.py,
# mmm_atmos_netcdf.py, and mmm_ocean_netcdf.py
# to generate the necessary NetCDF files. Also, the directory "cmip5/" must
# exist.
# Input: season = string containing the season key: 'djf', 'mam', 'jja', 'son',
#                 or 'annual'
def cmip5_all_plots (season):

    # All possible variable names
    var_names = ['Pair', 'Tair', 'Hair', 'cloud', 'Uwind', 'Vwind', 'precip', 'snow', 'evap', 'swrad', 'lwrad', 'temp', 'salt', 'u', 'v']

    # Make a list of all possible model names
    models = build_model_list()
    model_names = []
    for model in models:
        model_names.append(model.name)
    model_names.append('MMM')

    # Call cmip5_plot for each variable
    for i in range(len(var_names)):
        var = var_names[i]
        print 'Plotting ' + var
        cmip5_plot(var, season, model_names, True, 'cmip5/' + var + '_' + season + '.png')


# Command-line interface
if __name__ == "__main__":

    # Run for each season key at a time
    cmip5_all_plots('djf')
    cmip5_all_plots('mam')
    cmip5_all_plots('jja')
    cmip5_all_plots('son')
    cmip5_all_plots('annual')

    
                     
    
