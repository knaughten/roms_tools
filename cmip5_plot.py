from cmip5_paths import *
from cmip5_field import *
from eraint_field import *
from ecco2_field import *
from matplotlib.pyplot import *
from numpy import *
from netCDF4 import Dataset
from math import ceil
from matplotlib.font_manager import FontProperties

# Compare CMIP5 output from the given models to reanalyses (ERA-Interim for 
# atmosphere variables, ECCO2 for ocean variables) by creating a plot of the
# given variable, zonally averaged over the ROMS domain and time-averaged
# between 1992 and 2005, with the given variable on the x-axis and latitude
# (for atmospheric surface variables) or depth (for ocean northern boundary
# variables) on the y-axis.
# Input:
# var = string containing name of CMIP5 variable (ps, tas, huss, clt, uas, vas,
#       pr, prsn, evspsbl, rsds, rlds, thetao, so, uo, vo)
# season = string specifying which season to average over ('djf', 'mam', 'jja',
#          'son') or whether to average over all seasons ('annual')
# models = array of Model objects (see cmip5_paths.py)
# save = optional boolean flag indicating that the plot should be saved to a
#        file rather than displayed on the screen
# fig_name = if save=True, filename for figure
def cmip5_plot (var, season, models, save=False, fig_name=None):

    # Set parameters
    # Experiment name
    expt = 'historical'
    # Years to average over
    start_year = 1992
    end_year = 2005
    # Latitude of the northern boundary of the circmpolar ROMS domain
    nbdry = -30

    # Figure out whether it is an atmosphere or ocean variable
    if var in ['ps', 'tas', 'huss', 'clt', 'uas', 'vas', 'pr', 'prsn', 'evspsbl', 'rsds', 'rlds']:
        realm = 'atmos'
    elif var in ['thetao', 'so', 'uo', 'vo']:
        realm = 'ocean'
    else:
        print 'Unknown variable'
        # Exit early
        return None

    # Set plot title, plot units, and corresponding ERA-Interim or ECCO2
    # variable name
    if var == 'ps':
        var_era = 'sp'
        plot_title = 'Surface Pressure'
        plot_units = 'kPa'
    elif var == 'tas':
        var_era = 't2m'
        plot_title = 'Surface Air Temperature'
        plot_units = r'$^{\circ}$C'
    elif var == 'huss':
        var_era = 'd2m'
        plot_title = 'Surface Specific Humidity'
        plot_units = '1'
    elif var == 'clt':
        var_era = 'tcc'
        plot_title = 'Total Cloud Cover'
        plot_units = '%'
    elif var == 'uas':
        var_era = 'u10'
        plot_title = 'Eastward Wind Velocity'
        plot_units = 'm/s'
    elif var == 'vas':
        var_era = 'v10'
        plot_title = 'Northward Wind Velocity'
        plot_units = 'm/s'
    elif var == 'pr':
        var_era = 'tp'
        plot_title = 'Total Precipitation'
        plot_units = r'10$^6$ kg m$^{-2}$ s$^{-1}$'
    elif var == 'prsn':
        var_era = 'sf'
        plot_title = 'Snowfall'
        plot_units = r'10$^6$ kg m$^{-2}$ s$^{-1}$'
    elif var == 'evspsbl':
        var_era = 'e'
        plot_title = 'Evaporation'
        plot_units = r'10$^6$ kg m$^{-2}$ s$^{-1}$'
    elif var == 'rsds':
        var_era = 'ssrd'
        plot_title = 'Downgoing Shortwave Radiation'
        plot_units = r'W m$^{-2}$'
    elif var == 'rlds':
        var_era = 'strd'
        plot_title = 'Downgoing Longwave Radiation'
        plot_units = r'W m$^{-2}$'
    elif var == 'thetao':
        var_ecco2 = 'THETA'
        plot_title = 'Temperature'
        plot_units = r'$^{\circ}$C'
    elif var == 'so':
        var_ecco2 = 'SALT'
        plot_title = 'Salinity'
        plot_units = 'psu'
    elif var == 'uo':
        var_ecco2 = 'UVEL'
        plot_title = 'Eastward Velocity'
        plot_units = 'm/s'
    elif var == 'vo':
        var_ecco2 = 'VVEL'
        plot_title = 'Northward Velocity'
        plot_units = 'm/s'

    # RGB colours to be used for the max. 39 CMIP5 models
    # These 39 visually distinct (or as visually distinct as possible) colours
    # were generated from http://phrogz.net/css/distinct-colors.html
    all_cmip5_colours = [(0.83,0,0), (0.49,0,0), (0.32,0,0), (1,0.39,0.39), (0.66,0.26,0.26), (0.83,0.39,0), (0.49,0.23,0), (0.83,0.56,0.33), (0.49,0.33,0.19), (1,0.93,0), (0.66,0.61,0), (0.49,0.47,0.19), (0.19,0.32,0), (0.76,1,0.39), (0.31,0.66,0.26), (0,1,0.33), (0.13,0.32,0.19), (0,0.49,0.39), (0.39,1,0.88), (0,0.24,0.32), (0.33,0.69,0.83), (0.19,0.41,0.49), (0,0.18,0.66), (0.33,0.46,0.83), (0.19,0.27,0.49), (0.2,0.2,0.2), (0.16,0,0.83), (0.1,0,0.49), (0.51,0.39,1), (0.67,0,1), (0.21,0,0.32), (0.39,0.19,0.49), (1,0.39,0.92), (0.66,0.26,0.61), (0.32,0.13,0.29), (0.83,0,0.33), (0.83,0.33,0.53), (0.49,0.19,0.31), (1,0,0)]
    num_colours = len(all_cmip5_colours)
    num_models = len(models)
    cmip5_colours = []
    # Select num_models number of colours evenly from all_cmip5_colours
    for i in range(num_models):
        index = int(ceil(i*float(num_colours)/num_models))
        cmip5_colours.append(all_cmip5_colours[index])

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

    # Set up figure window
    figure(figsize=(12,9))
    ax = subplot(111)

    if realm == 'atmos':

        print 'Processing ERA-Interim'
        era_data, era_lat = eraint_field(var_era, start_year, end_year)
        # Zonally average - this is easy on a regular grid
        era_data_zonalavg = mean(era_data, axis=2)
        # Mask out months that we don't care about
        for t in range(size(era_data_zonalavg, 0)):
            curr_month = (t+1) % 12
            if curr_month not in months:
                era_data_zonalavg[t,:] = ma.masked
        # Time average (this will only capture the correct months)
        era_data_timeavg = mean(era_data_zonalavg, axis=0)
        # Plot ERA-Interim data
        ax.plot(era_data_timeavg, era_lat, label='ERA-Interim', color='black', linewidth=3, zorder=num_models+1)

    elif realm == 'ocean':

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
        # Plot ECCO2 data
        ax.plot(ecco2_data_timeavg, ecco2_depth, label='ECCO2', color='black', linewidth=3, zorder=num_models+1)

    # Loop through Model objects with a manual counter j
    j = 0
    for model in models:
        print 'Processing ' + model.name
        # Get the model output for this variable, the model's latitude or depth
        # axis (they are all on different grids), and the month number of each
        # time index
        model_data, model_axis, model_months = cmip5_field(model, expt, var, start_year, end_year)
        # Check for missing data
        if model_data is not None:
            # Zonally average - note all CMIP5 models have regular grids
            model_data_zonalavg = mean(model_data, axis=2)
            # Mask out months that we don't care about
            for t in range(size(model_data_zonalavg, 0)):
                if model_months[t] not in months:
                    model_data_zonalavg[t,:] = ma.masked
            # Time average (this will only capture the correct months)
            model_data_timeavg = mean(model_data_zonalavg, axis=0)
            # Add this model's output to the plot, selecting the correct colour
            ax.plot(model_data_timeavg, model_axis, label=model.name, color=cmip5_colours[j], linewidth=2)
        j += 1

    # Configure plot
    title(plot_title)
    xlabel(plot_units)
    grid(True)
    if realm == 'atmos':
        ylim([-90, nbdry])
        ylabel('Latitude')
    elif realm == 'ocean':
        ylabel('Depth (m)')
        gca().invert_yaxis()
    # Move the plot over a bit so there's room for the legend beside it
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    fontP = FontProperties()
    fontP.set_size('small')
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5), prop=fontP)

    if save:
        savefig(fig_name)
    else:
        show()    


# Command-line interface
if __name__ == "__main__":

    # Make a list of all possible Model objects
    all_models = build_model_list()
    # Make a list of corresponding model names
    all_model_names = []
    for model in all_models:
        all_model_names.append(model.name)
    models = []

    # Get variable name and season code
    var = raw_input("CMIP5 variable name: ")
    season = raw_input("Season (djf, mam, jja, son, or annual): ")

    # Get model name
    model_name = raw_input("Model name: ")
    try:
        # Try to find a match in master list of Models
        model_index = all_model_names.index(model_name)
        models.append(all_models[model_index])
    except ValueError:
        print 'No such model'
    # Keep asking for more models until the user is done
    while True:
        action = raw_input("Add another model (y/n)? ")
        if action == 'y':
            model_name = raw_input("Model name: ")
            try:
                model_index = all_model_names.index(model_name)
                models.append(all_models[model_index])
            except ValueError:
                print 'No such model'
        else:
            break

    # Get save/display choice
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        # Get file name for figure
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None

    # Make the plot
    cmip5_plot(var, season, models, save, fig_name)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the input parameters; repeat until the user is finished
                changes = raw_input("Enter a parameter to change: (1) CMIP5 variable, (2) season, (3) models, (4) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New variable name
                        var = raw_input("CMIP5 variable name: ")
                    elif int(changes) == 2:
                        # New season
                        season = raw_input("Season (djf, mam, jja, son, or annual): ")
                    elif int(changes) == 3:
                        # New list of models
                        models = []
                        model_name = raw_input("Model name: ")
                        try:
                            model_index = all_model_names.index(model_name)
                            models.append(all_models[model_index])
                        except ValueError:
                            print 'No such model'
                        while True:
                            action = raw_input("Add another model (y/n)? ")
                            if action == 'y':
                                model_name = raw_input("Model name: ")
                                try:
                                    model_index = all_model_names.index(model_name)
                                    models.append(all_models[model_index])
                                except ValueError:
                                    print 'No such model'
                            else:
                                break
                    elif int(changes) == 4:
                        # Change from display to save, or vice versa 
                        save = not save
            if save:
                # Get file name for figure
                fig_name = raw_input("File name for figure: ")
            # Make the plot
            cmip5_plot(var, season, models, save, fig_name)
        else:
            break
    
                        
