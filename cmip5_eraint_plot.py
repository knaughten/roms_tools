from cmip5_paths import *
from cmip5_field import *
from matplotlib.pyplot import *
from numpy import *
from netCDF4 import Dataset
from matplotlib.font_manager import FontProperties

# Compare ERA-Interim reanalyses and CMIP5 output by creating plots for 11
# atmospheric surface variables, zonally averaged over the ROMS domain and 
# time-averaged between 1995 and 2005, with the given variable on the 
# x-axis and latitude on the y-axis.
def cmip5_eraint_plot ():

    # Set parameters
    # Experiment name
    expt = 'historical'
    # Years to average over
    start_year = 1995
    end_year = 2005
    # Latitude of the northern boundary of the circumpolar ROMS domain
    nbdry = -38
    # Directory where ERA-Interim monthly averaged data is stored
    era_dir = '/short/y99/kaa561/FESOM/ERA_Interim_monthly/'
    # String that ERA-Interim files end with
    era_tail = '_monthly_orig.nc'

    # Variable names for CMIP5
    var_names_cmip5 = ['ps', 'tas', 'huss', 'clt', 'uas', 'vas', 'pr', 'prsn', 'evspsbl', 'rsds', 'rlds']
    # Corresponding variable names for ERA-Interim (with some conversions
    # that will be applied below)
    var_names_era = ['sp', 't2m', 'd2m', 'tcc', 'u10', 'v10', 'tp', 'sf', 'e', 'ssrd', 'strd']
    # Long names for variables to be used as plot titles
    plot_titles = ['Surface Pressure', 'Surface Air Temperature', 'Surface Specific Humidity', 'Total Cloud Cover', 'Eastward Wind Velocity', 'Northward Wind Velocity', 'Total Precipitation', 'Snowfall', 'Evaporation', 'Downgoing Shortwave Radiation', 'Downgoing Longwave Radiation']
    # Units for each variable
    plot_units = ['kPa', r'$^{\circ}$C', '1', '%', 'm/s', 'm/s', r'10$^6$ kg m$^{-2}$ s$^{-1}$', r'10$^6$ kg m$^{-2}$ s$^{-1}$', r'10$^6$ kg m$^{-2}$ s$^{-1}$', r'W m$^{-2}$', r'W m$^{-2}$']
    # Latent heat of vapourisation, J/kg
    Lv = 2.5e6
    # Ideal gas constant for water vapour, J/K/kg
    Rv = 461.5
    # Density of water, kg/m^3
    rho_w = 1e3
    # Conversion from K to C
    degKtoC = -273.15

    # RGB colours to be used for the 39 CMIP5 models
    # These 39 visually distinct (or as visually distinct as possible) colours
    # were generated from http://phrogz.net/css/distinct-colors.html
    cmip5_colours = [(0.83,0,0), (0.49,0,0), (0.32,0,0), (1,0.39,0.39), (0.66,0.26,0.26), (0.83,0.39,0), (0.49,0.23,0), (0.83,0.56,0.33), (0.49,0.33,0.19), (1,0.93,0), (0.66,0.61,0), (0.49,0.47,0.19), (0.19,0.32,0), (0.76,1,0.39), (0.31,0.66,0.26), (0,1,0.33), (0.13,0.32,0.19), (0,0.49,0.39), (0.39,1,0.88), (0,0.24,0.32), (0.33,0.69,0.83), (0.19,0.41,0.49), (0,0.18,0.66), (0.33,0.46,0.83), (0.19,0.27,0.49), (0.2,0.2,0.2), (0.16,0,0.83), (0.1,0,0.49), (0.51,0.39,1), (0.67,0,1), (0.21,0,0.32), (0.39,0.19,0.49), (1,0.39,0.92), (0.66,0.26,0.61), (0.32,0.13,0.29), (0.83,0,0.33), (0.83,0.33,0.53), (0.49,0.19,0.31), (1,0,0)]
    # FontProperties object to have smaller font size for legends
    fontP = FontProperties()
    fontP.set_size('small')

    # Build an array of Model objects, one for each of 39 CMIP5 models
    models = build_model_list()

    # Read ERA-Interim latitude and longitude from the first file
    id = Dataset(era_dir + 'AN_' + str(start_year) + '_monthly_orig.nc', 'r')
    era_lat = id.variables['lat'][:]
    era_lon = id.variables['lon'][:]
    id.close()
    # Latitude values are decreasing for some reason, i.e. 90 to -90
    # Find the first index south of nbdry, and subtract 1 to find the last
    # index north of nbdry
    j_min = nonzero(era_lat < nbdry)[0][0] - 1
    # Save only the values from j_min southward
    era_lat = era_lat[j_min:]

    # Loop over variables
    for i in range(len(var_names_cmip5)):

        print 'Plotting ' + plot_titles[i]
        # Select variable names
        var_cmip5 = var_names_cmip5[i]
        var_era = var_names_era[i]

        # Figure out how ERA-Interim filename will start; the atmospheric
        # variables for each year are split between 3 different files
        if var_era in ['sp', 't2m', 'd2m', 'tcc', 'u10', 'v10']:
            era_head = era_dir + 'AN_'
        elif var_era in ['tp', 'sf']:
            era_head = era_dir + 'FC_'
        elif var_era in ['e', 'ssrd', 'strd']:
            era_head = era_dir + 'ER_'

        print 'Processing ERA-Interim'

        # Create empty array of dimension time x latitude x longitude
        era_data = ma.empty([12*(end_year-start_year+1), size(era_lat), size(era_lon)])
        # Initialise next available time index in this array
        posn = 0
        # Loop over years
        for year in range(start_year, end_year+1):

            # Construct filename
            era_file = era_head + str(year) + era_tail
            # Read data
            id = Dataset(era_file, 'r')
            data = id.variables[var_era][:,j_min:,:]

            # Perform conversions if necessary
            if var_cmip5 == 'ps':
                # Convert pressure from Pa to kPa
                data = 1e-3*data
            elif var_cmip5 == 'tas':
                # Convert temperature from K to C
                data = data + degKtoC
            elif var_cmip5 == 'huss':
                # Calculate specific humidity from dewpoint temperature
                # and surface pressure
                d2m = data
                sp = id.variables['sp'][:,j_min:,:]
                # Intermediate step to calculate vapour pressure
                e = 611*exp(Lv/Rv*(1/273.0 - 1/d2m))
                data = 0.622*e/(sp - 0.378*e)
            elif var_cmip5 == 'clt':
                # Convert total cloud cover from fraction to percent
                data = data*100
            elif var_cmip5 in ['pr', 'prsn', 'evspsbl']:
                # Convert precip/snowfall/evap from kg/m^2/s to 1e-6 kg/m^2/s
                data = 1e6*data*rho_w/(12*60*60)
                if var_cmip5 == 'evspsbl':
                    # ERA-Interim defines evaporation as negative; fix this
                    data = -1*data
            elif var_cmip5 in ['rsds', 'rlds']:
                # Convert from J/m^2, integrated over 12 hours, to W/m^2
                data = data/(12*60*60)

            # Save each timestep to master era_data array
            for t in range(size(data, 0)):
                era_data[posn,:,:] = data[t,:,:]
                posn += 1

        # Zonally average - this is easy on a regular grid
        era_data_zonalavg = mean(era_data, axis=2)
        # Time average - also easy because equally spaced time indices
        era_data_timeavg = mean(era_data_zonalavg, axis=0)

        # Plot ERA-Interim data
        figure(figsize=(12,9))
        ax = subplot(111)
        # Keyword "zorder" makes sure this line will be on top of all the
        # thinner CMIP5 lines
        ax.plot(era_data_timeavg, era_lat, label='ERA-Interim', color='black', linewidth=3, zorder=len(models)+1)

        # Loop through Model objects with a manual counter j
        j = 0
        for model in models:

            print 'Processing ' + model.name

            # Get the model output for this variable, and the model's latitude
            # axis (they are all on different grids)
            model_data, model_lat = cmip5_field(model, expt, var_cmip5, start_year, end_year)

            # Check if the output actually exists (if not, cmip5_field will
            # return None)
            if model_data is not None:
                # Zonally average - note all CMIP5 models have regular grids
                model_data_zonalavg = mean(model_data, axis=2)
                # Time average
                model_data_timeavg = mean(model_data_zonalavg, axis=0)

                # Conversions if necessary
                if var_cmip5 in ['pr', 'prsn', 'evspsbl']:
                    # Convert precip/snowfall/evap from kg/m^2/s to
                    # 1e-6 kg/m^2/s
                    model_data_timeavg = 1e6*model_data_timeavg
                elif var_cmip5 == 'ps':
                    # Convert surface pressure from Pa to kPa
                    model_data_timeavg = 1e-3*model_data_timeavg
                elif var_cmip5 == 'tas':
                    # Convert temperature from K to C
                    model_data_timeavg = model_data_timeavg + degKtoC

                # Add this model's output to the plot, selecting the correct
                # colour
                if model.name == 'ACCESS1-0':
                    ax.plot(model_data_timeavg, model_lat, label=model.name, color=cmip5_colours[j], linewidth=3)
                else:
                    ax.plot(model_data_timeavg, model_lat, label=model.name, color=cmip5_colours[j])
            # Increment j whether or not this model had any output, so that
            # colour choices for each model will be consistent between variables
            j += 1

        # Configure plot
        title(plot_titles[i])
        grid(True)
        ylim([-90, nbdry])
        xlabel(plot_units[i])
        ylabel('Latitude')

        # Move the plot over a bit so there's room for the legend beside it
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        # Set font size to small for legend using fontP (previously created)
        # so the massive legend will fit
        ax.legend(loc='center left', bbox_to_anchor=(1,0.5), prop=fontP)
        savefig(var_cmip5 + '.png')


# Command-line interface
if __name__ == "__main__":

    cmip5_eraint_plot()
    
                     
    
