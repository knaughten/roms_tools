from cmip5_paths import *
from cmip5_field import *
from matplotlib.pyplot import *
from numpy import *
from netCDF4 import Dataset
from matplotlib.font_manager import FontProperties

# Compare ECCO2 reanalyses and CMIP5 output by creating plots for 4 ocean
# variables at the northern boundary of the ROMS grid, zonally averaged over
# the ROMS domain and time-averaged between 1995 and 2005, with the given
# variable on the x-axis and depth on the y-axis.
def cmip5_ecco2_plot ():

    # Set parameters
    # Experiment name
    expt = 'historical'
    # Years to average over
    start_year = 1995
    end_year = 2005
    # Latitude of the northern boundary of the circumpolar ROMS domain
    nbdry = -38
    # Directory where ECCO2 data is stored
    ecco2_dir = '/short/m68/kaa561/ROMS-CICE-MCT/data/ECCO2/raw/'
    # Middle of ECCO2 file names
    ecco2_mid = '.1440x720x50.'

    # Variable names for CMIP5
    var_names_cmip5 = ['thetao', 'so', 'uo', 'vo']
    # Corresponding variable names for ECCO2
    var_names_ecco2 = ['THETA', 'SALT', 'UVEL', 'VVEL']
    # Long names for variables to be used as plot titles
    plot_titles = ['Temperature', 'Salinity', 'Eastward Velocity', 'Northward Velocity']
    # Units for each variable
    plot_units = [r'$^{\circ}$C', 'psu', 'm/s', 'm/s']

    # RGB colours to be used for the 39 CMIP5 models
    # These 39 visually distinct (or as visually distinct as possible) colours
    # were generated from http://phrogz.net/css/distinct-colors.html
    cmip5_colours = [(0.83,0,0), (0.49,0,0), (0.32,0,0), (1,0.39,0.39), (0.66,0.26,0.26), (0.83,0.39,0), (0.49,0.23,0), (0.83,0.56,0.33), (0.49,0.33,0.19), (1,0.93,0), (0.66,0.61,0), (0.49,0.47,0.19), (0.19,0.32,0), (0.76,1,0.39), (0.31,0.66,0.26), (0,1,0.33), (0.13,0.32,0.19), (0,0.49,0.39), (0.39,1,0.88), (0,0.24,0.32), (0.33,0.69,0.83), (0.19,0.41,0.49), (0,0.18,0.66), (0.33,0.46,0.83), (0.19,0.27,0.49), (0.2,0.2,0.2), (0.16,0,0.83), (0.1,0,0.49), (0.51,0.39,1), (0.67,0,1), (0.21,0,0.32), (0.39,0.19,0.49), (1,0.39,0.92), (0.66,0.26,0.61), (0.32,0.13,0.29), (0.83,0,0.33), (0.83,0.33,0.53), (0.49,0.19,0.31), (1,0,0)]
    # FontProperties object to have smaller font size for legends
    fontP = FontProperties()
    fontP.set_size('small')

    # Build an array of Model objects, one for each of 39 CMIP5 models
    models = build_model_list()

    # Read ECCO2 latitude, longitude, and depth from the first file
    id = Dataset(ecco2_dir + 'THETA' + ecco2_mid + str(start_year) + '01.nc', 'r')
    ecco2_lat = id.variables['LATITUDE_T'][:]
    ecco2_lon = id.variables['LONGITUDE_T'][:]
    ecco2_depth = id.variables['DEPTH_T'][:]
    id.close()
    # Find the first index north of nbdry, and subtract 1 to find the last
    # index south of nbdry    
    j_max = nonzero(ecco2_lat > nbdry)[0][0]
    j_min = j_max - 1
    # Only save the ECCO2 latitude values at these indices
    ecco2_lat = ecco2_lat[j_min:j_max+1]

    # Loop over variables
    for i in range(len(var_names_cmip5)):

        print 'Plotting ' + plot_titles[i]
        # Select variable names
        var_cmip5 = var_names_cmip5[i]
        var_ecco2 = var_names_ecco2[i]

        print 'Processing ECCO2 data'
        # Create empty array of dimension time x depth x longitude
        ecco2_data = ma.empty([12*(end_year-start_year+1), size(ecco2_depth), size(ecco2_lon)])
        # Initialise next available time index in this array
        posn = 0
        # Loop over years and months
        for year in range(start_year, end_year+1):
            for month in range(12):
                print 'month '+str(month+1)+' of '+str(year)

                # Construct filename
                if (month+1) < 10:
                    month_str = '0' + str(month+1)
                else:
                    month_str = str(month+1)
                ecco2_file = ecco2_dir + str(var_ecco2) + ecco2_mid + str(year) + month_str + '.nc'
                # Read data
                id = Dataset(ecco2_file, 'r')
                data = id.variables[var_ecco2][0,:,j_min:j_max+1,:]
                id.close()

                # Linearly interpolate to nbdry and save to master array
                ecco2_data[posn,:,:] = (data[:,1,:]-data[:,0,:])/(ecco2_lat[1]-ecco2_lat[0])*(nbdry-ecco2_lat[0]) + data[:,0,:]
                posn +=1

        # Zonally average - this is easy on a regular grid
        ecco2_data_zonalavg = mean(ecco2_data, axis=2)
        # Time average - also easy because equally spaced time indices
        ecco2_data_timeavg = mean(ecco2_data_zonalavg, axis=0)

        # Plot ECCO2 data
        figure(figsize=(12,9))
        ax = subplot(111)
        # Keyword "zorder" makes sure this line will be on top of all the
        # thinner CMIP5 lines
        ax.plot(ecco2_data_timeavg, ecco2_depth, label='ECCO2', color='black', linewidth=3, zorder=len(models)+1)

        # Loop through Model objects with a manual counter j
        j = 0
        for model in models:
            print 'Processing ' + model.name

            # Get the model output for this variable, and the model's depth
            # axis (they are all on different grids)
            model_data, model_depth = cmip5_field(model, expt, var_cmip5, start_year, end_year)

            # Check if the output actually exists (if not, cmip5_field will
            # return None)
            if model_data is not None:
                # Zonally average - note all CMIP5 models have regular grids
                model_data_zonalavg = mean(model_data, axis=2)
                # Time average
                model_data_timeavg = mean(model_data_zonalavg, axis=0)

                # Add this model's output to the plot, selecting the correct
                # colour
                if model.name == 'ACCESS1-0':
                    ax.plot(model_data_timeavg, model_depth, label=model.name, color=cmip5_colours[j], linewidth=3)
                else:
                    ax.plot(model_data_timeavg, model_depth, label=model.name, color=cmip5_colours[j])
            # Increment j whether or not this model had any output, so that
            # colour choices for each model will be consistent between variables
            j += 1

        # Configure plot
        title(plot_titles[i])
        grid(True)
        xlabel(plot_units[i])
        ylabel('Depth (m)')
        gca().invert_yaxis()

        # Move the plot over a bit so there's room for the legend beside it
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        # Set font size to small for legend using fontP (previously created)
        # so the massive legend will fit
        ax.legend(loc='center left', bbox_to_anchor=(1,0.5), prop=fontP)
        savefig(var_cmip5 + '.png')


# Command-line interface
if __name__ == "__main__":

    cmip5_ecco2_plot()
    
                     
    
