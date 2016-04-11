from cmip5_paths import *
from cmip5_field_ss import *
from aviso_field import *
from matplotlib.pyplot import *
from numpy import *
from netCDF4 import Dataset

# Compare AVISO sea surface height reanalyses and CMIP5 output, time-averaged 
# between 1995 and 2005 and zonally averaged at the northern boundary of the 
# ROMS grid, with model on the x-axis and sea surface height on the y-axis.
def cmip5_aviso_plot ():

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

    # Start an array of values (one for each model)
    labels = ['AVISO']
    zeta = [aviso_data_timeavg]

    # Loop through Model objects
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
            # Append model name and value
            labels.append(model.name)
            zeta.append(model_data_timeavg)

    # Rerrange labels and zeta so that ACCESS1-0 is right after AVISO
    #if labels[-1] == 'ACCESS1-0':
        #labels = [labels[0]] + [labels[-1]] + labels[1:-1]

    # Plot
    figure(figsize=(12,9))
    ax = subplot(111)
    ax.plot(arange(len(zeta)), zeta, 'ko')

    # Configure plot
    title('Sea Surface Height (m)')
    grid(True)
    ylabel('Depth (m)')
    xlim([-1, len(labels)])
    # Move the plot upward so there's room for model labels below
    box = ax.get_position()
    ax.set_position([box.x0, box.y0+box.height*0.1, box.width, box.height*0.9])
    # Add model labels
    xticks(arange(len(labels)), labels, rotation=-90)    
    
    savefig('zos.png')


# Command-line interface
if __name__ == "__main__":

    cmip5_aviso_plot()
    
                     
    
