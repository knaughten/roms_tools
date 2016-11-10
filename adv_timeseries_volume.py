from numpy import *
from matplotlib.pyplot import *

# Plot timeseries of total sea ice volume for all the advection experiments.
# Before running this script, you must run timeseries_seaice for each
# experiment.
def adv_timeseries_volume ():

    num_simulations = 6
    # Number of output steps in the simulation (daily averages, 1 leap year)
    num_days = 366
    # Paths to simulation directories (in order of decreasing sea ice volume)
    paths = ['/short/m68/kaa561/ROMS-CICE-MCT/tmproms/run/advection/c4_lowdif/', '/short/m68/kaa561/ROMS-CICE-MCT/tmproms/run/advection/a4_lowdif/', '/short/m68/kaa561/ROMS-CICE-MCT/tmproms/run/advection/u3_lowdif/', '/short/m68/kaa561/ROMS-CICE-MCT/tmproms/run/advection/c4_highdif/', '/short/m68/kaa561/ROMS-CICE-MCT/tmproms/run/advection/a4_highdif/', '/short/m68/kaa561/ROMS-CICE-MCT/tmproms/run/advection/u3_highdif/']
    # Name of timeseries_seaice logfile in each directory
    logfile = 'seaice.log'
    # Abbreviations for each simulation to use in the legend
    labels = ['C4_L', 'A4_L', 'U3_L', 'C4_H', 'A4_H', 'U3_H']
    # Colours for plotting each simulation
    colours = ['k', 'm', 'b', 'r', 'g', 'c']

    # Set up array of volume timeseries for each simulation
    volume = zeros([num_simulations, num_days])
    # Read each logfile
    for sim in range(num_simulations):
        f = open(paths[sim] + logfile, 'r')
        # Skip the time header
        f.readline()
        # Skip the time values
        for line in f:
            try:
                tmp = float(line)
            except(ValueError):
                break
        # Skip the sea ice area values
        for line in f:
            try:
                tmp = float(line)
            except(ValueError):
                break
        # Save the sea ice volume values
        t = 0
        for line in f:
            volume[sim,t] = float(line)
            t += 1
        f.close()

    # Set up time array
    time = arange(num_days)
    # Plot
    fig, ax = subplots()
    for sim in range(num_simulations):
        ax.plot(time, volume[sim,:], label=labels[sim], linewidth=2, color=colours[sim])
    title(r'Total sea ice volume', fontsize=18)
    xlabel('days', fontsize=14)
    ylabel(r'million km$^3$', fontsize=14)
    grid(True)
    ax.legend(loc='lower right')

    fig.savefig('adv_timeseries_volume.png')


# Command-line interface
if __name__ == "__main__":
    adv_timeseries_volume()
    

    
