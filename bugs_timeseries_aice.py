from numpy import *
from matplotlib.pyplot import *

def bugs_timeseries_aice ():

    # File paths
    directories = ['/short/m68/kaa561/metroms_iceshelf/tmproms/run/bug_chapter/no_restoring/', '/short/m68/kaa561/metroms_iceshelf/tmproms/run/bug_chapter/no_restoring_evp/', '/short/m68/kaa561/metroms_iceshelf/tmproms/run/bug_chapter/no_restoring_no_kppmod/', '/short/m68/kaa561/metroms_iceshelf/tmproms/run/bug_chapter/no_restoring_eminusp/']
    log_file = 'seaice.log'
    num_expts = len(directories)
    # Title for each simulation
    expt_names = ['Baseline', 'EVP rheology', 'Original KPP', 'Internal evap']
    # Colours to plot
    expt_colours = ['black', 'red', 'blue', 'green']
    # First year of simulation
    year_start = 1992

    # Read the first logfile
    time = []
    total_area0 = []
    f = open(directories[0] + log_file, 'r')
    # Skip first line (header for time array)
    f.readline()
    # Read time values
    for line in f:
        try:
            time.append(float(line))
        except(ValueError):
            # Reached the header for the next variable
            break
    # Read sea ice area
    for line in f:
        try:
            total_area0.append(float(line))
        except(ValueError):
            break
    f.close()
    # Add start year to time array (already in years)
    time = array(time) + year_start
    # Figure out ticks and labels for time axis
    time_ticks = arange(year_start, time[-1], 2)
    time_labels = []
    for tick in time_ticks:
        time_labels.append(str(int(tick)))
    # Set up array of sea ice maxima in all simulations
    num_time = size(time)
    total_area = zeros([num_expts, num_time])
    total_area[0,:] = total_area0
    # Read other logfiles
    for expt in range(1, num_expts):
        f = open(directories[expt] + log_file, 'r')
        f.readline()
        # Skip time values
        for line in f:
            try:
                tmp = float(line)
            except(ValueError):
                break
        t = 0
        for line in f:
            try:
                total_area[expt,t] = float(line)
                t += 1
            except(ValueError):
                break
        f.close()

    # Plot
    fig, ax = subplots(figsize=(9,7))
    for expt in range(num_expts):
        plot(time, total_area[expt,:], label=expt_names[expt], color=expt_colours[expt], linewidth=1.2)
    xlim([time[0], time[-1]])
    ax.set_xticks(time_ticks)
    ax.set_xticklabels(time_labels)
    xlabel('Year', fontsize=16)
    ylim([0, 20])
    ax.set_yticks(arange(5,20+5,5))
    ylabel(r'million km$^2$', fontsize=16)
    title('Total sea ice area', fontsize=20)
    grid(True)
    ax.legend()
    fig.show()
    fig.savefig('bugs_timeseries_aice.png')


# Command-line interface
if __name__ == "__main__":

    bugs_timeseries_aice()
