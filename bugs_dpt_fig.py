from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from mip_timeseries import roms_annual_avg

def bugs_dpt_fig (dpt_log_laplacian, dpt_log_biharmonic):

    year_start = 1992
    
    # Read Laplacian timeseries
    time = []
    dpt_laplacian = []
    f = open(dpt_log_laplacian, 'r')
    # Skip first line (header for time array)
    f.readline()
    for line in f:
        try:
            time.append(float(line))
        except(ValueError):
            # Reached the header for the next variable
            break
    for line in f:
        dpt_laplacian.append(float(line))
    f.close()
    # Read biharmonic timeseries
    dpt_biharmonic = []
    f = open(dpt_log_biharmonic, 'r')
    f.readline()
    # Skip the time array, as it's the same
    for line in f:
        try:
            tmp = float(line)
        except(ValueError):
            break
    for line in f:
        dpt_biharmonic.append(float(line))
    f.close()
    # Make sure they're the same length
    if len(dpt_laplacian) != len(dpt_biharmonic):
        print 'These logfiles do not cover the same time periods!'
        exit
    # Add start year to time array
    time = array(time) + year_start
    # Annually average DPT
    dpt_laplacian_avg = roms_annual_avg(time, dpt_laplacian, year_start)
    dpt_biharmonic_avg = roms_annual_avg(time, dpt_biharmonic, year_start)
    # Make new time array
    time_annual = arange(len(dpt_laplacian_avg)) + year_start
    # Plot
    fig, ax = subplots(figsize=(10,6))
    ax.plot(time_annual, dpt_laplacian_avg, label='Laplacian', color='green', linewidth=2)
    ax.plot(time_annual, dpt_biharmonic_avg, label='Biharmonic', color='blue', linewidth=2)
    title('Drake Passage Transport (annually averaged)', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel('Sv', fontsize=14)
    xlim([time_annual[0], time_annual[-1]])
    grid(True)
    # Move plot over to make room for legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    # Make legend
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.show()
    fig.savefig('bugs_dpt.png')


# Command-line interface
if __name__ == "__main__":

    dpt_log_laplacian = raw_input("Path to Drake Passage Transport logfile for Laplacian viscosity simulation: ")
    dpt_log_biharmonic = raw_input("Path to Drake Passage Transport logfile for baseline simulation: ")
    bugs_dpt_fig(dpt_log_laplacian, dpt_log_biharmonic)

    
