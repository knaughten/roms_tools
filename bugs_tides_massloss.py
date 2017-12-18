from numpy import *
from matplotlib.pyplot import *

def bugs_tides_massloss (tides_logfile, notides_logfile):

    # Read timeseries with tides
    time = []
    massloss_tides = []
    f = open(tides_logfile, 'r')
    f.readline()
    # Read time values
    for line in f:
        try:
            time.append(float(line))
        except(ValueError):
            break
    # Read total mass loss values
    for line in f:
        try:
            massloss_tides.append(float(line))
        except(ValueError):
            break
    f.close()
    massloss_tides = array(massloss_tides)
    # Get number of time indices
    num_time = len(time)

    # Read timeseries without tides
    massloss_notides = []
    f = open(notides_logfile, 'r')
    f.readline()
    # Skip the time values
    for line in f:
        try:
            tmp = float(line)
        except(ValueError):
            break
    # Read total mass loss values
    for line in f:
        try:
            massloss_notides.append(float(line))
        except(ValueError):
            break
    f.close()
    massloss_notides = array(massloss_notides)
    # Trim to be the same size as the tides simulation
    massloss_notides = massloss_notides[:num_time]
    # Calculate percent difference due to tides
    percent_change = (massloss_tides - massloss_notides)/massloss_notides*100
    print 'Initial change in total mass loss: ' + str(percent_change[0]) + '%'
    print 'Change over first year (excluding initial): ' + str(mean(percent_change[1:74])) + '%'
    print 'Change over last year: ' + str(mean(percent_change[-73:])) + '%'

    # Plot both timeseries
    fig, ax = subplots()
    ax.plot(time, massloss_tides, color='blue', linewidth=1.5, label='RMS tides')
    ax.plot(time, massloss_notides, color='green', linewidth=1.5, label='No tides')
    ax.legend()
    grid(True)
    xlim([time[0]-0.05, time[-1]])
    xlabel('Years', fontsize=14)
    ylabel('Gt/y', fontsize=14)
    title('Total basal mass loss', fontsize=18)
    fig.show()
    fig.savefig('bugs_tides_massloss.png')

    # Plot percent difference
    fig, ax = subplots()
    ax.plot(time, percent_change, linewidth=1.5)
    grid(True)
    xlim([time[0]-0.05, time[-1]])
    xlabel('Years', fontsize=14)
    ylabel('%', fontsize=14)
    title('Percent change in total basal mass loss due to\nRMS tide parameterisation in ice shelf cavities', fontsize=17)
    fig.show()
    fig.savefig('bugs_tides_percent_change.png')


# Command-line interface
if __name__ == "__main__":

    tides_logfile = raw_input("Path to mass loss logfile from simulation with tides: ")
    notides_logfile = raw_input("Path to mass loss logfile from simulation without tides: ")
    bugs_tides_massloss(tides_logfile, notides_logfile)
