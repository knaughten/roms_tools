from re import split
from numpy import *
from matplotlib.pyplot import plot,xlabel,ylabel,clf,show

# Read the volume values written in ocean.log each timestep
# and return them as a 1D array. 
def volume (file_path):

    file = open(file_path, 'r')
    vol_his = []

    # Skip to the part of the file where timestepping begins, i.e.
    # to the line in the file starting with "STEP"
    while True:
        # Split the next line into words based on white space
        words = (file.readline()).split()
        # Check if the first word is STEP
        # If so, exit the loop; otherwise, discard this line and continue 
        if size(words) > 0:
            if words[0] == 'STEP':
                break

    # For every following line, check if it is the first line of output for
    # a given timestep (note each timestep has 2 lines of output)
    # If so, extract and save the volume value            
    for line in file:
        # Split the line into words based on white space
        words = line.split()
        try:
            # If the next line is empty, this will throw an IndexError
            # If the next line doesn't start with a number, this will
            # throw a ValueError
            # This cuts out second lines of output (which start with indices
            # in brackets) and non-output lines (eg messages about reading
            # boundary conditions or writing to history file)
            test = float(words[0])
            # If we are still in the try block, the line is one we want
            # Extract the volume and append to array
            vol = float(words[6])
            vol_his.append(vol)
        except (ValueError, IndexError):
            # It was a line we didn't want; just skip over it
            pass

    file.close()
    return vol_his
        

# Read the volume values from any number of ocean.log files (designed
# for a long simulation split up into several runs) and plot the percent
# anomalies (with respect to the initial value) against time.
def plot_volume (files, dt, freq):

    vol_all = []
    seconds_per_year = 365.0*24.0*60.0*60.0

    # Extract the volume values for each ocean.log file
    for filename in files:
        vol_curr = volume(filename)
        # Append to existing array
        vol_all = concatenate([vol_all, vol_curr])
        print 'Completed file ' + filename

    # Remove singleton dimension
    vol_all = squeeze(vol_all)
    # Calculate percent anomalies
    vol_all = 100*(vol_all - vol_all[0])/vol_all[0]
    # Calculate time array in years
    time = arange(size(vol_all))*dt/seconds_per_year*freq

    # Plot the results
    clf()
    plot(time, vol_all)
    xlabel('Years')
    ylabel('Volume Percent Anomalies')
    grid(True)
    show()


# Command-line interface
if __name__ == "__main__":

    files = []
    filename = raw_input("Path to first ocean.log file: ")
    files.append(filename)

    while True:
        filename = raw_input("Path to next ocean.log file, or enter if finished: ")
        if len(filename) == 0:
            break
        files.append(filename)
        

    dt = double(raw_input("Timestep in seconds: "))
    freq = double(raw_input("Output frequency in timesteps: "))

    plot_volume(files, dt, freq)
        
