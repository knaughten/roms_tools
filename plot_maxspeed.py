from re import split
from numpy import *
from matplotlib.pyplot import *

# Read the maximum speed values written in ocean.log each timestep
# and return them as a 1D array. 
def maxspeed (file_path):

    file = open(file_path, 'r')
    maxv_his = []

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

    # For every following line, check if it is the second line of output for
    # a given timestep (note each timestep has 2 lines of output)
    # If so, extract and save the maximum speed value
    save_next = False
    for line in file:
        # Split the line into words based on white space
        words = line.split()
        if save_next:
            # The previous line was a first line of output
            # This is a second line of output
            maxv = float(words[4])
            maxv_his.append(maxv)
            save_next = False
        else:
            try:
                # If the next line is empty, this will throw an IndexError
                # If the next line doesn't start with a number, this will
                # throw a ValueError
                # This cuts out non-output lines (eg messages about reading
                # boundary conditions or writing to history file)
                test = float(words[0])
                # If we are still in the try block, the line is one we want
                # Boolean flag to save the next line of output
                save_next = True
            except (ValueError, IndexError):
                # It was a line we didn't want; just skip over it
                pass

    file.close()
    return maxv_his

# Read the kinetic energy values from any number of ocean.log files (designed
# for a long simulation split up into several runs) and plot them against time
def plot_maxspeed (files, dt, freq):

    maxv_all = []
    seconds_per_year = 365.0*24.0*60.0*60.0

    # Extract the kinetic energy values for each ocean.log file
    for filename in files:
        maxv_curr = maxspeed(filename)
        # Append to existing array
        maxv_all = concatenate([maxv_all, maxv_curr])
        print 'Completed file ' + filename

    # Remove singleton dimension
    maxv_all = squeeze(maxv_all)
    # Calculate time array in years
    time = arange(size(maxv_all))*dt/seconds_per_year*freq

    # Plot the results
    figure()
    plot(time, maxv_all)
    xlabel('Years')
    ylabel('Maximum speed')
    grid(True)
    show(block=False)


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

    plot_maxspeed(files, dt, freq)
