# Plot timeseries of the maximum |u| and |v| in a ROMS history file
# Input: file_path = path to history file (in single quotes)

from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import plot,xlabel,ylabel,legend,clf,show

def max_vel (file_path):

    # Read variables from the history file
    file = Dataset(file_path, 'r')
    time = file.variables['ocean_time'][:]
    # Convert time from seconds to years
    time = time/(365*24*60*60)
    max_u = []
    max_v = []

    for l in range(size(time)):
        print 'Processing timestep ' + str(l+1) + ' of ' + str(size(time))
        u = file.variables['u'][l,:,:-15,:-1]
        v = file.variables['v'][l,:,:-15,1:-1]
        # Find the maximum |u| and |v| at this timestep
        max_u.append(amax(abs(u)))
        max_v.append(amax(abs(v)))

    file.close()

    # Plot these two timeseries
    clf()
    plot(time, max_u, label='max u')
    plot(time, max_v, label='max v')
    xlabel('Years')
    ylabel('m/s')
    legend(loc='lower right')
    show()


if __name__ == "__main__":

    file_path = raw_input("Path to ocean history file: ")
    max_vel(file_path)
