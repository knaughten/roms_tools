from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from cartesian_grid_2d import *
from timeseries_massloss import calc_grid

# Plot timeseries of total basal mass loss and area-averaged ice shelf melt
# rates split up into 3 different depth classes for the ice shelf draft. 
# Input:
# file_path = path to ocean history/averages file
# log_path = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
def timeseries_massloss_depth (file_path, log_path):

    # Bounds on depth classes
    draft_min = array([0, 250, 500])
    draft_max = array([250, 500, 3000])
    num_classes = size(draft_min)
    # Labels for legend
    labels = ['<'+str(draft_max[0])+' m']
    for n in range(1, num_classes-1):
        labels.append(str(draft_min[n])+'-'+str(draft_max[n])+' m')
    labels.append('>'+str(draft_min[-1])+' m')
    # Density of ice in kg/m^3
    rho_ice = 916

    old_time = []
    # Check if the log file exists
    if exists(log_path):
        print 'Reading previously calculated values'
        f = open(log_path, 'r')
        # Skip the first line (header for time array)
        f.readline()
        for line in f:
            try:
                old_time.append(float(line))
            except(ValueError):
                # Reached the header for the next variable
                break
        # Set up array for mass loss values for each depth class
        old_massloss = empty([num_classes, len(old_time)])
        n = 0
        # Loop over depth classes
        while n < num_classes:
            t = 0
            for line in f:
                try:
                    old_massloss[n,t] = float(line)
                    t += 1
                except(ValueError):
                    # Reached the header for the next depth class
                    break
            n += 1
        f.close()

    # Calculate dA (masked with ice shelf mask) and lon and lat coordinates
    print 'Analysing grid'
    dA, lon, lat = calc_grid(file_path)

    # Read time data and convert from seconds to years
    id = Dataset(file_path, 'r')
    new_time = id.variables['ocean_time'][:]/(365.25*24*60*60)
    if exists(log_path):
        # Concatenate with time values from log file
        start_t = len(old_time)
        time = old_time
        for t in range(size(new_time)):
            time.append(new_time[t])
        time = array(time)
    else:
        start_t = 0
        time = new_time
    # Set up array of mass loss values
    massloss = empty([num_classes, size(time)])
    if exists(log_path):
        # Fill first start_t timesteps with existing values
        massloss[:,0:start_t] = old_massloss[:,:]

    print 'Reading data'
    # Read ice shelf draft, make it positive
    zice = -1*id.variables['zice'][:-15,1:-1]
    # Read melt rate and convert from m/s to m/y
    ismr = id.variables['m'][:,:-15,1:-1]*365.25*24*60*60
    id.close()

    print 'Setting up arrays'
    # Set up array of masked area values for each depth class
    dA_masked = ma.empty([num_classes, size(dA,0), size(dA,1)])
    # Set up array of conversion factors from mass loss to area-averaged melt
    # rate for each depth class
    factors = empty(num_classes)
    for n in range(num_classes):
        # Mask dA for the current depth class
        dA_tmp = ma.masked_where((zice <= draft_min[n]) + (zice > draft_max[n]), dA)
        dA_masked[n,:,:] = dA_tmp
        # Calculate area of this depth class
        area_tmp = sum(dA_masked[n,:,:])
        print 'Area of ice shelf draft between '+str(draft_min[n])+' and '+str(draft_max[n])+'m: '+str(area_tmp)+' m^2'
        # Calculate conversion factor
        factors[n] = 1e12/(rho_ice*area_tmp)

    # Build timeseries
    print 'Calculating timeseries'
    for t in range(start_t, size(time)):
        # Loop over depth classes
        for n in range(num_classes):
            # Integrate ice shelf melt rate over area to get volume loss
            volumeloss = sum(ismr[t-start_t,:,:]*dA_masked[n,:,:])
            # Convert to mass loss in Gt/y
            massloss[n,t] = 1e-12*rho_ice*volumeloss

    print 'Plotting'

    # Start with mass loss
    fig, ax = subplots(figsize=(10,5))
    # One line for each depth class
    for n in range(num_classes):
        ax.plot(time, massloss[n,:], label=labels[n], linewidth=2)
    # Configure plot
    title('Basal Mass Loss', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel('Gt/y', fontsize=14)
    xlim([time[0], time[-1]])
    grid(True)
    # Move the plot over to make room for legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    # Make legend
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('massloss_depth.png')

    # Repeat for average melt rate
    fig, ax = subplots(figsize=(10,6))
    for n in range(num_classes):
        ax.plot(time, massloss[n,:]*factors[n], label=labels[n], linewidth=2)
    # Configure plot
    title('Area-Averaged Ice Shelf Melt Rate', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel('m/y', fontsize=14)
    grid(True)
    # Move the plot over to make room for legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    # Make legend
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    fig.savefig('ismr_depth.png')

    print 'Saving results to log file'
    f = open(log_path, 'w')
    f.write('Time (years):\n')
    for t in range(size(time)):
        f.write(str(time[t]) + '\n')
    f.write('Basal Mass Loss for ice shelf drafts <' + str(draft_max[0]) + ' m:\n')
    for t in range(size(time)):
        f.write(str(massloss[0,t]) + '\n')
    for n in range(1, num_classes-1):
        f.write('Basal Mass Loss for ice shelf drafts ' + str(draft_min[n]) + '-' + str(draft_max[n]) + ' m:\n')
        for t in range(size(time)):
            f.write(str(massloss[n,t]) + '\n')
    f.write('Basal Mass Loss for ice shelf drafts >' + str(draft_min[-1]) + 'm:\n')
    for t in range(size(time)):
        f.write(str(massloss[-1,t]) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input('Enter path to ocean history/averages file: ')
    log_path = raw_input('Enter path to log file to save values and/or read previously calculated values: ')
    timeseries_massloss_depth(file_path, log_path)
