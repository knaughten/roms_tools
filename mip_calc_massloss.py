from numpy import *

# Calculate the average mass loss for each ice shelf over 2002-2016 for MetROMS,
# FESOM low-res, and FESOM high-res. Print to the screen along with Rignot's
# estimates for 2003-2008.
# Input:
# roms_logfile = path to ROMS logfile from timeseries_massloss.py
# fesom_logfile_lr, fesom_logfile_hr = paths to FESOM logfiles from
#                   timeseries_massloss.py in the fesomtools repository, for
#                   low-res and high-res respectively
def mip_calc_massloss (roms_logfile, fesom_logfile_lr, fesom_logfile_hr):

    # Year simulations start
    year_start = 1992
    # Years to avearge over
    calc_start = 2002
    calc_end = 2016
    # Number of output steps per year in FESOM
    peryear = 365/5
    # Name of each ice shelf
    names = ['Total Mass Loss', 'Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    # Observed mass loss (Rignot 2013) and uncertainty for each ice shelf, in Gt/y
    obs_massloss = [1325, 1.4, 20.7, 135.4, 155.4, 51.8, 101.2, 97.5, 45.2, 144.9, 4.2, 18.2, 7.9, 90.6, 72.6, 27.2, 35.5, -2, 21.6, 6.3, 3.9, 26.8, 9.7, 47.7]
    obs_massloss_error = [235, 14, 67, 40, 45, 19, 8, 7, 4, 14, 2, 3, 3, 8, 15, 10, 23, 3, 18, 2, 2, 14, 16, 34]
    num_shelves = len(obs_massloss)

    # Read ROMS logfile
    roms_time = []
    f = open(roms_logfile, 'r')
    # Skip the first line (header for time array)
    f.readline()
    for line in f:
        try:
            roms_time.append(float(line))
        except(ValueError):
            # Reached the header for the next variable
            break       
    # Set up array for mass loss values at each ice shelf
    roms_massloss_ts = empty([num_shelves, len(roms_time)])
    index = 0
    # Loop over ice shelves
    while index < num_shelves:
        t = 0
        for line in f:
            try:
                roms_massloss_ts[index, t] = float(line)
                t += 1
            except(ValueError):
                # Reached the header for the next ice shelf
                break
        index +=1
    f.close()
    # Add start year to ROMS time array
    roms_time = array(roms_time) + year_start
    # Average between given years
    t_start = nonzero(roms_time >= calc_start)[0][0]
    if calc_end == 2016:
        t_end = size(roms_time)
    else:
        t_end = nonzero(roms_time >= calc_end+1)[0][0]
    roms_massloss = mean(roms_massloss_ts[:,t_start:t_end], axis=1)

    # Read FESOM timeseries
    # Low-res
    tmp = []
    f = open(fesom_logfile_lr, 'r')
    # Skip the first line (header)
    f.readline()
    # Read total mass loss
    num_time = 0
    for line in f:
        try:
            tmp.append(float(line))
            num_time += 1
        except(ValueError):
            # Reached the header for the next variable
            break
    # Set up array for mass loss values at each ice shelf
    fesom_massloss_ts_lr = empty([num_shelves, num_time])
    # Save the total values in the first index
    fesom_massloss_ts_lr[0,:] = array(tmp)
    # Loop over ice shelves
    index = 1
    while index < num_shelves:
        t = 0
        for line in f:
            try:
                fesom_massloss_ts_lr[index,t] = float(line)
                t += 1
            except(ValueError):
                # Reached the header for the next ice shelf
                break
        index += 1
    f.close()
    # Average between given years
    fesom_massloss_lr = mean(fesom_massloss_ts_lr[:,peryear*(calc_start-year_start):peryear*(calc_end+1-year_start)], axis=1)
    # Repeat for high-res
    tmp = []
    f = open(fesom_logfile_hr, 'r')
    f.readline()
    num_time = 0
    for line in f:
        try:
            tmp.append(float(line))
            num_time += 1
        except(ValueError):
            break
    fesom_massloss_ts_hr = empty([num_shelves, num_time])
    fesom_massloss_ts_hr[0,:] = array(tmp)
    index = 1
    while index < num_shelves:
        t = 0
        for line in f:
            try:
                fesom_massloss_ts_hr[index,t] = float(line)
                t += 1
            except(ValueError):
                break
        index += 1
    f.close()
    fesom_massloss_hr = mean(fesom_massloss_ts_hr[:,peryear*(calc_start-year_start):peryear*(calc_end+1-year_start)], axis=1)

    # Loop over ice shelves
    for index in range(num_shelves):
        print names[index]
        print 'MetROMS: ' + str(roms_massloss[index])
        print 'FESOM low-res: ' + str(fesom_massloss_lr[index])
        print 'FESOM high-res: ' + str(fesom_massloss_hr[index])
        # Find the range of observations
        massloss_low = obs_massloss[index] - obs_massloss_error[index]
        massloss_high = obs_massloss[index] + obs_massloss_error[index]
        print 'Rignot: ' + str(massloss_low) + '-' + str(massloss_high)


# Command-line interface
if __name__ == "__main__":

    roms_logfile = raw_input("Path to ROMS logfile from timeseries_massloss.py: ")
    fesom_logfile_lr = raw_input("Path to FESOM low-res logfile from timeseries_massloss.py: ")
    fesom_logfile_hr = raw_input("Path to FESOM high-res logfile from timeseries_massloss.py: ")
    mip_calc_massloss(roms_logfile, fesom_logfile_lr, fesom_logfile_hr)



