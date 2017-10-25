from numpy import *

# Calculate the average amplitude of the seasonal cycle in total basal mass
# loss over 2002-2016, in MetROMS, low-res FESOM, and high-res FESOM, for the
# entire continent as well as the Amery. Print results to the screen.
def mip_seasonal_cycle (roms_logfile, fesom_logfile_lr, fesom_logfile_hr):

    # Year simulations start
    year_start = 1992
    # Years to avearge over
    calc_start = 2002
    calc_end = 2016
    # Number of output steps per year in FESOM
    peryear = 365/5
    # Name of each ice shelf
    names = ['Total Mass Loss', 'Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    num_shelves = len(names)
    i_total = 0
    i_amery = 16

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
    roms_massloss = empty([num_shelves, len(roms_time)])
    index = 0
    # Loop over ice shelves
    while index < num_shelves:
        t = 0
        for line in f:
            try:
                roms_massloss[index, t] = float(line)
                t += 1
            except(ValueError):
                # Reached the header for the next ice shelf
                break
        index +=1
    f.close()
    # Add start year to ROMS time array
    roms_time = array(roms_time) + year_start
    # Calculate amplitude for each year
    #roms_amplitude_total = []
    roms_min_total = []
    roms_max_total = []
    #roms_amplitude_amery = []
    roms_min_amery = []
    roms_max_amery = []
    for year in range(calc_start, calc_end):
        # Find the first index after the beginning of this year
        t_start = nonzero(roms_time >= year)[0][0]
        # Find the first index after the beginning of next year
        tmp2 = nonzero(roms_time >= year+1)[0]
        if len(tmp2)==0:
            # No such index, but we might not have run out of data
            # eg simulation that ends on 31 December 2016
            t_end = len(roms_time)
        else:
            t_end = tmp2[0]
        # Total Antarctica
        # Get min and max between these bounds
        curr_min = amin(roms_massloss[i_total,t_start:t_end])
        curr_max = amax(roms_massloss[i_total,t_start:t_end])
        roms_min_total.append(curr_min)
        roms_max_total.append(curr_max)
        # Save amplitude
        #roms_amplitude_total.append(curr_max-curr_min)
        # Amery
        curr_min = amin(roms_massloss[i_amery,t_start:t_end])
        curr_max = amax(roms_massloss[i_amery,t_start:t_end])
        #roms_amplitude_amery.append(curr_max-curr_min)
        roms_min_amery.append(curr_min)
        roms_max_amery.append(curr_max)
    # Calculate average amplitude
    #roms_val_total = mean(array(roms_amplitude_total))
    #roms_val_amery = mean(array(roms_amplitude_amery))
    #print 'Average amplitude of seasonal cycle in total basal mass loss: '
    #print 'ROMS (Total Antarctica): ' + str(roms_val_total)
    #print 'ROMS (Amery): ' + str(roms_val_amery)
    print 'Average min/max in total basal mass loss: '
    print 'ROMS (Total Antarctica): ' + str(mean(array(roms_min_total))) + ' Gt/y min, ' + str(mean(array(roms_max_total))) + ' Gt/y max'
    print 'ROMS (Amery): ' + str(mean(array(roms_min_amery))) + ' Gt/y min, ' + str(mean(array(roms_max_amery))) + ' Gt/y max'

    # Read FESOM logfiles
    # Low-res
    f = open(fesom_logfile_lr, 'r')
    # Skip the first line (header)
    f.readline()
    # Read total mass loss
    tmp = []
    num_time = 0
    for line in f:
        try:
            tmp.append(float(line))
            num_time += 1
        except(ValueError):
            # Reached the header for the first individual ice shelf
            break
    # Set up array for mass loss values at each ice shelf
    fesom_massloss_lr = empty([num_shelves, num_time])
    # Save the total values in the first index
    fesom_massloss_lr[0,:] = array(tmp)
    # Loop over ice shelves
    index = 1
    while index < num_shelves:
        t = 0
        for line in f:
            try:
                fesom_massloss_lr[index,t] = float(line)
                t += 1
            except(ValueError):
                # Reached the header for the next ice shelf
                break
        index += 1
    f.close()
    # Calculate amplitude for each year
    #fesom_amplitude_total_lr = []
    fesom_min_total_lr = []
    fesom_max_total_lr = []
    #fesom_amplitude_amery_lr = []
    fesom_min_amery_lr = []
    fesom_max_amery_lr = []
    for year in range(calc_start, calc_end):
        t_start = peryear*(year-year_start)
        t_end = peryear*(year+1-year_start)
        # Total Antarctica
        # Get min and max between these bounds
        curr_min = amin(fesom_massloss_lr[i_total,t_start:t_end])
        curr_max = amax(fesom_massloss_lr[i_total,t_start:t_end])
        # Save amplitude
        #fesom_amplitude_total_lr.append(curr_max-curr_min)
        fesom_min_total_lr.append(curr_min)
        fesom_max_total_lr.append(curr_max)
        # Amery
        curr_min = amin(fesom_massloss_lr[i_amery,t_start:t_end])
        curr_max = amax(fesom_massloss_lr[i_amery,t_start:t_end])
        #fesom_amplitude_amery_lr.append(curr_max-curr_min)
        fesom_min_amery_lr.append(curr_min)
        fesom_max_amery_lr.append(curr_max)
    # Calculate average amplitude
    #fesom_val_total_lr = mean(array(fesom_amplitude_total_lr))
    #fesom_val_amery_lr = mean(array(fesom_amplitude_amery_lr))
    #print 'FESOM low-res (Total Antarctica): ' + str(fesom_val_total_lr)
    #print 'FESOM low-res (Amery): ' + str(fesom_val_amery_lr)
    print 'FESOM low-res (Total Antarctica): ' + str(mean(array(fesom_min_total_lr))) + ' Gt/y min, ' + str(mean(array(fesom_max_total_lr))) + ' Gt/y max'
    print 'FESOM low-res (Amery): ' + str(mean(array(fesom_min_amery_lr))) + ' Gt/y min, ' + str(mean(array(fesom_max_amery_lr))) + ' Gt/y max'

    # High-res
    f = open(fesom_logfile_hr, 'r')
    f.readline()
    tmp = []
    num_time = 0
    for line in f:
        try:
            tmp.append(float(line))
            num_time += 1
        except(ValueError):
            break
    fesom_massloss_hr = empty([num_shelves, num_time])
    fesom_massloss_hr[0,:] = array(tmp)
    index = 1
    while index < num_shelves:
        t = 0
        for line in f:
            try:
                fesom_massloss_hr[index,t] = float(line)
                t += 1
            except(ValueError):
                break
        index += 1
    f.close()
    #fesom_amplitude_total_hr = []
    fesom_min_total_hr = []
    fesom_max_total_hr = []
    #fesom_amplitude_amery_hr = []
    fesom_min_amery_hr = []
    fesom_max_amery_hr = []
    for year in range(calc_start, calc_end):
        t_start = peryear*(year-year_start)
        t_end = peryear*(year+1-year_start)
        curr_min = amin(fesom_massloss_hr[i_total,t_start:t_end])
        curr_max = amax(fesom_massloss_hr[i_total,t_start:t_end])
        fesom_min_total_hr.append(curr_min)
        fesom_max_total_hr.append(curr_max)
        #fesom_amplitude_total_hr.append(curr_max-curr_min)
        curr_min = amin(fesom_massloss_hr[i_amery,t_start:t_end])
        curr_max = amax(fesom_massloss_hr[i_amery,t_start:t_end])
        #fesom_amplitude_amery_hr.append(curr_max-curr_min)
        fesom_min_amery_hr.append(curr_min)
        fesom_max_amery_hr.append(curr_max)
    #fesom_val_total_hr = mean(array(fesom_amplitude_total_hr))
    #fesom_val_amery_hr = mean(array(fesom_amplitude_amery_hr))
    #print 'FESOM high-res (Total Antarctica): ' + str(fesom_val_total_hr)
    #print 'FESOM high-res (Amery): ' + str(fesom_val_amery_hr)
    print 'FESOM high-res (Total Antarctica): ' + str(mean(array(fesom_min_total_hr))) + ' Gt/y min, ' + str(mean(array(fesom_max_total_hr))) + ' Gt/y max'
    print 'FESOM high-res (Amery): ' + str(mean(array(fesom_min_amery_hr))) + ' Gt/y min, ' + str(mean(array(fesom_max_amery_hr))) + ' Gt/y max'
        

# Command-line interface
if __name__ == "__main__":

    roms_logfile = raw_input("Path to ROMS logfile from timeseries_massloss.py: ")
    fesom_logfile_lr = raw_input("Path to FESOM low-res logfile from timeseries_massloss.py: ")
    fesom_logfile_hr = raw_input("Path to FESOM high-res logfile from timeseries_massloss.py: ")
    mip_seasonal_cycle(roms_logfile, fesom_logfile_lr, fesom_logfile_hr)
    
