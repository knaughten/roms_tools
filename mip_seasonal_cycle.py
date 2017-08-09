from numpy import *

# Calculate the average amplitude of the seasonal cycle in total basal mass
# loss over 2002-2016, in MetROMS, low-res FESOM, and high-res FESOM. Print
# results to the screen.
def mip_seasonal_cycle (roms_logfile, fesom_logfile_lr, fesom_logfile_hr):

    # Year simulations start
    year_start = 1992
    # Years to avearge over
    calc_start = 2002
    calc_end = 2016
    # Number of output steps per year in FESOM
    peryear = 365/5

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
    # Read total mass loss
    roms_massloss = []
    for line in f:
        try:
            roms_massloss.append(float(line))
        except(ValueError):
            # Reached the header for the first individual ice shelf
            break
    f.close()
    # Add start year to ROMS time array
    roms_time = array(roms_time) + year_start
    # Calculate amplitude for each year
    roms_amplitude = []
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
        # Get min and max between these bounds
        curr_min = amin(roms_massloss[t_start:t_end])
        curr_max = amax(roms_massloss[t_start:t_end])
        # Save amplitude
        roms_amplitude.append(curr_max-curr_min)
    # Calculate average amplitude
    roms_val = mean(array(roms_amplitude))
    print 'Average amplitude of seasonal cycle in total basal mass loss: '
    print 'ROMS: ' + str(roms_val)

    # Read FESOM logfiles
    # Low-res
    f = open(fesom_logfile_lr, 'r')
    # Skip the first line (header)
    f.readline()
    # Read total mass loss
    fesom_massloss_lr = []
    for line in f:
        try:
            fesom_massloss_lr.append(float(line))
        except(ValueError):
            # Reached the header for the first individual ice shelf
            break
    f.close()
    # Calculate amplitude for each year
    fesom_amplitude_lr = []
    for year in range(calc_start, calc_end):
        t_start = peryear*(year-year_start)
        t_end = peryear*(year+1-year_start)
        # Get min and max between these bounds
        curr_min = amin(fesom_massloss_lr[t_start:t_end])
        curr_max = amax(fesom_massloss_lr[t_start:t_end])
        # Save amplitude
        fesom_amplitude_lr.append(curr_max-curr_min)
    # Calculate average amplitude
    fesom_val_lr = mean(array(fesom_amplitude_lr))
    print 'FESOM low-res: ' + str(fesom_val_lr)

    # High-res
    f = open(fesom_logfile_hr, 'r')
    f.readline()
    fesom_massloss_hr = []
    for line in f:
        try:
            fesom_massloss_hr.append(float(line))
        except(ValueError):
            break
    f.close()
    fesom_amplitude_hr = []
    for year in range(calc_start, calc_end):
        t_start = peryear*(year-year_start)
        t_end = peryear*(year+1-year_start)
        curr_min = amin(fesom_massloss_hr[t_start:t_end])
        curr_max = amax(fesom_massloss_hr[t_start:t_end])
        fesom_amplitude_hr.append(curr_max-curr_min)
    fesom_val_hr = mean(array(fesom_amplitude_hr))
    print 'FESOM high-res: ' + str(fesom_val_hr)
        

# Command-line interface
if __name__ == "__main__":

    roms_logfile = raw_input("Path to ROMS logfile from timeseries_massloss.py: ")
    fesom_logfile_lr = raw_input("Path to FESOM low-res logfile from timeseries_massloss.py: ")
    fesom_logfile_hr = raw_input("Path to FESOM high-res logfile from timeseries_massloss.py: ")
    mip_seasonal_cycle(roms_logfile, fesom_logfile_lr, fesom_logfile_hr)
    
