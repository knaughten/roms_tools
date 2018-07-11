from numpy import *

def mip_std_massloss (roms_logfile, roms_logfile_bs, fesom_logfile_lr, fesom_logfile_bs_lr, fesom_logfile_hr, fesom_logfile_bs_hr):

    year_start = 1992
    year_end = 2016
    # First year to consider
    calc_start = 2002
    # Days per output in FESOM
    days_per_output = 5
    # Name of each ice shelf
    names = ['Total Mass Loss', 'Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange & Bach Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    num_shelves = len(names)
    # Some Bellingshausen ice shelves were split up later
    names_bs = ['Wilkins Ice Shelf', 'Stange Ice Shelf', 'George VI Ice Shelf']
    num_shelves_bs = len(names_bs)

    num_years = year_end-year_start+1
    num_years_calc = year_end-calc_start+1

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
    # Add start year to time array
    roms_time = array(roms_time) + year_start
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
    # Repeat for Bellingshausen
    f = open(roms_logfile_bs, 'r')
    f.readline()
    # Skip the time values (should be the same)
    for line in f:
        try:
            tmp = float(line)
        except(ValueError):
            # Reached the header for the next variable
            break       
    roms_massloss_bs = empty([num_shelves_bs, len(roms_time)])
    index = 0
    while index < num_shelves_bs:
        t = 0
        for line in f:
            try:
                roms_massloss_bs[index, t] = float(line)
                t += 1
            except(ValueError):
                break
        index +=1

    # Read FESOM logfiles
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
    # Repeat for Bellingshausen
    f = open(fesom_logfile_bs_lr, 'r')
    f.readline()
    fesom_massloss_bs_lr = empty([num_shelves_bs, num_time])
    index = 0
    while index < num_shelves_bs:
        t = 0
        for line in f:
            try:
                fesom_massloss_bs_lr[index,t] = float(line)
                t += 1
            except(ValueError):
                break
        index += 1
    f.close()

    # High-res
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
    f = open(fesom_logfile_bs_hr, 'r')
    f.readline()
    fesom_massloss_bs_hr = empty([num_shelves_bs, num_time])
    index = 0
    while index < num_shelves_bs:
        t = 0
        for line in f:
            try:
                fesom_massloss_bs_hr[index,t] = float(line)
                t += 1
            except(ValueError):
                break
        index += 1
    f.close()

    # Annually average
    # ROMS
    roms_massloss_avg = empty([num_shelves, num_years])
    for index in range(num_shelves):
        roms_massloss_avg[index,:] = roms_annual_avg(roms_time, roms_massloss[index,:], year_start)
    roms_massloss_bs_avg = empty([num_shelves_bs, num_years])
    for index in range(num_shelves_bs):
        roms_massloss_bs_avg[index,:] = roms_annual_avg(roms_time, roms_massloss_bs[index,:], year_start)
    # Low-res FESOM
    fesom_massloss_lr_avg = empty([num_shelves, num_years])
    for index in range(num_shelves):
        fesom_massloss_lr_avg[index,:] = fesom_annual_avg(fesom_massloss_lr[index,:], days_per_output)
    fesom_massloss_bs_lr_avg = empty([num_shelves_bs, num_years])
    for index in range(num_shelves_bs):
        fesom_massloss_bs_lr_avg[index,:] = fesom_annual_avg(fesom_massloss_bs_lr[index,:], days_per_output)
    # High-res FESOM
    fesom_massloss_hr_avg = empty([num_shelves, num_years])
    for index in range(num_shelves):
        fesom_massloss_hr_avg[index,:] = fesom_annual_avg(fesom_massloss_hr[index,:], days_per_output)
    fesom_massloss_bs_hr_avg = empty([num_shelves_bs, num_years])
    for index in range(num_shelves_bs):
        fesom_massloss_bs_hr_avg[index,:] = fesom_annual_avg(fesom_massloss_bs_hr[index,:], days_per_output)
    # Slice off the years we don't care about
    roms_massloss_avg = roms_massloss_avg[:,calc_start-year_start:]
    roms_massloss_bs_avg = roms_massloss_bs_avg[:,calc_start-year_start:]
    fesom_massloss_lr_avg = fesom_massloss_lr_avg[:,calc_start-year_start:]
    fesom_massloss_bs_lr_avg = fesom_massloss_bs_lr_avg[:,calc_start-year_start:]
    fesom_massloss_hr_avg = fesom_massloss_hr_avg[:,calc_start-year_start:]
    fesom_massloss_bs_hr_avg = fesom_massloss_bs_hr_avg[:,calc_start-year_start:]

    # Loop over ice shelves
    for index in range(num_shelves):
        print names[index]
        # Calculate standard deviation as percent of annual mean
        print 'MetROMS: ' + str(std(roms_massloss_avg[index])/mean(roms_massloss_avg[index])*100) + '%'
        print 'FESOM low-res: ' + str(std(fesom_massloss_lr_avg[index])/mean(fesom_massloss_lr_avg[index])*100) + '%'
        print 'FESOM high-res: ' + str(std(fesom_massloss_hr_avg[index])/mean(fesom_massloss_hr_avg[index])*100) + '%'
        #print 'MetROMS: ' + str((median(roms_massloss_avg[index])-mean(roms_massloss_avg[index]))/mean(roms_massloss_avg[index])*100)
        #print 'FESOM low-res: ' + str((median(fesom_massloss_lr_avg[index])-mean(fesom_massloss_lr_avg[index]))/mean(fesom_massloss_lr_avg[index])*100)
        #print 'FESOM high-res: ' + str((median(fesom_massloss_hr_avg[index])-mean(fesom_massloss_hr_avg[index]))/mean(fesom_massloss_hr_avg[index])*100)
        #print 'MetROMS: mean ' + str(mean(roms_massloss_avg[index])) + ', median ' + str(median(roms_massloss_avg[index]))
        #print 'FESOM low-res: mean ' + str(mean(fesom_massloss_lr_avg[index])) + ', median ' + str(median(fesom_massloss_lr_avg[index]))
        #print 'FESOM high-res: mean ' + str(mean(fesom_massloss_hr_avg[index])) + ', median ' + str(median(fesom_massloss_hr_avg[index]))

    # Repeat for Bellingshausen ice shelves
    for index in range(num_shelves_bs):
        print names_bs[index]
        print 'MetROMS: ' + str(std(roms_massloss_bs_avg[index])/mean(roms_massloss_bs_avg[index])*100) + '%'
        print 'FESOM low-res: ' + str(std(fesom_massloss_bs_lr_avg[index])/mean(fesom_massloss_bs_lr_avg[index])*100) + '%'
        print 'FESOM high-res: ' + str(std(fesom_massloss_bs_hr_avg[index])/mean(fesom_massloss_bs_hr_avg[index])*100) + '%'
        #print 'MetROMS: ' + str((median(roms_massloss_bs_avg[index])-mean(roms_massloss_bs_avg[index]))/mean(roms_massloss_bs_avg[index])*100)
        #print 'FESOM low-res: ' + str((median(fesom_massloss_bs_lr_avg[index])-mean(fesom_massloss_bs_lr_avg[index]))/mean(fesom_massloss_bs_lr_avg[index])*100)
        #print 'FESOM high-res: ' + str((median(fesom_massloss_bs_hr_avg[index])-mean(fesom_massloss_bs_hr_avg[index]))/mean(fesom_massloss_bs_hr_avg[index])*100)
        #print 'MetROMS: mean ' + str(mean(roms_massloss_bs_avg[index])) + ', median ' + str(median(roms_massloss_bs_avg[index]))
        #print 'FESOM low-res: mean ' + str(mean(fesom_massloss_bs_lr_avg[index])) + ', median ' + str(median(fesom_massloss_bs_lr_avg[index]))
        #print 'FESOM high-res: mean ' + str(mean(fesom_massloss_bs_hr_avg[index])) + ', median ' + str(median(fesom_massloss_bs_hr_avg[index]))

    
# Annually average the given ROMS data. Note that ROMS takes leap years into
# account for its output, i.e. the 5-day averaging period holds true throughout
# the entire simulation with no days skipped.
# Input:
# time = 1D array of ROMS time values, in years, with year_start added (eg
#        1992.001, 1993.50) assuming year length of 365.25 days
# data = 1D array of ROMS data corresponding to time array
# year_start = integer containing the first year to process
# Output: data_avg = 1D array of annually averaged data
def roms_annual_avg (time, data, year_start):

    data_avg = []
    year = year_start
    # Loop over years until we run out of data
    while True:
        # Find the first index after the beginning of this year
        tmp1 = nonzero(time >= year)[0]
        if len(tmp1)==0:
            # No such index, we have run out of data
            break
        t_start = tmp1[0]
        # Find the first index after the beginning of next year
        tmp2 = nonzero(time >= year+1)[0]
        if len(tmp2)==0:
            # No such index, but we might not have run out of data
            # eg simulation that ends on 31 December 2016
            t_end = len(time)
        else:
            t_end = tmp2[0]
        if t_end - t_start < 73:
            # This is a partial year, don't count it
            break
        # Calculate mean between these bounds
        data_avg.append(mean(array(data[t_start:t_end])))
        # Increment year
        year += 1
    return data_avg


# Annually average the given FESOM data. Note that FESOM neglects leap years
# in its output, i.e. the 5-day averaging period will skip 1 day every 4 years,
# so every year has exactly 365/5 = 73 batches of output.
# Input:
# data = 1D array of FESOM data 
# days_per_output = averaging period for data
# Output: data_avg = 1D array of annually averaged data
def fesom_annual_avg (data, days_per_output):

    peryear = 365/days_per_output
    data_avg = []
    year = 0
    # Loop over years
    while True:
        if len(data) < peryear*(year+1):
            # Run out of data
            # FESOM output comes in annual files so we don't have to worry about
            # partial years like in ROMS
            break
        # Calculate mean for this year
        data_avg.append(mean(array(data[peryear*year:peryear*(year+1)])))
        # Increment year
        year += 1
    return data_avg


# Command-line interface
if __name__ == "__main__":

    roms_logfile = raw_input("Path to ROMS logfile from timeseries_massloss.py: ")
    roms_logfile_bs = raw_input("Path to ROMS logfile from timeseries_massloss_bellingshausen.py: ")
    fesom_logfile_lr = raw_input("Path to FESOM low-res logfile from timeseries_massloss.py: ")
    fesom_logfile_bs_lr = raw_input("Path to FESOM low-res logfile from timeseries_massloss_bellingshausen.py: ")
    fesom_logfile_hr = raw_input("Path to FESOM high-res logfile from timeseries_massloss.py: ")
    fesom_logfile_bs_hr = raw_input("Path to FESOM high-res logfile from timeseries_massloss_bellingshausen.py: ")
    mip_std_massloss(roms_logfile, roms_logfile_bs, fesom_logfile_lr, fesom_logfile_bs_lr, fesom_logfile_hr, fesom_logfile_bs_hr)
