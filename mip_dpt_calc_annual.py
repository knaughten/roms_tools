from numpy import *
from scipy.stats import linregress

# Calculate the mean Drake Passage transport over 2002-2016, as well as the
# linear trend and standard deviation of annual averages, for MetROMS, low-res
# FESOM, and high-res FESOM. Print the results to the screen.
# Input:
# roms_log = logfile from timeseries_dpt.py for MetROMS
# fesom_log_low, fesom_log_high = logfiles from timeseries_dpt.py in the
#                                 fesomtools repository, for FESOM low-res and
#                                 high-res respectively
def mip_dpt_calc_annual (roms_log, fesom_log_low, fesom_log_high):

    # Averaging period (days)
    days_per_output = 5
    # Years of simulation
    year_start = 1992
    year_end = 2016 
    # Year calculation starts
    calc_start = 2002

    num_years = year_end-year_start+1
    num_years_calc = year_end-calc_start+1

    # Read ROMS timeseries
    roms_time = []
    roms_dpt = []
    f = open(roms_log, 'r')
    # Skip first line (header for time array)
    f.readline()
    for line in f:
        try:
            roms_time.append(float(line))
        except(ValueError):
            # Reached the header for the next variable
            break
    for line in f:
        roms_dpt.append(float(line))
    f.close()
    # Add start year to ROMS time array
    roms_time = array(roms_time) + year_start
    roms_dpt = array(roms_dpt)

    # Read FESOM low-res timeseries
    fesom_dpt_low = []
    f = open(fesom_log_low, 'r')
    # Skip header
    f.readline()
    for line in f:
        fesom_dpt_low.append(float(line))
    f.close()
    # Read FESOM high-res timeseries
    fesom_dpt_high = []
    f = open(fesom_log_high, 'r')
    f.readline()
    for line in f:
        fesom_dpt_high.append(float(line))
    f.close()
    # Make FESOM time array (note that FESOM neglects leap years in its output)
    fesom_time = arange(len(fesom_dpt_low))*days_per_output/365. + year_start
    fesom_dpt_low = array(fesom_dpt_low)
    fesom_dpt_high = array(fesom_dpt_high)

    # Annually average
    roms_dpt_avg = roms_annual_avg(roms_time, roms_dpt, year_start)
    fesom_dpt_low_avg = fesom_annual_avg(fesom_dpt_low, days_per_output)
    fesom_dpt_high_avg = fesom_annual_avg(fesom_dpt_high, days_per_output)
    # Slice off the years we don't care about
    roms_dpt_avg = roms_dpt_avg[calc_start-year_start:]
    fesom_dpt_low_avg = fesom_dpt_low_avg[calc_start-year_start:]
    fesom_dpt_high_avg = fesom_dpt_high_avg[calc_start-year_start:]

    # Calculate and print averages and standard deviations
    print 'Average Drake Passage Transport'
    print 'MetROMS: ' + str(mean(roms_dpt_avg)) + ' +/- ' + str(std(roms_dpt_avg))
    print 'FESOM low-res: ' + str(mean(fesom_dpt_low_avg)) + ' +/- ' + str(std(fesom_dpt_low_avg))
    print 'FESOM high-res: ' + str(mean(fesom_dpt_high_avg)) + ' +/- ' + str(std(fesom_dpt_high_avg))

    # Calculate and print trends
    # Also print p-values so we can see if it's statistically significant
    annual_time = arange(calc_start, year_end+1)
    print 'Trends in Drake Passage Transport'
    slope, intercept, r_value, p_value, std_err = linregress(annual_time, roms_dpt_avg)
    print 'MetROMS: ' + str(slope) + ' Sv/y, p=' + str(p_value)
    slope, intercept, r_value, p_value, std_err = linregress(annual_time, fesom_dpt_low_avg)
    print 'FESOM low-res: ' + str(slope) + ' Sv/y, p=' + str(p_value)
    slope, intercept, r_value, p_value, std_err = linregress(annual_time, fesom_dpt_high_avg)
    print 'FESOM high-res: ' + str(slope) + ' Sv/y, p=' + str(p_value)


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

    roms_log = raw_input("Path to ROMS logfile from timeseries_dpt.py: ")
    fesom_log_low = raw_input("Path to FESOM low-res logfile from timeseries_dpt.py: ")
    fesom_log_high = raw_input("Path to FESOM high-res logfile from timeseries_dpt.py: ")
    mip_dpt_calc_annual(roms_log, fesom_log_low, fesom_log_high)
