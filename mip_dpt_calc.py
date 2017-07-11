from numpy import *
from scipy.stats import linregress

# Calculate the mean Drake Passage transport over 2002-2016, as well as the
# linear trend, for MetROMS, low-res FESOM, and high-res FESOM. Print the
# results to the screen.
# Input:
# roms_log = logfile from timeseries_dpt.py for MetROMS
# fesom_log_low, fesom_log_high = logfiles from timeseries_dpt.py in the
#                                 fesomtools repository, for FESOM low-res and
#                                 high-res respectively
def mip_dpt_calc (roms_log, fesom_log_low, fesom_log_high):

    # Averaging period (days)
    days_per_output = 5
    # Year simulation starts
    year_start = 1992
    # Year calculation starts
    calc_start = 2002

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

    # Find range of time indices to consider
    # ROMS
    t_start_roms = nonzero(roms_time >= calc_start)[0][0]
    t_end_roms = len(roms_time)
    # FESOM
    t_start_fesom = (calc_start-year_start)*365/days_per_output
    t_end_fesom = len(fesom_time)
    # Slice off the indices we don't care about
    roms_time = roms_time[t_start_roms:t_end_roms]
    roms_dpt = roms_dpt[t_start_roms:t_end_roms]
    fesom_time = fesom_time[t_start_fesom:t_end_fesom]
    fesom_dpt_low = fesom_dpt_low[t_start_fesom:t_end_fesom]
    fesom_dpt_high = fesom_dpt_high[t_start_fesom:t_end_fesom]

    # Calculate and print averages
    print 'Average Drake Passage Transport'
    print 'MetROMS: ' + str(mean(roms_dpt))
    print 'FESOM low-res: ' + str(mean(fesom_dpt_low))
    print 'FESOM high-res: ' + str(mean(fesom_dpt_high))

    # Calculate and print trends
    # Also print p-values so we can see if it's statistically significant
    print 'Trends in Drake Passage Transport'
    slope, intercept, r_value, p_value, std_err = linregress(roms_time, roms_dpt)
    print 'MetROMS: ' + str(slope) + ' Sv/y, p=' + str(p_value)
    slope, intercept, r_value, p_value, std_err = linregress(fesom_time, fesom_dpt_low)
    print 'FESOM low-res: ' + str(slope) + ' Sv/y, p=' + str(p_value)
    slope, intercept, r_value, p_value, std_err = linregress(fesom_time, fesom_dpt_high)
    print 'FESOM high-res: ' + str(slope) + ' Sv/y, p=' + str(p_value)
    

# Command-line interface
if __name__ == "__main__":

    roms_log = raw_input("Path to ROMS logfile from timeseries_dpt.py: ")
    fesom_log_low = raw_input("Path to FESOM low-res logfile from timeseries_dpt.py: ")
    fesom_log_high = raw_input("Path to FESOM high-res logfile from timeseries_dpt.py: ")
    mip_dpt_calc(roms_log, fesom_log_low, fesom_log_high)
