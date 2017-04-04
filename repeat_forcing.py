from netCDF4 import Dataset
from numpy import *

# Convert a 1-year dataset into an annually repeating dataset for ROMS-CICE.
# The easiest way to do this is by making 4 identical files, and then for
# the year in that 4 which is a leap year, make Feb 29th the average of
# Feb 28th and March 1st. Set the cycle length attribute to 1461 days (4 years
# where one is a leap year).
# NB: This script assumes the first file has been copied three times to the
#     correct names of the next three files.
# Sort of NB: This script uses the basic definition of leap year = year
#             divisible by 4. This does not hold for all years ending in 00,
#             e.g. 2000 was a leap year but 1900 wasn't. If you have sub-daily
#             data for 1900 I am impressed.
# Input:
# directory = string containing path to the directory where the files exist
# head = string containing the first part of each filename
# tail = string containing the last part of each filename
#        NB: This script assumes the files differ only by the year, e.g.
#            AN_1995_unlim.nc through AN_1998_unlim.nc
# year_start = integer containing the first year, which is also the year of
#              the original data
# perday = number of records per day (for ERA-Interim, either 2 or 4)
# var_list = list of strings containing variable names for which Feb 29th
#            will need to be interpolated (i.e. every NetCDF variable which
#            depends on time, except for time itself).

def process (directory, head, tail, year_start, perday, var_list):

    # Loop through the four years
    for year in range(year_start,year_start+4):

        file = directory + head + str(year) + tail
        print 'Processing ' + file
        id = Dataset(file, 'a')
        # Set the cycle_length to 4 years
        id.variables['time'].cycle_length = 365.0*4 + 1

        if year > year_start:
            # For every year except the first year, alter the time axis
            # to add on the correct number of days
            days_to_add = 365.0*(year-year_start)            
            for year_tmp in range(year_start, year):
                if year_tmp % 4 == 0:
                    # A leap year has occurred since year_start
                    days_to_add = days_to_add+1
            id.variables['time'][:] = id.variables['time'][:] + days_to_add

            if year % 4 == 0:                
                print 'This is a leap year'
                # Add Feb 29th to the time axis
                num_time = id.variables['time'].shape[0]
                for i in range(perday):
                    id.variables['time'][num_time+i] = id.variables['time'][num_time+i-perday] + 1.0
                # Save indices of Feb 29th and March 1st in the new time axis
                feb29 = (31+28)*perday
                mar1 = feb29 + perday

                for var in var_list:
                    # Interpolate each variable on Feb 29th to be the mean of
                    # the Feb 28th and March 1st values at the same time of day
                    print 'Processing variable ' + var
                    var_feb29 = zeros((perday, id.variables[var].shape[1], id.variables[var].shape[2]))
                    for i in range(perday):
                        var_feb29[i,:,:] = 0.5*(id.variables[var][feb29-perday+i,:,:] + id.variables[var][feb29+i,:,:])
                    # Shift the existing March-December data along by perday
                    # indices
                    id.variables[var][mar1:,:,:] = id.variables[var][feb29:num_time,:,:]
                    # Insert interpolated Feb 29th data into this space
                    id.variables[var][feb29:mar1,:,:] = var_feb29
    id.close()


# Command-line interface
if __name__ == '__main__':

    # User parameters to edit here

    # Data where files <an_head>yyyy<tail> and <fc_head>yyyy<tail> exist
    directory = '/short/m68/kaa561/metroms_iceshelf/data/ERA_Interim/subdaily/30day_smoothed/'
    # First part of filename for AN and FC files
    an_head = 'AN_'
    fc_head = 'FC_'
    # Last part of filename (common to both AN and FC)
    tail = '_subdaily.nc'
    # Year to build data from
    year_start = 1995
    # Number of records per day
    an_perday = 4
    fc_perday = 2
    # Variables to interpolate Feb 29th for each file
    an_var = ['Pair', 'Tair', 'Qair', 'cloud', 'Uwind', 'Vwind']
    fc_var = ['rain', 'snow']

    if year_start % 4 == 0:
        # This script assumes year_start has 365 days and then the leap year
        # 1, 2, or 3 years after that will have Feb 29th added in by interpolation.
        # However you could rework this script to remove Feb 29th from year_start
        # and every following year except the leap year.
        print 'year_start cannot be a leap year. Either choose a different year_start or rework this script.'
        exit

    # Run the actual script
    process(directory, an_head, tail, year_start, an_perday, an_var)
    process(directory, fc_head, tail, year_start, fc_perday, fc_var)
                
                
                
                
    
