from netCDF4 import Dataset
from numpy import *

# Convert a 1-year monthly dataset into an annually repeating dataset for 
# ROMS-CICE. The easiest way to do this is by making 4 identical files, altering
# the time axis so data is always on the 15th of a month, and setting the cycle
# length attribute to 1461 days (4 years where one is a leap year).
# NB: This script assumes the first file has been copied three times to the
#     correct names of the next three files.
# Sort of NB: This script uses the basic definition of leap year = year
#             divisible by 4. This does not hold for all years ending in 00,
#             e.g. 2000 was a leap year but 1900 wasn't.
# Input:
# directory = string containing path to the directory where the files exist
# head = string containing the first part of each filename
# tail = string containing the last part of each filename
#        NB: This script assumes the files differ only by the year, e.g.
#            AN_1995_unlim.nc through AN_1998_unlim.nc
# year_start = integer containing the first year, which is also the year of
#              the original data

def process (directory, head, tail, year_start):

    days_per_month = array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

    # Loop through the four years
    for year in range(year_start,year_start+4):

        # Check if this is a leap year and update days in February
        if year % 4 == 0:
            days_per_month[1] = 29
        else:
            days_per_month[1] = 28

        # Make a time axis with data on the 15th of every month
        # First get the number of days between year_start and the current year
        start_day = 365*(year-year_start)
        if year > year_start:
            for year_tmp in range(year_start, year):
                if year_tmp % 4 == 0:
                    # A leap year has occurred 
                    start_day += 1
        # Start on Jan 15th at midnight
        time = [start_day + 14]
        # Loop over months
        for month in range(1,12):
            time.append(start_day + sum(days_per_month[0:month]) + 14)        

        file = directory + head + str(year) + tail
        print 'Processing ' + file
        id = Dataset(file, 'a')
        id.variables['time'][:] = time
        # Set the cycle_length to 4 years
        id.variables['time'].cycle_length = 365.0*4 + 1
        id.close()


# Command-line interface
if __name__ == '__main__':

    # User parameters to edit here

    # Data where files <an_head>yyyy<tail> and <fc_head>yyyy<tail> exist
    directory = '/short/m68/kaa561/metroms_iceshelf/data/ERA_Interim/'
    # First part of filename for AN and FC files
    an_head = 'AN_'
    fc_head = 'FC_'
    # Last part of filename (common to both AN and FC)
    tail = '_monthly.nc'
    # Year to build data from
    year_start = 1995

    # Run the actual script
    process(directory, an_head, tail, year_start)
    process(directory, fc_head, tail, year_start)
                
                
                
                
    
