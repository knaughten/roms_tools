from numpy import *
from netCDF4 import Dataset
from os import listdir
from datetime import date

# Read CMIP5 output for the given model, experiment, and variable name,
# interpolated to the northern boundary of ROMS.
# This script is for ocean surface variables.
# Input:
# model = Model object (class definition in cmip5_paths.py)
# expt = string containing name of experiment, eg 'historical'
# var_name = string containing name of variable, eg 'zos'
# start_year, end_year = integers containing years to average over, from the
#                        beginning of start_year to the end of end_year. 
#                        Therefore, if start_year = end_year, this script will
#                        average over one year of model output.
# Output:
# data_trimmed = 2D array of model output, with dimension time x longitude
#                (interpolated to the northern boundary of ROMS)
def cmip5_field_ss (model, expt, var_name, start_year, end_year):

    # Northern boundary of ROMS
    nbdry = -38

    # Build Date objects for 1 Jan on start_year and 31 Dec on end_year
    time_start = date(start_year, 1, 1)
    time_end = date(end_year, 12, 31)

    # Get the directory where the model output is stored
    path = model.get_directory(expt, var_name)
    # If the string is empty, this output doesn't exist
    if path == '':
        print 'Warning: no data found for model ' + model.name + ', experiment ' + expt + ', variable ' + var_name
        # Exit early
        return None

    # 1D array of time values; initialise as None and then add to it with
    # each file
    time = None
    # Similarly, a 2D array of data values (time x lon, interpolated to nbdry)
    data = None
    # Date object containing the reference time, i.e. the date when t=0 in
    # model output
    time_ref = None

    # Loop over all files in this directory
    for file in listdir(path):

        # Check every netCDF file
        if file.endswith('.nc'):

            # Read the time values
            id = Dataset(path + file, 'r')
            time_tmp = id.variables['time'][:]
            # Parse the units to figure out the reference date
            # Units will always be in the form 'days since year-month-day' or
            # 'days since year-month-day 00:00:00' or 'days since year-day'
            # First split along the hyphens
            time_units = (id.variables['time'].units).split('-')
            # Split the first segment ('days since year') along the spaces,
            # and select the last sub-segment ('year')
            year_ref = int(time_units[0].split()[-1])
            if year_ref == 0:
                # datetime can't handle year 0
                # Set to year 1 and adjust time data accordingly
                year_ref = 1
                time_tmp = time_tmp - 365
            if len(time_units) == 3:
                # The middle segment is just 'month'
                mon_ref = int(time_units[1])
                # Split the last segment ('day' or 'day 00:00:00') along the
                # spaces, and select the first sub-segment ('day')
                day_ref = int(time_units[2].split()[0])
            elif len(time_units) == 2:
                # No month, just year and day
                # This only happens for day=1 so set month=1 also
                mon_ref = 1
                day_ref = int(time_units[1].split()[0])

            if time_ref is None:
                # This is the first file we've read
                # Initialise time_ref
                time_ref = date(year_ref, mon_ref, day_ref)
            else:
                # Compare the new reference time to the old reference time
                # and convert the time values to be with respect to the
                # existing time_ref
                new_time_ref = date(year_ref, mon_ref, day_ref)
                # Just add on the days between the old reference time and the
                # new reference time
                time_tmp = time_tmp + (new_time_ref - time_ref).days

            # Check if the time values in this file actually contain any
            # dates we're interested in
            if time_tmp[0] > (time_end-time_ref).days or time_tmp[-1] < (time_start-time_ref).days:
                # No overlap, skip this file and go to the next one
                id.close()
            else:
                # Initialise master time array if it doesn't exist yet,
                # otherwise add the new time values to the end.
                if time is None:
                    time = time_tmp[:]
                else:
                    time = concatenate((time, time_tmp), axis=0)

                # Read the latitude array
                if len(id.variables['lat'].shape) == 2:
                    # Latitude is 2D; average over longitude
                    lat = mean(id.variables['lat'][:,:], axis=1)
                else:
                    lat = id.variables['lat'][:]
                if lat[0] > lat[1]:
                    # Latitude decreases
                    # Find the first index south of nbdry
                    j_max = nonzero(lat <= nbdry)[0][0]
                elif lat[0] < lat[1]:
                    # Latitude increases
                    # Find the first index where latitude exceeds nbdry
                    j_max = nonzero(lat >= nbdry)[0][0]                
                # Trim latitude values
                lat = lat[0:j_max+1]

                # Read model output
                # Linearly interpolate to nbdry
                data_tmp1 = id.variables[var_name][:,j_max-1,:]
                data_tmp2 = id.variables[var_name][:,j_max,:]
                data_tmp = (data_tmp2 - data_tmp1)/(lat[j_max]-lat[j_max-1])*(nbdry-lat[j_max-1]) + data_tmp1
                # Initialise or add to data array as before
                if data is None:
                    data = data_tmp[:,:]
                else:
                    data = ma.concatenate((data, data_tmp), axis=0)

                id.close()

    # Check if we actually read any data
    if time is None or data is None:
        print 'No files found in specified date range'
        # Exit early
        return None

    # Figure out how many time indices are between the dates we're interested in
    num_time = count_nonzero((time >= (time_start-time_ref).days)*(time <= (time_end-time_ref).days))

    # Sort the data chronologically
    # We won't necessarily keep all of the data, so first just find the indices
    # of the sorted time array, eg [1 7 3 6] would have sorted indices [0 2 3 1]
    sort_index = argsort(time)
    # Set up data array with the correct number of time indices
    data_trimmed = ma.empty([num_time, size(data,1)])

    # Initialise next available time index in data_trimmed
    posn = 0
    # Loop over each time index
    for index in sort_index:
        # Figure out if it falls between the dates we're interested in
        if time[index] >= (time_start-time_ref).days and time[index] <= (time_end-time_ref).days:
            # Save model output at this time index to the new array
            data_trimmed[posn,:] = data[index,:]
            posn += 1

    return data_trimmed
    

    
