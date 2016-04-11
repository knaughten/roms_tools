from numpy import *
from netCDF4 import Dataset, num2date, date2num
from os import listdir

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
# month_indices = 1D array containing the month indices (1 to 12) of each time 
#                 index in data_trimmed.
def cmip5_field_ss (model, expt, var_name, start_year, end_year):

    # Northern boundary of ROMS
    nbdry = -30

    # Build Date objects for 1 Jan on start_year and 31 Dec on end_year
    time_start = date(start_year, 1, 1)
    time_end = date(end_year, 12, 31)

    # Get the directory where the model output is stored
    path = model.get_directory(expt, var_name)
    # If the string is empty, this output doesn't exist
    if path == '':
        print 'Warning: no data found for model ' + model.name + ', experiment ' + expt + ', variable ' + var_name
        # Exit early
        return None, None

    # 1D array of time values (as datetime objects); initialise as None and
    # then add to it with each file
    time = None
    # Similarly, a 2D array of data values (time x lon, interpolated to nbdry)
    data = None

    # Loop over all files in this directory
    for file in listdir(path):

        # Check every netCDF file
        if file.endswith('.nc'):

            # Read the time values
            id = Dataset(path + file, 'r')
            time_id = id.variables['time']
            if amin(time_id[:]) < 0:
                # Missing values here; this occurs for one 1900-1949 file
                # We can just skip it
                break
            # Convert to datetime objects
            time_tmp = num2date(time_id[:], units=time_id.units, calendar=time_id.calendar)

            # Check if the time values in this file actually contain any
            # dates we're interested in
            if time_tmp[0].year > end_year or time_tmp[-1].year < start_year:
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
        return None, None

    # Figure out how many time indices are between the dates we're interested in
    num_time = 0
    for t in time:
        if t.year >= start_year and t.year <= end_year:
            num_time +=1
    # Set up data array with the correct number of time indices
    data_trimmed = ma.empty([num_time, size(data,1)])
    # Also set up array of corresponding months
    month_indices = []

    # Sort the data chronologically
    # First convert the datetime objects back to floats, just for the purpose
    # of sorting
    time_floats = date2num(time, units='days since 0001-01-01 00:00:00', calendar='standard')
    # We won't necessarily keep all of the data, so first just find the indices
    # of the sorted time array, eg [1 7 3 6] would have sorted indices [0 2 3 1]
    sort_index = argsort(time_floats)

    # Initialise next available time index in data_trimmed
    posn = 0
    # Loop over each time index
    for index in sort_index:
        # Figure out if it falls between the dates we're interested in
        if time[index].year >= start_year and time[index].year <= end_year:
            # Save model output at this time index to the new array
            data_trimmed[posn,:] = data[index,:]
            # Save month index
            month_indices.append(time[index].month)
            posn += 1

    return data_trimmed, month_indices
    

    
