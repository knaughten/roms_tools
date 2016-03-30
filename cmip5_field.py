from numpy import *
from netCDF4 import Dataset
from os import listdir
from datetime import date
from scipy.interpolate import interp1d

# Read CMIP5 output for the given model, experiment, and variable name.
# Return the data over the Southern Ocean (for atmosphere variables) 
# or at the northern boundary of ROMS (for ocean variables).
# Input:
# model = Model object (class definition in cmip5_paths.py)
# expt = string containing name of experiment, eg 'historical'
# var_name = string containing name of variable, eg 'tas'
# start_year, end_year = integers containing years to average over, from the
#                        beginning of start_year to the end of end_year. 
#                        Therefore, if start_year = end_year, this script will
#                        average over one year of model output.
# Output:
# data_trimmed = 3D array of model output, with dimension time x latitude x
#                longitude (for atmosphere variables, at the surface) or
#                time x depth x longitude (for ocean variables, interpolated
#                to the northern boundary of ROMS)
# axis = 1D array containing latitude values (for atmosphere variables) or
#        depth values (for ocean variables).
def cmip5_field (model, expt, var_name, start_year, end_year):

    # Northern boundary of ROMS
    nbdry = -38

    # Figure out whether it is an atmosphere or ocean variable
    if var_name in ['ps', 'tas', 'huss', 'clt', 'uas', 'vas', 'pr', 'prsn', 'evspsbl', 'rsds', 'rlds']:
        realm = 'atmos'
    elif var_name in ['thetao', 'so', 'uo', 'vo', 'zos']:
        realm = 'ocean'
    else:
        print 'Unknown variable'
        # Exit early
        return None, None

    # There is something weird going on with evaporation in the Norwegian GCMs;
    # they claim to have units of kg/m^2/s but there's no way that's correct.
    # Exclude them for now, until we figure out what the units actually are.
    if var_name == 'evspsbl' and model.name in ['NorESM1-M', 'NorESM1-ME']:
        print 'Skipping ' + model.name + ' because evaporation units are screwy'
        # Exit early
        return None, None

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

    # 1D array of time values; initialise as None and then add to it with
    # each file
    time = None
    # Similarly, a 3D array of data values (time x lat x lon for atmosphere
    # variables, time x depth x lon for ocean variables interpolated to nbdry)
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
                if realm == 'atmos':
                    # The data is already 3D (surface variable) so this is easy
                    data_tmp = id.variables[var_name][:,0:j_max+1,:]
                    # Initialise the master data array if it doesn't exist yet,
                    # otherwise add the new data values to the end
                    if data is None:
                        data = data_tmp[:,:,:]
                    else:
                        data = ma.concatenate((data, data_tmp), axis=0)
                    # Save latitude as the axis this function will return
                    axis = lat

                elif realm == 'ocean':
                    # Read data immediately south and north of nbdry
                    data_tmp1 = id.variables[var_name][:,:,j_max-1,:]
                    data_tmp2 = id.variables[var_name][:,:,j_max,:]
                    # Some of the CMIP5 ocean models are not saved as masked
                    # arrays, but rather as regular arrays with the value 0
                    # at all land points. Catch these with a try-except block
                    # and convert to masked arrays.
                    try:
                        mask = data_tmp1.mask
                        mask = data_tmp2.mask
                    except (AttributeError):
                        data_tmp1 = ma.masked_where(data_tmp1==0, data_tmp1)
                        data_tmp2 = ma.masked_where(data_tmp2==0, data_tmp2)
                                        
                    if model.name == 'inmcm4':
                        # inmcm4 has sigma coordinates
                        # Calculate z-coordinates at points on either side of
                        # northern boundary
                        sigma = id.variables['lev'][:]
                        h1 = id.variables['depth'][j_max-1,:]
                        h2 = id.variables['depth'][j_max,:]
                        depth1 = ma.empty([size(sigma), size(h1)])
                        depth2 = ma.empty([size(sigma), size(h1)])
                        for k in range(size(sigma)):
                            depth1[k,:] = -sigma[k]*h1[:]
                            depth2[k,:] = -sigma[k]*h2[:]

                        # Interpolate data to a regular z-grid
                        axis = linspace(0, max(amax(depth1),amax(depth2)), 50)
                        # Replace mask with NaNs
                        data_unmask1 = data_tmp1.data
                        data_unmask1[data_tmp1.mask] = NaN
                        data_unmask2 = data_tmp2.data
                        data_unmask2[data_tmp2.mask] = NaN
                        data_interp1 = empty([size(data_unmask1,0), size(axis), size(data_unmask1,2)])
                        data_interp2 = empty([size(data_unmask1,0), size(axis), size(data_unmask1,2)])
                        # Loop over longitude and time
                        for i in range(size(data_unmask1,2)):
                            for t in range(size(data_unmask1,0)):
                                # Fill out-of-bounds values (eg deeper than
                                # the seafloor) with NaNs
                                interp_function1 = interp1d(depth1[:,i], data_unmask1[t,:,i], bounds_error=False, fill_value=NaN)
                                data_interp1[t,:,i] = interp_function1(axis)
                                interp_function2 = interp1d(depth2[:,i], data_unmask2[t,:,i], bounds_error=False, fill_value=NaN)
                                data_interp2[t,:,i] = interp_function2(axis)
                        # Replace NaNs with mask
                        data_tmp1 = ma.masked_where(isnan(data_interp1), data_interp1)
                        data_tmp2 = ma.masked_where(isnan(data_interp2), data_interp2)

                    else:
                        # This is a z-coordinate model
                        # Save depth as the axis this function will return
                        axis = id.variables['lev'][:]

                    # Linearly interpolate to nbdry
                    data_tmp = (data_tmp2 - data_tmp1)/(lat[j_max]-lat[j_max-1])*(nbdry-lat[j_max-1]) + data_tmp1
                    # Initialise or add to data array as before
                    if data is None:
                        data = data_tmp[:,:,:]
                    else:
                        data = ma.concatenate((data, data_tmp), axis=0)

                id.close()

    # Check if we actually read any data
    if time is None or data is None:
        print 'No files found in specified date range'
        # Exit early
        return None, None

    # Figure out how many time indices are between the dates we're interested in
    num_time = count_nonzero((time >= (time_start-time_ref).days)*(time <= (time_end-time_ref).days))

    # Sort the data chronologically
    # We won't necessarily keep all of the data, so first just find the indices
    # of the sorted time array, eg [1 7 3 6] would have sorted indices [0 2 3 1]
    sort_index = argsort(time)
    # Set up data array with the correct number of time indices
    data_trimmed = ma.empty([num_time, size(data,1), size(data,2)])

    # Initialise next available time index in data_trimmed
    posn = 0
    # Loop over each time index
    for index in sort_index:
        # Figure out if it falls between the dates we're interested in
        if time[index] >= (time_start-time_ref).days and time[index] <= (time_end-time_ref).days:
            # Save model output at this time index to the new array
            data_trimmed[posn,:,:] = data[index,:,:]
            posn += 1

    return data_trimmed, axis
    

    
