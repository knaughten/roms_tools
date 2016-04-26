from numpy import *
from netCDF4 import Dataset, num2date, date2num
from os import listdir
from scipy.interpolate import interp1d

# Read CMIP5 output for the given model, experiment, and variable name. Return
# the monthly climatology as well as the grid.
# Input:
# model = Model object (class definition in cmip5_paths.py)
# expt = string containing name of experiment, eg 'historical'
# var_name = string containing name of variable, eg 'tas'
# start_year, end_year = integers containing years over which to calculate
#                        monthly climatology, from the beginning of start_year
#                        to the end of end_year. Therefore, if 
#                        start_year = end_year, this script will read one year
#                        of output with no climatological averaging.
# Output:
# data = array of model output, with dimension month x latitude x longitude
#        (for atmosphere variables, at the surface) or month x depth x latitude
#        x longitude (for ocean variables), possibly with units converted to be
#        more comparable to other models and/or reanalyses
# lon, lat = longitude and latitude arrays (1D or 2D)
# depth = depth array (1D or 3D) if ocean variable, otherwise None
def cmip5_field (model, expt, var_name, start_year, end_year):

    # Conversion from K to C
    degKtoC = -273.15

    # Figure out whether it is an atmosphere or ocean variable
    if var_name in ['ps', 'tas', 'huss', 'clt', 'uas', 'vas', 'pr', 'prsn', 'evspsbl', 'rsds', 'rlds']:
        realm = 'atmos'
    elif var_name in ['thetao', 'so', 'uo', 'vo', 'zos']:
        realm = 'ocean'
    else:
        print 'Unknown variable'
        # Exit early
        return None, None, None, None

    # There is something weird going on with evaporation in the Norwegian GCMs;
    # they claim to have units of kg/m^2/s but there's no way that's correct.
    # Exclude them for now, until we figure out what the units actually are.
    if var_name == 'evspsbl' and model.name in ['NorESM1-M', 'NorESM1-ME']:
        print 'Skipping ' + model.name + ' because evaporation units are screwy'
        # Exit early
        return None, None, None, None

    # Get the directory where the model output is stored
    path = model.get_directory(expt, var_name)
    # If the string is empty, this output doesn't exist
    if path == '':
        print 'Warning: no data found for model ' + model.name + ', experiment ' + expt + ', variable ' + var_name
        # Exit early
        return None, None, None, None

    data = None
    # Number of records for each month
    num_months = zeros(12)

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
            curr_units = time_id.units
            if curr_units == 'days since 0001-01':
                curr_units = 'days since 0001-01-01'
            if curr_units == 'days since 0000-01-01 00:00:00':
                curr_units = 'days since 0001-01-01 00:00:00'
            time = num2date(time_id[:], units=curr_units, calendar=time_id.calendar)

            # Check if the time values in this file actually contain any
            # dates we're interested in
            if time[0].year > end_year or time[-1].year < start_year:
                # No overlap, skip this file and go to the next one
                id.close()
            else:
                # Read the grid
                lat = id.variables['lat'][:]
                lon = id.variables['lon'][:]
                if realm == 'ocean':
                    if model.name == 'inmcm4':
                        # inmcm4 has sigma coordinates; convert to z
                        sigma = id.variables['lev'][:]
                        h = id.variables['depth'][:,:]
                        depth = ma.empty([size(sigma), size(h,0), size(h,1)])
                        for k in range(size(sigma)):
                            depth[k,:,:] = -sigma[k]*h
                    else:
                        depth = id.variables['lev'][:]
                else:
                    depth = None                

                # Read model output, one timestep at a time to save memory
                for t in range(size(time)):
                    if realm == 'atmos':
                        data_tmp = id.variables[var_name][t,:,:]
                    elif realm == 'ocean':
                        data_tmp = id.variables[var_name][t,:,:,:]
                    # Some of the CMIP5 ocean models are not saved as masked
                    # arrays, but rather as regular arrays with the value 0
                    # at all land points. Catch these with a try-except block
                    # and convert to masked arrays.
                    try:
                        mask = data_tmp.mask
                    except (AttributeError):
                        data_tmp = ma.masked_where(data_tmp==0, data_tmp)

                    if data is None:
                        # Set up data array of correct size
                        if realm == 'atmos':
                            data = ma.empty([12, size(data_tmp,0), size(data_tmp,1)])
                            for tt in range(12):
                                # Initialise as zeros masked with land mask
                                data[tt,:,:] = data_tmp[:,:]*0.0
                        elif realm == 'ocean':
                            data = ma.empty([12, size(data_tmp,0), size(data_tmp,1), size(data_tmp,2)])
                            for tt in range(12):
                                data[tt,:,:,:] = data_tmp[:,:,:]*0.0

                    # Figure out if it falls between the dates we want
                    if time[t].year >= start_year and time[t].year <= end_year:
                        # Find the month index, increment curr_month, and add
                        # to data array
                        curr_month = time[t].month-1
                        num_months[curr_month] += 1
                        if realm == 'atmos':
                            data[curr_month,:,:] += data_tmp[:,:]
                        elif realm == 'ocean':
                            data[curr_month,:,:,:] += data_tmp[:,:,:]

    # Check if we actually read any data
    if data is None:
        print 'No files found in specified date range'
        # Exit early
        return None, None, None, None

    # Convert from monthly sums to monthly averages
    for t in range(12):
        if num_months[t] == 0:
            # None for this month; mask it
            if realm == 'atmos':
                data[t,:,:] = ma.masked
            elif realm == 'ocean':
                data[t,:,:,:] = ma.masked
        else:
            if realm == 'atmos':
                data[t,:,:] /= num_months[t]
            elif realm == 'ocean':
                data[t,:,:,:] /= num_months[t]

    # Conversions if necessary
    if var_name in ['pr', 'prsn', 'evspsbl']:
        # Convert precip/snowfall/evap from kg/m^2/s to
        # 1e-6 kg/m^2/s
        data = 1e6*data
    elif var_name == 'ps':
        # Convert surface pressure from Pa to kPa
        data = 1e-3*data
    elif var_name == 'tas':
        # Convert temperature from K to C
        data = data + degKtoC
    elif var_name == 'thetao' and amin(data) > 100:
        # Convert ocean temperature from K to C if needed
        data = data + degKtoC
    elif var_name == 'so' and amax(data) < 1:
        # Convert salinity from fraction to psu if needed
        data = 1e3*data

    return data, lon, lat, depth
    

    
