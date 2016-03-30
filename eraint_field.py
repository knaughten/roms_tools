from numpy import *
from netCDF4 import Dataset

# Read ERA-Interim monthly data for the given variable name, between the given
# start and end years. Return the data over the Southern Ocean.
# Input:
# var_name = string containing name of variable, eg 't2m'
# start_year, end_year = integers containing years to average over, from the
#                        beginning of start_year to the end of end_year.
#                        Therefore, if start_year = end_year, this script will
#                        average over one year of model output.
# Output:
# era_data = 3D array of ERA-Interim data, with dimension time x latitude x
#            longitude, possibly with units converted to be more comparable
#            to CMIP5 models
# era_lat = 1D array containing latitude values.
def eraint_field (var_name, start_year, end_year):

    # Latitude of the northern boundary of the circumpolar ROMS domain
    nbdry = -38
    # Directory where ERA-Interim monthly averaged data is stored
    era_dir = '/short/y99/kaa561/FESOM/ERA_Interim_monthly/'
    # String that ERA-Interim files end with
    era_tail = '_monthly_orig.nc'

    # Latent heat of vapourisation, J/kg
    Lv = 2.5e6
    # Ideal gas constant for water vapour, J/K/kg
    Rv = 461.5
    # Density of water, kg/m^3
    rho_w = 1e3
    # Conversion from K to C
    degKtoC = -273.15

    # Read ERA-Interim latitude and longitude
    id = Dataset(era_dir + 'AN_' + str(start_year) + era_tail, 'r')
    era_lat = id.variables['lat'][:]
    era_lon = id.variables['lon'][:]
    id.close()
    # Latitude values are decreasing for some reason, i.e. 90 to -90
    # Find the first index south of nbdry, and subtract 1 to find the last
    # index north of nbdry
    j_min = nonzero(era_lat < nbdry)[0][0] - 1
    # Save only the values from j_min southward
    era_lat = era_lat[j_min:]

    # Figure out how ERA-Interim filename will start; the atmospheric
    # variables for each year are split between 3 different files
    if var_name in ['sp', 't2m', 'd2m', 'tcc', 'u10', 'v10']:
        era_head = era_dir + 'AN_'
    elif var_name in ['tp', 'sf']:
        era_head = era_dir + 'FC_'
    elif var_name in ['e', 'ssrd', 'strd']:
        era_head = era_dir + 'ER_'

    # Create empty array of dimension time x latitude x longitude
    era_data = ma.empty([12*(end_year-start_year+1), size(era_lat), size(era_lon)])
    # Initialise next available time index in this array
    posn = 0
    # Loop over years
    for year in range(start_year, end_year+1):

        # Construct filename
        era_file = era_head + str(year) + era_tail
        # Read data
        id = Dataset(era_file, 'r')
        data = id.variables[var_name][:,j_min:,:]

        # Perform conversions if necessary
        if var_name == 'sp':
            # Convert pressure from Pa to kPa
            data = 1e-3*data
        elif var_name == 't2m':
            # Convert temperature from K to C
            data = data + degKtoC
        elif var_name == 'd2m':
            # Calculate specific humidity from dewpoint temperature
            # and surface pressure
            d2m = data
            sp = id.variables['sp'][:,j_min:,:]
            # Intermediate step to calculate vapour pressure
            e = 611*exp(Lv/Rv*(1/273.0 - 1/d2m))
            data = 0.622*e/(sp - 0.378*e)
        elif var_name == 'tcc':
            # Convert total cloud cover from fraction to percent
            data = data*100
        elif var_name in ['tp', 'sf', 'e']:
            # Convert precip/snowfall/evap from kg/m^2/s to 1e-6 kg/m^2/s
            data = 1e6*data*rho_w/(12*60*60)
            if var_name == 'e':
                # ERA-Interim defines evaporation as negative; fix this
                data = -1*data
        elif var_name in ['ssrd', 'strd']:
            # Convert from J/m^2, integrated over 12 hours, to W/m^2
            data = data/(12*60*60)

        # Save each timestep to master era_data array
        for t in range(size(data, 0)):
            era_data[posn,:,:] = data[t,:,:]
            posn += 1

    return era_data, era_lat
