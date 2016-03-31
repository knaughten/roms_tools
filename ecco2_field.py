from numpy import *
from netCDF4 import Dataset

# Read ECCO2 data for the given variable name, between the given start and end
# years. Interpolate to the northern boundary of the circumpolar ROMS domain.
# Input:
# var_name = string containing name of variable, eg 'THETA'
# start_year, end_year = integers containing years to average over, from the
#                        beginning of start_year to the end of end_year.
#                        Therefore, if start_year = end_year, this script will
#                        average over one year of model output.
# Output:
# ecco2_data = 3D array of ECCO2 data, with dimension time x depth x longitude
# ecco2_depth = 1D array containing depth values
def ecco2_field (var_name, start_year, end_year):

    # Latitude of the northern boundary of the circumpolar ROMS domain
    nbdry = -38
    # Directory where ECCO2 data is stored
    ecco2_dir = '/short/m68/kaa561/ROMS-CICE-MCT/data/ECCO2/raw/'
    # Middle of ECCO2 file names
    ecco2_mid = '.1440x720x50.'

    # Read ECCO2 latitude, longitude, and depth from the first file
    id = Dataset(ecco2_dir + 'THETA' + ecco2_mid + str(start_year) + '01.nc', 'r')
    ecco2_lat = id.variables['LATITUDE_T'][:]
    ecco2_lon = id.variables['LONGITUDE_T'][:]
    ecco2_depth = id.variables['DEPTH_T'][:]
    id.close()
    # Find the first index north of nbdry, and subtract 1 to find the last
    # index south of nbdry    
    j_max = nonzero(ecco2_lat > nbdry)[0][0]
    j_min = j_max - 1
    # Only save the ECCO2 latitude values at these indices
    ecco2_lat = ecco2_lat[j_min:j_max+1]

    # Create empty array of dimension time x depth x longitude
    ecco2_data = ma.empty([12*(end_year-start_year+1), size(ecco2_depth), size(ecco2_lon)])
    # Initialise next available time index in this array
    posn = 0
    # Loop over years and months
    for year in range(start_year, end_year+1):
        for month in range(12):
            # Construct filename
            if (month+1) < 10:
                month_str = '0' + str(month+1)
            else:
                month_str = str(month+1)
            ecco2_file = ecco2_dir + str(var_name) + ecco2_mid + str(year) + month_str + '.nc'
            # Read data
            id = Dataset(ecco2_file, 'r')
            data = id.variables[var_name][0,:,j_min:j_max+1,:]
            id.close()

            # Linearly interpolate to nbdry and save to master array
            ecco2_data[posn,:,:] = (data[:,1,:]-data[:,0,:])/(ecco2_lat[1]-ecco2_lat[0])*(nbdry-ecco2_lat[0]) + data[:,0,:]
            posn +=1

    return ecco2_data, ecco2_depth

