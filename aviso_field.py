from numpy import *
from netCDF4 import Dataset

# Read AVISO sea surface height data between the given start and end years.
# Interpolate to the northern boundary of the circumpolar ROMS domain.
# Input:
# start_year, end_year = integers containing years to average over, from the
#                        beginning of start_year to the end of end_year.
#                        Therefore, if start_year = end_year, this script will
#                        average over one year of model output.
# Output:
# aviso_data = 2D array of AVISO data interpolated to the northern boundary,
#              with dimension time x longitude
def aviso_field (start_year, end_year):

    # Latitude of the northern boundary of the circumpolar ROMS domain
    nbdry = -38
    # Beginning of AVISO file names
    aviso_head = '/short/m68/kaa561/ROMS-CICE-MCT/data/AVISO/dt_global_allsat_msla_h_y'

    # Read AVISO latitude and longitude from the first file
    id = Dataset(aviso_head + str(start_year) + '_m01.nc', 'r')
    aviso_lat = id.variables['lat'][:]
    aviso_lon = id.variables['lon'][:]
    id.close()
    # Find the first index north of nbdry, and subtract 1 to find the last
    # index south of nbdry    
    j_max = nonzero(aviso_lat > nbdry)[0][0]
    j_min = j_max - 1
    # Only save the AVISO latitude values at these indices
    aviso_lat = aviso_lat[j_min:j_max+1]

    # Create empty array of dimension time x longitude
    aviso_data = ma.empty([12*(end_year-start_year+1), size(aviso_lon)])
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
            aviso_file = aviso_head + str(year) + '_m' + month_str + '.nc'
            # Read data
            id = Dataset(aviso_file, 'r')
            data = id.variables['sla'][0,j_min:j_max+1,:]
            id.close()

            # Linearly interpolate to nbdry and save to master array
            aviso_data[posn,:] = (data[1,:]-data[0,:])/(aviso_lat[1]-aviso_lat[0])*(nbdry-aviso_lat[0]) + data[0,:]
            posn +=1

    return aviso_data
