from netCDF4 import Dataset, num2date
from numpy import *

def calc_ice_prod (in_file, start_year, end_year, out_file):

    # cm/day to m/day conversion
    cm_to_m = 1e-2

    # Read the grid
    id = Dataset(in_file, 'r')
    lon = id.variables['TLON'][:,:]
    lat = id.variables['TLAT'][:,:]
    num_lon = size(lon,1)
    num_lat = size(lat,0)
    # Set up array to hold output data
    ice_prod = ma.empty([num_lat, num_lon])
    ice_prod[:,:] = 0.0

    # Read the time values
    time_id = id.variables['time']
    cice_time = num2date(time_id[:], units=time_id.units, calendar='standard')
    num_time = time_id.size
    # Find the first timestep and how many days in it we care about
    if start_year == 1992:
        # Starting at the beginning of the simulation
        start_t = 0
        start_days = 5
    else:
        for t in range(num_time):
            if cice_time[t].year == start_year and cice_time[t].month-1 == 0 and cice_time[t].day in range(2,6+1):
                start_t = t
                start_days = cice_time[t].day-1
                break
    # Find the last timestep and how many days in it we care about
    if end_year == 2016:
        # Ending at the end of the simulation
        end_t = num_time-1
        end_days = 5
    else:
        for t in range(num_time):
            if cice_time[t].year == end_year+1 and cice_time[t].month-1 == 0 and cice_time[t].day in range(1,5+1):
                end_t = t
                end_days = 6-cice_time[t].day
                break
    print 'Starting at index ' + str(start_t) + ' (' + str(cice_time[start_t].year) + '-' + str(cice_time[start_t].month) + '-' + str(cice_time[start_t].day) + '), ' + str(start_days) + ' days'
    print 'Ending at index ' + str(end_t) + ' (' + str(cice_time[end_t].year) + '-' + str(cice_time[end_t].month) + '-' + str(cice_time[end_t].day) + '), ' + str(end_days) + ' days'

    # Integrate the start days
    thdgr_start = id.variables['frazil'][start_t,:,:] + id.variables['congel'][start_t,:,:] - id.variables['meltt'][start_t,:,:] - id.variables['meltb'][start_t,:,:] - id.variables['meltl'][start_t,:,:]
    index = thdgr_start < 0
    thdgr_start[index] = 0
    ice_prod += thdgr_start*cm_to_m*start_days
    # Integrate the middle days
    thdgr_mid = id.variables['frazil'][start_t+1:end_t,:,:] + id.variables['congel'][start_t+1:end_t,:,:] - id.variables['meltt'][start_t+1:end_t,:,:] - id.variables['meltb'][start_t+1:end_t,:,:] - id.variables['meltl'][start_t+1:end_t,:,:]
    index = thdgr_mid < 0
    thdgr_mid[index] = 0
    ice_prod += sum(thdgr_mid, axis=0)*cm_to_m*5
    # Integrate the end days
    thdgr_end = id.variables['frazil'][end_t,:,:] + id.variables['congel'][end_t,:,:] - id.variables['meltt'][end_t,:,:] - id.variables['meltb'][end_t,:,:] - id.variables['meltl'][end_t,:,:]
    index = thdgr_end < 0
    thdgr_end[index] = 0
    ice_prod += thdgr_end*cm_to_m*end_days
    # Get annual average in m/y
    ice_prod = ice_prod/(end_year-start_year+1)

    # Write to file
    id = Dataset(out_file, 'w')
    id.createDimension('ni', size(lon,1))
    id.createDimension('nj', size(lon,0))
    id.createDimension('time', 4)
    id.createVariable('TLON', 'f8', ('nj', 'ni'))
    id.variables['TLON'][:,:] = lon
    id.createVariable('TLAT', 'f8', ('nj', 'ni'))
    id.variables['TLAT'][:,:] = lat
    id.createVariable('ice_prod', 'f8', ('nj', 'ni'))
    id.variables['ice_prod'].units = 'm/y'
    id.variables['ice_prod'][:,:] = ice_prod
    id.close()


# Command-line interface
if __name__ == "__main__":

    in_file = raw_input("Path to CICE iceh_tot.nc file for entire simulation: ")
    start_year = int(raw_input("First year to process: "))
    end_year = int(raw_input("End year to process: "))
    out_file = raw_input("Path to desired output file: ")
    calc_ice_prod (in_file, start_year, end_year, out_file)
                
            
        
            
    
