from netCDF4 import Dataset, num2date
from numpy import *

# Calculate the seasonal climatology (DJF, MAM, JJA, SON) of ocean temperature
# and salinity during a ROMS simulation and save to a NetCDF file.
# Input:
# directory = path to ROMS output directory containing ocean averages files,
#             assuming 5-day averages
# start_index, end_index = integers containing range of files to process. For
#                          example, start_index=1 and end_index=102 will
#                          process files ocean_avg_0001.nc through
#                          ocean_avg_0102.nc.
# out_file = path to desired output file
# start_year = optional integer containing the first year to consider
def seasonal_climatology_roms (directory, start_index, end_index, out_file, start_year=1992):

    # Starting and ending months (1-based) for each season
    start_month = [12, 3, 6, 9]
    end_month = [2, 5, 8, 11]
    # Starting and ending days of the month (1-based) for each season
    # Assume no leap years, we'll fix this later if we need
    start_day = [1, 1, 1, 1]
    end_day = [28, 31, 31, 30]

    # Read grid from the first file
    id = Dataset(directory + index_to_file(start_index), 'r')
    lon = id.variables['lon_rho'][:,:]
    lat = id.variables['lat_rho'][:,:]
    # Read number of vertical levels
    num_depth = id.variables['temp'].shape[1]
    id.close()

    # Set up arrays to integrate seasonal climatology of temp and salt
    seasonal_temp = ma.empty([4, num_depth, size(lon,0), size(lon,1)])
    seasonal_salt = ma.empty([4, num_depth, size(lon,0), size(lon,1)])
    seasonal_temp[:,:,:,:] = 0.0
    seasonal_salt[:,:,:,:] = 0.0
    # Also integrate number of days in each season
    ndays = zeros(4)

    # Loop over files
    for index in range(start_index, end_index+1):
        filename = directory + index_to_file(index)
        print 'Processing ' + filename
        id = Dataset(filename, 'r')
        print 'Reading data'
        # Read temperature and salinity
        temp = id.variables['temp'][:,:,:,:]
        salt = id.variables['salt'][:,:,:,:]
        # Read the time values and convert to Date objects
        time_id = id.variables['ocean_time']
        time = num2date(time_id[:], units=time_id.units, calendar=time_id.calendar.lower())
        id.close()
        print 'Integrating climatology'
        # Loop over timesteps
        for t in range(size(time)):
            # Make sure we are past start_year
            if time[t].year >= start_year:
                print '...time index ' + str(t+1) + ' of ' + str(size(time))
                # 5-day averages marked with middle day's date
                year = time[t].year
                month = time[t].month
                day = time[t].day
                # Get the season of this middle day
                if month in [12, 1, 2]:
                    # DJF
                    season = 0
                elif month in [3, 4, 5]:
                    # MAM
                    season = 1
                elif month in [6, 7, 8]:
                    # JJA
                    season = 2
                elif month in [9, 10, 11]:
                    # SON
                    season = 3
                # Check for leap years
                leap_year = False
                if mod(year, 4) == 0:
                    leap_year = True
                    if mod(year, 100) == 0:
                        leap_year = False
                        if mod(year, 400) == 0:
                            leap_year = True
                # Update last day in February
                if leap_year:
                    end_day[0] = 29
                else:
                    end_day[0] = 28
                if month == start_month[season]:
                    # We are in the first month of the season
                    if day-2 < start_day[season]:
                        # Partially spills over into the previous season
                        prev_season = mod(season-1, 4)
                        # How many days does it spill over by?
                        spill_days = start_day[season]-day+2
                        # Should be either 1 or 2
                        if spill_days not in [1,2]:
                            print 'Problem: spill_days is ' + str(spill_days)
                            print 'Timestep ' + str(t+1)
                            print 'Year ' + str(year)
                            print 'Month ' + str(month+1)
                            print 'Day ' + str(day)
                            return
                        # Split between previous season and this season
                        seasonal_temp[prev_season,:,:,:] += temp[t,:,:,:]*spill_days
                        seasonal_salt[prev_season,:,:,:] += salt[t,:,:,:]*spill_days
                        ndays[prev_season] += spill_days
                        seasonal_temp[season,:,:,:] += temp[t,:,:,:]*(5-spill_days)
                        seasonal_salt[season,:,:,:] += salt[t,:,:,:]*(5-spill_days)
                        ndays[season] += 5-spill_days
                    else:
                        # Entirely within the season
                        seasonal_temp[season,:,:,:] += temp[t,:,:,:]*5
                        seasonal_salt[season,:,:,:] += salt[t,:,:,:]*5
                        ndays[season] += 5
                elif month == end_month[season]:
                    # We are in the last month of the season
                    if day+2 > end_day[season]:
                        # Partially spills over into the next season
                        next_season = mod(season+1, 4)
                        # How many days does it spill over by?
                        spill_days = day+2-end_day[season]
                        # Should be either 1 or 2
                        if spill_days not in [1,2]:
                            print 'Problem: spill_days is ' + str(spill_days)
                            print 'Timestep ' + str(t+1)
                            print 'Year ' + str(year)
                            print 'Month ' + str(month+1)
                            print 'Day ' + str(day)
                            return
                        # Split between this season and next season
                        seasonal_temp[next_season,:,:,:] += temp[t,:,:,:]*spill_days
                        seasonal_salt[next_season,:,:,:] += salt[t,:,:,:]*spill_days
                        ndays[next_season] += spill_days
                        seasonal_temp[season,:,:,:] += temp[t,:,:,:]*(5-spill_days)
                        seasonal_salt[season,:,:,:] += salt[t,:,:,:]*(5-spill_days)
                        ndays[season] += 5-spill_days
                    else:
                        # Entirely within the season
                        seasonal_temp[season,:,:,:] += temp[t,:,:,:]*5
                        seasonal_salt[season,:,:,:] += salt[t,:,:,:]*5
                        ndays[season] += 5
                else:
                    # We are in the middle month of the season
                    # The 5 days in this index are entirely within the season
                    seasonal_temp[season,:,:,:] += temp[t,:,:,:]*5
                    seasonal_salt[season,:,:,:] += salt[t,:,:,:]*5
                    ndays[season] += 5            
    # Convert from sums to averages
    for season in range(4):
        seasonal_temp[season,:,:,:] = seasonal_temp[season,:,:,:]/ndays[season]
        seasonal_salt[season,:,:,:] = seasonal_salt[season,:,:,:]/ndays[season]

    # Write to file
    print 'Writing ' + out_file
    id = Dataset(out_file, 'w')
    id.createDimension('xi_rho', size(lon,1))
    id.createDimension('eta_rho', size(lon,0))
    id.createDimension('s_rho', num_depth)    
    id.createDimension('time', 4)
    id.createVariable('lon_rho', 'f8', ('eta_rho', 'xi_rho'))
    id.variables['lon_rho'].long_name = 'longitude of rho-points'
    id.variables['lon_rho'].units = 'degree_east'
    id.variables['lon_rho'][:,:] = lon
    id.createVariable('lat_rho', 'f8', ('eta_rho', 'xi_rho'))
    id.variables['lat_rho'].long_name = 'latitude of rho-points'
    id.variables['lat_rho'].units = 'degree_north'
    id.variables['lat_rho'][:,:] = lat
    id.createVariable('time', 'f8', ('time'))
    id.variables['time'].units = 'season'
    id.variables['time'].description = 'DJF, MAM, JJA, SON'
    id.variables['time'][:] = arange(1,4+1)
    id.createVariable('temp', 'f8', ('time', 's_rho', 'eta_rho', 'xi_rho'))
    id.variables['temp'].units = 'degC'
    id.variables['temp'][:,:,:,:] = seasonal_temp
    id.createVariable('salt', 'f8', ('time', 's_rho', 'eta_rho', 'xi_rho'))
    id.variables['salt'].units = 'psu'
    id.variables['salt'][:,:,:,:] = seasonal_salt
    id.close()
        

# Given an integer, return the filename for the corresponding ocean averages
# file. For example, index_to_file(1) = 'ocean_avg_0001.nc', and
# index_to_file(95) = 'ocean_avg_0095.nc'.
def index_to_file (index):

    if index < 10:
        return 'ocean_avg_000' + str(index) + '.nc'
    elif index < 100:
        return 'ocean_avg_00' + str(index) + '.nc'
    elif index < 1000:
        return 'ocean_avg_0' + str(index) + '.nc'
    else:
        return 'ocean_avg_' + str(index) + '.nc'


# Command-line interface
if __name__ == "__main__":

    directory = raw_input("Path to ROMS output directory: ")
    start_index = int(raw_input("Index of first ocean averages file (e.g. ocean_avg_0001.nc has index 1): "))
    end_index = int(raw_input("Index of last ocean averages file: "))
    out_file = raw_input("Path to desired output file: ")
    seasonal_climatology_roms(directory, start_index, end_index, out_file)
