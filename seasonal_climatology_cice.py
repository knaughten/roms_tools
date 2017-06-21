from netCDF4 import Dataset, num2date
from numpy import *
from os import listdir

# Calculate the seasonal climatology (DJF, MAM, JJA, SON) of sea ice
# concentration and thickness during a CICE simulation and save to a NetCDF
# file.
# Input:
# directory = path to CICE output history directory containing files of the
#             form iceh.*.nc, assuming 5-day averages with one time index per
#             file. The script will process all these files.
# out_file = path to desired output file
def seasonal_climatology_cice (directory, out_file):

    # Starting and ending months (1-based) for each season
    start_month = [12, 3, 6, 9]
    end_month = [2, 5, 8, 11]
    # Starting and ending days of the month (1-based) for each season
    # Assume no leap years, we'll fix this later if we need
    start_day = [1, 1, 1, 1]
    end_day = [28, 31, 31, 30]

    # Find the first output file
    for file in listdir(directory):
        if file.startswith('iceh.') and file.endswith('.nc'):
            # Read the grid
            id = Dataset(directory + file, 'r')
            lon = id.variables['TLON'][:,:]
            lat = id.variables['TLAT'][:,:]
            id.close()
            break

    # Set up arrays to integrate seasonal climatology of aice and hi
    seasonal_aice = ma.empty([4, size(lon,0), size(lon,1)])
    seasonal_hi = ma.empty([4, size(lon,0), size(lon,1)])
    seasonal_aice[:,:,:] = 0.0
    seasonal_hi[:,:,:] = 0.0
    # Also integrate number of days in each season
    ndays = zeros(4)
    # Total number of days for a final check
    total_days = 0

    # Loop over files
    for file in listdir(directory):
        if file.startswith('iceh.') and file.endswith('.nc'):
            print 'Processing ' + directory + file
            # Read aice and hi
            id = Dataset(directory+file, 'r')
            aice = id.variables['aice'][0,:,:]
            hi = id.variables['hi'][0,:,:]
            # Read the time value and convert to Date object
            time_id = id.variables['time']
            time = num2date(time_id[0], units=time_id.units, calendar=time_id.calendar.lower())
            id.close()
            total_days += 5
            # 5-day averages marked with next day's date
            year = time.year
            month = time.month
            day = time.day
            # Get the season of this day
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
                if day-5 < start_day[season]:
                    # Spills over into the previous season
                    prev_season = mod(season-1, 4)
                    # How many days does it spill over by?
                    spill_days = start_day[season]-day+5
                    # Should be between 1 and 5
                    if spill_days < 1 or spill_days > 5:
                        print 'Problem: spill_days is ' + str(spill_days)
                        print 'Timestep ' + str(t+1)
                        print 'Year ' + str(year)
                        print 'Month ' + str(month+1)
                        print 'Day ' + str(day)
                        return
                    # Split between previous season and this season
                    seasonal_aice[prev_season,:,:] += aice*spill_days
                    seasonal_hi[prev_season,:,:] += hi*spill_days
                    ndays[prev_season] += spill_days
                    seasonal_aice[season,:,:] += aice*(5-spill_days)
                    seasonal_hi[season,:,:] += hi*(5-spill_days)
                    ndays[season] += 5-spill_days
                else:
                    # Entirely within the season
                    seasonal_aice[season,:,:] += aice*5
                    seasonal_hi[season,:,:] += hi*5
                    ndays[season] += 5
            else:
                # Not in the first month of the season
                # The 5 days will be entirely within the season
                seasonal_aice[season,:,:] += aice*5
                seasonal_hi[season,:,:] += hi*5
                ndays[season] += 5
    # Convert from sums to averages
    for season in range(4):
        seasonal_aice[season,:,:] = seasonal_aice[season,:,:]/ndays[season]
        seasonal_hi[season,:,:] = seasonal_hi[season,:,:]/ndays[season]
    if sum(ndays) != total_days:
        print 'Problem: files have ' + str(total_days) + ' days, but climatology has ' + str(sum(ndays))
        return

    # Write to file
    print 'Writing ' + out_file
    id = Dataset(out_file, 'w')
    id.createDimension('ni', size(lon,1))
    id.createDimension('nj', size(lon,0))
    id.createDimension('time', 4)
    id.createVariable('TLON', 'f8', ('nj', 'ni'))
    id.variables['TLON'].long_name = 'T grid center longitude'
    id.variables['TLON'].units = 'degrees_east'
    id.variables['TLON'][:,:] = lon
    id.createVariable('TLAT', 'f8', ('nj', 'ni'))
    id.variables['TLAT'].long_name = 'T grid center latitude'
    id.variables['TLAT'].units = 'degrees_north'
    id.variables['TLAT'][:,:] = lat
    id.createVariable('time', 'f8', ('time'))
    id.variables['time'].units = 'season'
    id.variables['time'].long_name = 'DJF, MAM, JJA, SON'
    id.variables['time'][:] = arange(1,4+1)
    id.createVariable('aice', 'f8', ('time', 'nj', 'ni'))
    id.variables['aice'].units = '1'
    id.variables['aice'][:,:,:] = seasonal_aice
    id.createVariable('hi', 'f8', ('time', 'nj', 'ni'))
    id.variables['hi'].units = 'm'
    id.variables['hi'][:,:,:] = seasonal_hi
    id.close()


# Command-line interface
if __name__ == "__main__":

    directory = raw_input("Path to CICE output history directory: ")
    out_file = raw_input("Path to desired output file: ")
    seasonal_climatology_cice(directory, out_file)
    
            

    

    
