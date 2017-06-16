from netCDF4 import Dataset, num2date
from numpy import *

# Build a monthly climatology of the extra surface salt flux due to surface
# salinity restoring in a long ROMS run. Save as a NetCDF file that can be
# used as a forcing file for another ROMS run (option SSFLUX_EXTRA) with
# surface salinity restoring off.
# Input:
# in_file = path to ROMS ocean averages file concatenated for the entire
#           simulation (or just the variable ssflux_restoring extracted)
# out_file = desired path to output ROMS forcing file
def romscice_flux_correction (in_file, out_file):

    # Starting and ending days for each month
    # Ignore leap years, they will be dealt with later
    start_day = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    end_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    print 'Reading ' + in_file
    id = Dataset(in_file, 'r')
    lon = id.variables['lon_rho'][:,:]
    lat = id.variables['lat_rho'][:,:]
    # Convert time to Date objects
    time_id = id.variables['ocean_time']
    time = num2date(time_id[:], units=time_id.units, calendar=time_id.calendar.lower())
    restoring = id.variables['ssflux_restoring'][:,:,:]
    id.close()    

    print 'Calculating climatology'
    num_lat = size(lat,0)
    num_lon = size(lon,1)
    # Initialise integral arrays
    climatology = ma.empty([12, num_lat, num_lon])
    climatology[:,:,:] = 0.0
    num_days = zeros(12)
    # Loop over timesteps
    for t in range(size(time)):
        print '...time index ' + str(t+1) + ' of ' + str(size(time))
        # 5-day averages marked with middle day's date
        year = time[t].year
        month = time[t].month-1
        day = time[t].day
        # Check for leap years
        leap_year = False
        if mod(year, 4) == 0:
            leap_year = True
            if mod(year, 100) == 0:
                leap_year = False
                if mod(year, 400) == 0:
                    leap_year = True
        if leap_year:
            end_day[1] = 29
        else:
            end_day[1] = 28
        # Integrate ssflux weighted by days for the correct month(s)
        if day-2 >= start_day[month] and day+2 <= end_day[month]:
            # The 5 days in this index are entirely within 1 month
            climatology[month,:,:] += restoring[t,:,:]*5
            # Also integrate the number of days
            num_days[month] += 5
        elif day-2 < start_day[month]:
            # Partially spills over into the previous month
            prev_month = mod(month-1,12)
            # How many days does it spill over by?
            spill_days = start_day[month]-day+2
            # Should be either 1 or 2
            if spill_days < 1 or spill_days > 2:
                print 'Problem: spill_days is ' + str(spill_days)
                print 'Timestep ' + str(t+1)
                print 'Year ' + str(year)
                print 'Month ' + str(month+1)
                print 'Day ' + str(day)
                return
            # Split between previous month and this month appropriately
            climatology[prev_month,:,:] += restoring[t,:,:]*spill_days
            num_days[prev_month] += spill_days
            climatology[month,:,:] += restoring[t,:,:]*(5-spill_days)
            num_days[month] += 5-spill_days
        elif day+2 > end_day[month]:
            # Partially spills over into the next month
            next_month = mod(month+1,12)
            # How many days does it spill over by?
            spill_days = day+2-end_day[month]
            # Should be either 1 or 2
            if spill_days < 1 or spill_days > 2:
                print 'Problem: spill_days is ' + str(spill_days)
                print 'Timestep ' + str(t+1)
                print 'Year ' + str(year)
                print 'Month ' + str(month+1)
                print 'Day ' + str(day)
                return
            # Split between this month and next month appropriately
            climatology[next_month,:,:] += restoring[t,:,:]*spill_days
            num_days[next_month] += spill_days
            climatology[month,:,:] += restoring[t,:,:]*(5-spill_days)
            num_days[month] += 5-spill_days
    # Make sure we counted days correctly
    if sum(num_days) != size(time)*5:
        print 'Problem with the number of days: found ' + str(num_days) + ' instead of ' + sum(size(time)*5)
    # Convert from sum to averages
    for month in range(12):
        climatology[month,:,:] /= num_days[month]    

    print 'Writing ' + out_file
    id = Dataset(out_file, 'w')
    id.createDimension('xi_rho', num_lon)
    id.createDimension('eta_rho', num_lat)
    id.createDimension('time', None)
    id.createVariable('lon_rho', 'f8', ('eta_rho', 'xi_rho'))
    id.variables['lon_rho'].long_name = 'longitude of rho-points'
    id.variables['lon_rho'].units = 'degree_east'
    id.variables['lon_rho'][:,:] = lon
    id.createVariable('lat_rho', 'f8', ('eta_rho', 'xi_rho'))
    id.variables['lat_rho'].long_name = 'latitude of rho-points'
    id.variables['lat_rho'].units = 'degree_north'
    id.variables['lat_rho'][:,:] = lat
    id.createVariable('time', 'f8', ('time'))
    id.variables['time'].units = 'days since 1992-01-01 00:00:0.0'
    # Add a cycle length so it repeats every year
    id.variables['time'].cycle_length = 365.25
    # Time values in the middle of each month
    id.variables['time'][:] = 365.25/12*(arange(12)+0.5)
    id.createVariable('ssflux_extra', 'f8', ('time', 'eta_rho', 'xi_rho'))
    id.variables['ssflux_extra'].long_name = 'additional surface salt flux'
    id.variables['ssflux_extra'].units = 'psu m/s'
    id.variables['ssflux_extra'][:,:,:] = climatology
    id.close()


# Command-line interface
if __name__ == "__main__":

    in_file = raw_input("Path to ROMS ocean averages file containing ssflux_restoring for the entire simulation: ")
    out_file = raw_input("Path to desired output forcing file: ")
    romscice_flux_correction(in_file, out_file)
