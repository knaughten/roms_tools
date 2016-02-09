from netCDF4 import Dataset
from numpy import *
from scipy.interpolate import LinearNDInterpolator, RectBivariateSpline

# Convert two ERA-Interim files:
# AN_yyyy_monthly_orig.nc: one year of monthly averaged measurements for
#                          surface pressure (sp), 2-metre temperature (t2m)
#                          and dew point (d2m), total cloud cover (tcc), and
#                          10-metre winds (u10, v10)
# FC_yyyy_monthly_orig.nc: one year of monthly averaged measurements for
#                          total precipitation (tp) and snowfall (sf)
# to two ROMS-CICE input forcing files with the correct units and naming
# conventions:
# AN_yyyy_monthly.nc: one year of monthly averaged measurements for surface
#                     pressure (Pair), temperature (Tair), specific humidity
#                     (Qair), cloud fraction (cloud), and winds (Uwind, Vwind)
# FC_yyyy_monthly.nc: one year of monthly averaged measurements for rainfall
#                     (rain) and snowfall (snow)
# Input: year = integer containing the year to process
def convert_file (year):

    # Paths of ROMS grid file, input ERA-Interim files, and output ROMS-CICE
    # files; other users will need to change these
    grid_file = '../ROMS-CICE-MCT/apps/common/grid/circ38S_quarterdegree.nc'
    input_atm_file = '../ROMS-CICE-MCT/data/ERA_Interim/monthly/originals/AN_' + str(year) + '_monthly_orig.nc'
    input_ppt_file = '../ROMS-CICE-MCT/data/ERA_Interim/monthly/originals/FC_' + str(year) + '_monthly_orig.nc'
    output_atm_file = '../ROMS-CICE-MCT/data/ERA_Interim/monthly/AN_' + str(year) + '_monthly.nc'
    output_ppt_file = '../ROMS-CICE-MCT/data/ERA_Interim/monthly/FC_' + str(year) + '_monthly.nc'

    Lv = 2.5e6 # Latent heat of vapourisation, J/kg
    Rv = 461.5 # Ideal gas constant for water vapour, J/K/kg

    print 'Reading grids'

    # Read ROMS latitude and longitude
    grid_fid = Dataset(grid_file, 'r')
    lon_roms = grid_fid.variables['lon_rho'][:,:]
    lat_roms = grid_fid.variables['lat_rho'][:,:]
    grid_fid.close()
    num_lon = size(lon_roms, 1)
    num_lat = size(lon_roms, 0)

    # Open input AN file and read ERA-Interim grid
    iatm_fid = Dataset(input_atm_file, 'r')
    lon_era = iatm_fid.variables['longitude'][:]
    lat_era = iatm_fid.variables['latitude'][:]
    iatm_fid.close()

    # Create time axis: 12 equally spaced values throughout the year,
    # units of 'days since 1992-01-01 00:00:0.0'
    time = (year-1992)*365.25 + (arange(12) + 0.5)/12*365.25

    print 'Setting up ' + output_atm_file
    oatm_fid = Dataset(output_atm_file, 'w')
    # Define dimensions (note unlimited time dimension)
    oatm_fid.createDimension('xi_rho', num_lon)
    oatm_fid.createDimension('eta_rho', num_lat)
    oatm_fid.createDimension('time', None)
    # Define variables; write latitude and longitude since they are not
    # time-dependent
    oatm_fid.createVariable('lon_rho', 'f8', ('eta_rho', 'xi_rho'))
    oatm_fid.variables['lon_rho'].long_name = 'longitude of rho-points'
    oatm_fid.variables['lon_rho'].units = 'degree_east'
    oatm_fid.variables['lon_rho'][:,:] = lon_roms
    oatm_fid.createVariable('lat_rho', 'f8', ('eta_rho', 'xi_rho'))
    oatm_fid.variables['lat_rho'].long_name = 'latitude of rho-points'
    oatm_fid.variables['lat_rho'].units = 'degree_north'
    oatm_fid.variables['lat_rho'][:,:] = lat_roms
    oatm_fid.createVariable('time', 'f8', ('time'))
    oatm_fid.variables['time'].units = 'days since 1992-01-01 00:00:0.0'
    oatm_fid.createVariable('Pair', 'f8', ('time', 'eta_rho', 'xi_rho'))
    oatm_fid.variables['Pair'].long_name = 'surface air pressure'
    oatm_fid.variables['Pair'].units = 'Pascal'
    oatm_fid.createVariable('Tair', 'f8', ('time', 'eta_rho', 'xi_rho'))
    oatm_fid.variables['Tair'].long_name = 'surface air temperature'
    oatm_fid.variables['Tair'].units = 'Celsius'
    oatm_fid.createVariable('Qair', 'f8', ('time', 'eta_rho', 'xi_rho'))
    oatm_fid.variables['Qair'].long_name = 'surface relative humidity'
    oatm_fid.variables['Qair'].units = 'kg/kg'
    oatm_fid.createVariable('cloud', 'f8', ('time', 'eta_rho', 'xi_rho'))
    oatm_fid.variables['cloud'].long_name = 'cloud fraction'
    oatm_fid.variables['cloud'].units = 'nondimensional'
    oatm_fid.createVariable('Uwind', 'f8', ('time', 'eta_rho', 'xi_rho'))
    oatm_fid.variables['Uwind'].long_name = 'surface u-wind component'
    oatm_fid.variables['Uwind'].units = 'm/s'
    oatm_fid.createVariable('Vwind', 'f8', ('time', 'eta_rho', 'xi_rho'))
    oatm_fid.variables['Vwind'].long_name = 'surface v-wind component'
    oatm_fid.variables['Vwind'].units = 'm/s'
    oatm_fid.close()

    # Process one timestep at a time to minimise memory use
    for t in range(size(time)):
        oatm_fid = Dataset(output_atm_file, 'a')
        print 'Processing record ' + str(t+1) + ' of ' + str(size(time))
        # Write the current time value to output AN file
        oatm_fid.variables['time'][t] = time[t]
        # Read variables for this timestep
        iatm_fid = Dataset(input_atm_file, 'r')
        sp = transpose(iatm_fid.variables['sp'][t,:,:])
        t2m = transpose(iatm_fid.variables['t2m'][t,:,:])
        d2m = transpose(iatm_fid.variables['d2m'][t,:,:])
        tcc = transpose(iatm_fid.variables['tcc'][t,:,:])
        u10 = transpose(iatm_fid.variables['u10'][t,:,:])
        v10 = transpose(iatm_fid.variables['v10'][t,:,:])
        iatm_fid.close()
        # Calculate relative humidity from temperature and dew point
        rh = exp(Lv/Rv*(t2m**(-1) - d2m**(-1)))        
        # Interpolate each variable to ROMS grid and write to output AN file
        pair = interp_era2roms(sp, lon_era, lat_era, lon_roms, lat_roms)
        oatm_fid.variables['Pair'][t,:,:] = pair
        tair = interp_era2roms(t2m, lon_era, lat_era, lon_roms, lat_roms)
        oatm_fid.variables['Tair'][t,:,:] = tair-273.15
        qair = interp_era2roms(rh, lon_era, lat_era, lon_roms, lat_roms)
        # Constrain humidity to be between 0 and 1
        qair[qair < 0] = 0.0
        qair[qair > 1] = 1.0
        oatm_fid.variables['Qair'][t,:,:] = qair
        cloud = interp_era2roms(tcc, lon_era, lat_era, lon_roms, lat_roms)
        # Constrain cloud fractions to be between 0 and 1
        cloud[cloud < 0] = 0.0
        cloud[cloud > 1] = 1.0
        oatm_fid.variables['cloud'][t,:,:] = cloud
        uwind = interp_era2roms(u10, lon_era, lat_era, lon_roms, lat_roms)
        oatm_fid.variables['Uwind'][t,:,:] = uwind
        vwind = interp_era2roms(v10, lon_era, lat_era, lon_roms, lat_roms)
        oatm_fid.variables['Vwind'][t,:,:] = vwind
        oatm_fid.close()

    print 'Setting up ' + output_ppt_file
    oppt_fid = Dataset(output_ppt_file, 'w')
    # Define dimensions
    oppt_fid.createDimension('xi_rho', num_lon)
    oppt_fid.createDimension('eta_rho', num_lat)
    oppt_fid.createDimension('time', None)
    # Define variables
    oppt_fid.createVariable('lon_rho', 'f8', ('eta_rho', 'xi_rho'))
    oppt_fid.variables['lon_rho'].long_name = 'longitude of rho-points'
    oppt_fid.variables['lon_rho'].units = 'degree_east'
    oppt_fid.variables['lon_rho'][:,:] = lon_roms
    oppt_fid.createVariable('lat_rho', 'f8', ('eta_rho', 'xi_rho'))
    oppt_fid.variables['lat_rho'].long_name = 'latitude of rho-points'
    oppt_fid.variables['lat_rho'].units = 'degree_north'
    oppt_fid.variables['lat_rho'][:,:] = lat_roms
    oppt_fid.createVariable('time', 'f8', ('time'))
    oppt_fid.variables['time'].units = 'days since 1992-01-01 00:00:0.0'
    oppt_fid.createVariable('rain', 'f8', ('time', 'eta_rho', 'xi_rho'))
    oppt_fid.variables['rain'].long_name = 'rain fall rate'
    oppt_fid.variables['rain'].units = 'm_per_12hr'
    oppt_fid.createVariable('snow', 'f8', ('time', 'eta_rho', 'xi_rho'))
    oppt_fid.variables['snow'].long_name = 'snow fall rate'
    oppt_fid.variables['snow'].units = 'm_per_12hr'
    oppt_fid.close()


    for t in range(size(time)):
        oppt_fid = Dataset(output_ppt_file, 'a')
        print 'Processing record ' + str(t+1) + ' of ' + str(size(time))
        # Write the current time value to output FC file
        oppt_fid.variables['time'][t] = time[t]
        # Read data for this timestep
        ippt_fid = Dataset(input_ppt_file, 'r')
        tp = transpose(ippt_fid.variables['tp'][t,:,:])
        sf = transpose(ippt_fid.variables['sf'][t,:,:])
        ippt_fid.close()
        # Interpolate to ROMS grid and write to output FC file
        rain = interp_era2roms(tp, lon_era, lat_era, lon_roms, lat_roms)
        snow = interp_era2roms(sf, lon_era, lat_era, lon_roms, lat_roms)
        # Make sure there are no negative values
        rain[rain < 0] = 0.0
        oppt_fid.variables['rain'][t,:,:] = rain
        snow[snow < 0] = 0.0
        oppt_fid.variables['snow'][t,:,:] = snow
        oppt_fid.close()



# Given an array on the ERA-Interim grid, interpolate any missing values, and
# then interpolate to the ROMS grid.
# Input:
# A = array of size nxm containing values on the ERA-Interim grid (dimension
#     longitude x latitude)
# lon_era = array of length n containing ERA-Interim longitude values
# lat_era = array of length m containing ERA-Interim latitude values
# lon_roms = array of size pxq containing ROMS longitude values
# lat_roms = array of size pxq containing ROMS latitude values
# Output:
# B = array of size pxq containing values on the ROMS grid (dimension latitude x
#     longitude)
def interp_era2roms (A, lon_era, lat_era, lon_roms, lat_roms):

    # Save the sizes of ROMS axes
    num_lon = size(lon_roms, 1)
    num_lat = size(lon_roms, 0)

    # Missing values are something <<0, but these can change when the offset
    # and scale factor attributes are automatically applied. Either way,
    # missing values will be the minimum values in A.
    flag = amin(A)

    # Interpolate missing values
    # I got this bit of code from Stack Exchange
    # It seems to work, not exactly sure how
    valid_mask = A > flag
    coords = array(nonzero(valid_mask)).T
    values = A[valid_mask]
    fill_function = LinearNDInterpolator(coords, values)
    Afill = fill_function(list(ndindex(A.shape))).reshape(A.shape)
    # Fill any still-missing values with the mean
    Afill[isnan(Afill)] = nanmean(Afill)

    # Now interpolate from ERA-Interim grid to ROMS grid
    # First build a function to approximate A with 2D splines
    # Note that latitude has to be multiplied by -1 so both axes are strictly
    # ascending
    interp_function = RectBivariateSpline(lon_era, -lat_era, Afill)
    B = zeros(shape(lon_roms))
    # Call this function for each grid point individually - if you try to do
    # it all at once it throws a MemoryError
    for i in range(num_lon):
        for j in range(num_lat):
            B[j,i] = interp_function(lon_roms[j,i], -lat_roms[j,i])

    return B


# Command-line interface
if __name__ == "__main__":

    # Start and end years; other users will need to change these
    first_year = 1992
    last_year = 2005
    for year in range(first_year, last_year+1):
        print 'Processing '+str(year)
        convert_file(year)
        
        
    
    
        
        
    
    
    
