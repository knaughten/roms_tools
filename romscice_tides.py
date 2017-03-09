from netCDF4 import Dataset
from numpy import *
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import KDTree

# Create a ROMS tide file containing the first 10 tidal components interpolated
# from TPXO 7.2. 

# NB for raijin users: RegularGridInterpolator needs python/2.7.6 but the
# default is 2.7.3. Before running this script, switch them as follows:
# module unload python/2.7.3
# module unload python/2.7.3-matplotlib
# module load python/2.7.6
# module load python/2.7.6-matplotlib

def romscice_tides ():

    # Path to ROMS grid file
    grid_file = '../ROMS-CICE-MCT/apps/common/grid/circ30S_quarterdegree_10m.nc'
    # Path to TPXO file
    tpxo_file = '../ROMS-CICE-MCT/data/h_tpxo7.2.nc'
    # Desired path to output file
    out_file = '../ROMS-CICE-MCT/data/tides_tpxo72.nc' #_1yr.nc'
    # Bounds on latitude indices to read from TPXO2 (as close together as
    # possible while containing all the latitudes ROMS needs, so that there
    # aren't too many land points to fill with nearest neighbours which is
    # quite a slow process)
    tpxo_sbdry = 11
    tpxo_nbdry = 241

    # Number of tidal components to read
    num_cmp = 10
    # String containing names of tidal components (just for output file global
    # attribute)
    cmp_names = 'm2,s2,n2,k2,k1,o1,p1,q1,mf,mm'
    # Periods of each tidal component in seconds (lifted from old ROMS tide
    # file caisom001_tides.nc)
    period = array([44714.165191868, 43200.0012869521, 86164.0770050671, 92949.6357005365, 45570.0535117177, 86637.1997716528, 43082.0503185947, 96726.0857029666, 2380715.86358729, 1180295.54554976])    
    # Tweak all these periods so they evenly divide 1 year
    #period_1yr = zeros(num_cmp)
    #sec_per_year = 365.25*24*60*60
    #for n in range(num_cmp):
        #period_orig = period[n]
        #freq = round(sec_per_year/period_orig)
        #period_1yr[n] = sec_per_year/freq

    # Read ROMS grid
    id = Dataset(grid_file, 'r')
    lon_roms = id.variables['lon_rho'][:,:]
    lat_roms = id.variables['lat_rho'][:,:]
    mask = id.variables['mask_rho'][:,:]
    id.close()
    num_lon = size(lon_roms, 1)
    num_lat = size(lon_roms, 0)

    # Read TPXO grid, tidal phase, and tidal amplitude
    id = Dataset(tpxo_file, 'r')
    lon_tpxo = id.variables['lon_z'][:,0]
    lat_tpxo = id.variables['lat_z'][0,tpxo_sbdry:tpxo_nbdry]
    Ephase_tpxo = id.variables['hp'][0:num_cmp,:,tpxo_sbdry:tpxo_nbdry]
    Eamp_tpxo = id.variables['ha'][0:num_cmp,:,tpxo_sbdry:tpxo_nbdry]    
    id.close()

    # Interpolate phase and amplitude to ROMS grid
    Ephase_interp = ma.empty([num_cmp, num_lat, num_lon])
    Eamp_interp = ma.empty([num_cmp, num_lat, num_lon])
    # Loop over components
    for n in range(num_cmp):
        print 'Component ' + str(n+1) + ' of ' + str(num_cmp)
        print 'Processing tidal elevation phase angle'
        tmp = interp_tpxo2roms(Ephase_tpxo[n,:,:], lon_tpxo, lat_tpxo, lon_roms, lat_roms)
        # Fill ROMS land mask with zeros
        index = mask==0
        tmp[index] = 0.0
        Ephase_interp[n,:,:] = tmp
        print 'Processing tidal elevation amplitude'
        tmp = interp_tpxo2roms(Eamp_tpxo[n,:,:], lon_tpxo, lat_tpxo, lon_roms, lat_roms)
        tmp[index] = 0.0
        Eamp_interp[n,:,:] = tmp

    # Output to NetCDF file
    print 'Writing ' + out_file
    id = Dataset(out_file, 'w')
    id.createDimension('xi_rho', num_lon)
    id.createDimension('eta_rho', num_lat)
    id.createDimension('tide_period', num_cmp)
    id.createVariable('tide_period', 'f8', ('tide_period'))
    id.variables['tide_period'].long_name = 'tide angular period'
    id.variables['tide_period'].units = 'seconds'
    id.variables['tide_period'][:] = period #_1yr
    id.createVariable('tide_Ephase', 'f8', ('tide_period', 'eta_rho', 'xi_rho'))
    id.variables['tide_Ephase'].long_name = 'tidal elevation phase angle'
    id.variables['tide_Ephase'].units = 'degrees, time of maximum elevation with respect to chosen time origin'
    id.variables['tide_Ephase'][:,:,:] = Ephase_interp
    id.createVariable('tide_Eamp', 'f8', ('tide_period', 'eta_rho', 'xi_rho'))
    id.variables['tide_Eamp'].long_name = 'tidal elevation amplitude'
    id.variables['tide_Eamp'].units = 'meter'
    id.variables['tide_Eamp'][:,:,:] = Eamp_interp
    id.components = cmp_names
    id.close()


# Given a field on the TPXO grid, interpolate to the ROMS grid.
# Input:
# data_tpxo = 2D array (size mxn) containing data on the TPXO grid
# lon_tpxo = 1D array (size m) containing TPXO longitude values
# lat_tpxo = 1D array (size n) containing TPXO latitude values
# lon_roms = 2D array (size pxq) containing ROMS longitude values
# lat_roms = 2D array (size pxq) containing ROMS latitude values
# Output:
# data_interp = 2D array (size pxq) containing data interpolate to the ROMS grid
def interp_tpxo2roms(data_tpxo, lon_tpxo, lat_tpxo, lon_roms, lat_roms):

    # TPXO data is stored sideways in NetCDF files (lon x lat instead of
    # lat x lon); fix this
    data_tpxo = transpose(data_tpxo)
    # TPXO data is masked with zeros; convert this to a proper mask
    data_tpxo = ma.masked_where(data_tpxo==0, data_tpxo)

    # Fill in masked values with nearest neighbours so they don't screw up
    # the interpolation
    # I got this code from Stack Exchange, not really sure how it works
    print 'Filling land mask with nearest neighbours'
    j,i = mgrid[0:data_tpxo.shape[0], 0:data_tpxo.shape[1]]
    jigood = array((j[~data_tpxo.mask], i[~data_tpxo.mask])).T
    jibad = array((j[data_tpxo.mask], i[data_tpxo.mask])).T
    data_tpxo[data_tpxo.mask] = data_tpxo[~data_tpxo.mask][KDTree(jigood).query(jibad)[1]]

    # TPXO axis doesn't wrap around; there is a gap between longitudes 0.25
    # and 360, but ROMS needs to interpolate in this gap.
    # So copy the last longitude value (mod 360) to the beginning, and the
    # first longitude value (mod 360) to the end.
    lon_tpxo_wrap = zeros(size(lon_tpxo)+2)
    lon_tpxo_wrap[0] = lon_tpxo[-1] - 360
    lon_tpxo_wrap[1:-1] = lon_tpxo
    lon_tpxo_wrap[-1] = lon_tpxo[0] + 360
    lon_tpxo = lon_tpxo_wrap
    # Copy the westernmost and easternmost data points to match
    data_tpxo_wrap = ma.array(zeros((size(lat_tpxo), size(lon_tpxo))))
    data_tpxo_wrap[:,1:-1] = data_tpxo
    data_tpxo_wrap[:,0] = data_tpxo[:,-1]
    data_tpxo_wrap[:,-1] = data_tpxo[:,0]
    data_tpxo = data_tpxo_wrap

    print 'Interpolating to ROMS grid'
    # Build an interpolation function for TPXO data
    interp_function = RegularGridInterpolator((lat_tpxo, lon_tpxo), data_tpxo)
    # Call it for the ROMS grid
    data_interp = interp_function((lat_roms, lon_roms))

    return data_interp


# Command-line interface
if __name__ == "__main__":

    romscice_tides()
    
