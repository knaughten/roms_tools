from netCDF4 import Dataset
from numpy import *
from scipy.interpolate import griddata

# Make a CICE initial conditions file (non-standard option ice_ic='mask') based
# on NSIDC satellite observations of sea ice concentration. Wherever 
# concentration exceeds 15%, set it to 100%; everywhere else, set it to 0%. 
# The CICE code will read this "mask" and set a constant sea ice thickness of 
# 1 m and snow thickness of 0.2 m.
# Input:
# nsidc_file = path to NSIDC file: CDR monthly average of sea ice concentration
# roms_grid = path to ROMS grid file
# out_file = desired path to output file
def cice_ini (nsidc_file, roms_grid, out_file):

    # Read the NSIDC grid and data
    id = Dataset(nsidc_file, 'r')
    nsidc_lon = id.variables['longitude'][:,:]
    nsidc_lat = id.variables['latitude'][:,:]
    nsidc_aice = id.variables['seaice_conc_monthly_cdr'][0,:,:]
    id.close()

    # Make sure longitude goes from 0 to 360 to fit with ROMS convention
    index = nsidc_lon < 0
    nsidc_lon[index] = nsidc_lon[index] + 360

    # Make sure sea ice concentration values go from 0 to 1
    index = nsidc_aice < 0
    nsidc_aice[index] = 0.0
    index = nsidc_aice > 1
    nsidc_aice[index] = 1.0

    # Read the ROMS grid
    id = Dataset(roms_grid, 'r')
    roms_lon = id.variables['lon_rho'][:,:]
    roms_lat = id.variables['lat_rho'][:,:]
    id.close()

    # Make an nx2 array of the NSIDC coordinates, flattened
    points = empty([size(nsidc_lon), 2])
    points[:,0] = ravel(nsidc_lon)
    points[:,1] = ravel(nsidc_lat)
    # Also flatten the NSIDC data
    values = ravel(nsidc_aice)
    # Now make an mx2 array of the ROMS coordinates
    xi = empty([size(roms_lon), 2])
    xi[:,0] = ravel(roms_lon)
    xi[:,1] = ravel(roms_lat)
    # Linear interpolation to the ROMS coordinates
    result = griddata(points, values, xi, method='linear', fill_value=-999)
    # Also do a nearest-neighbour interpolation
    result2 = griddata(points, values, xi, method='nearest')
    # Replace missing values (from land mask) in linear interpolation with
    # nearest-neighbour interpolation so there are no coastal artifacts
    index = result==-999
    result[index] = result2[index]
    # Reshape to the correct dimensions
    aice = reshape(result, shape(roms_lon))

    # Cut-off concentration of 15%
    index = aice > 0.15
    aice[index] = 1.0
    aice[~index] = 0.0

    # Remove the edges (CICE only wants interior points)
    aice = aice[1:-1,1:-1]
    num_lon = size(aice, 1)
    num_lat = size(aice, 0)

    # Write output
    out_id = Dataset(out_file, 'w')
    out_id.createDimension('xi_rho', num_lon)
    out_id.createDimension('eta_rho', num_lat)
    out_id.createDimension('time', None)
    out_id.createVariable('time', 'f8', ('time'))
    out_id.variables['time'].units = 'days since 1992-01-01 00:00:0.0'
    out_id.variables['time'][0] = 0.0
    out_id.createVariable('aice', 'f8', ('time', 'eta_rho', 'xi_rho'))
    out_id.variables['aice'].long_name = 'sea ice concentration'
    out_id.variables['aice'].units = '1'
    out_id.variables['aice'][0,:,:] = aice
    out_id.close()


# Command-line interface
if __name__ == "__main__":

    nsidc_file = '/short/m68/kaa561/nsidc_aice/seaice_conc_monthly_sh_f11_199201_v02r00.nc'
    roms_grid = '/short/m68/kaa561/metroms_iceshelf/apps/common/grid/circ30S_quarterdegree.nc'
    out_file = '/short/m68/kaa561/metroms_iceshelf/data/cice_ini.nc'
    cice_ini(nsidc_file, roms_grid, out_file)
