from netCDF4 import Dataset
from numpy import *
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import KDTree

def romscice_sss_nudging ():

    grid_file = '../metroms_iceshelf/apps/common/grid/circ30S_quarterdegree.nc'
    woa_head = '../metroms_iceshelf/data/originals/WOA_2013/woa13_95A4_s'
    woa_tail = '_01v2.nc'
    output_file = '../metroms_iceshelf/data/sss_nudge.nc'
    nbdry_woa = 61

    print 'Reading ECCO2 grid'
    id = Dataset(woa_head + '01' + woa_tail, 'r')
    lon_woa_raw = id.variables['lon'][:]
    lat_woa = id.variables['lat'][0:nbdry_woa]
    id.close()

    index = lon_woa_raw > 180
    lon_woa_raw[index] -= 360

    # The WOA longitude axis doesn't wrap around; there is a gap between
    # almost-180W and almost-180E, and the ROMS grid has points in this gap.
    # So copy the last longitude value (mod 360) to the beginning, and the
    # first longitude value (mod 360) to the end.
    lon_woa = zeros(size(lon_woa_raw)+2)
    lon_woa[0] = lon_woa_raw[-1]-360
    lon_woa[1:-1] = lon_woa_raw
    lon_woa[-1] = lon_woa_raw[0]+360

    # Read ROMS grid
    print 'Reading ROMS grid'
    id = Dataset(grid_file, 'r')
    lon_roms = id.variables['lon_rho'][:,:]
    lat_roms = id.variables['lat_rho'][:,:]
    h = id.variables['h'][:,:]    
    mask_rho = id.variables['mask_rho'][:,:]
    mask_zice = id.variables['mask_zice'][:,:]
    id.close()
    num_lon = size(lon_roms, 1)
    num_lat = size(lon_roms, 0)

    index = lon_roms > 180
    lon_roms[index] -= 360

    print 'Setting up ' + output_file
    out_id = Dataset(output_file, 'w')
    out_id.createDimension('xi_rho', num_lon)
    out_id.createDimension('eta_rho', num_lat)
    out_id.createDimension('sss_time', None)
    out_id.createVariable('sss_time', 'f8', ('sss_time'))
    out_id.variables['sss_time'].long_name = 'time since initialization'
    out_id.variables['sss_time'].units = 'days'
    out_id.variables['sss_time'].cycle_length = 365.25
    out_id.createVariable('SSS', 'f8', ('sss_time', 'eta_rho', 'xi_rho'))
    out_id.variables['SSS'].long_name = 'surface salinity'
    out_id.variables['SSS'].units = 'psu'

    for month in range(12):
        print 'Processing month ' + str(month+1) + ' of 12'
        if month+1 < 10:
            filename = woa_head + '0' + str(month+1) + woa_tail
        else:
            filename = woa_head + str(month+1) + woa_tail
        id = Dataset(filename, 'r')
        sss_raw = transpose(id.variables['s_an'][0,0,0:nbdry_woa,:])
        id.close()
        sss = ma.array(zeros((size(lon_woa), size(lat_woa))))
        sss[1:-1,:] = ma.copy(sss_raw)
        sss[0,:] = ma.copy(sss_raw[-1,:])
        sss[-1,:] = ma.copy(sss_raw[0,:])
        sss_interp = interp_woa2roms_sfc(sss, lon_woa, lat_woa, lon_roms, lat_roms)

        time = 365.25/12*(month+0.5)
        out_id.variables['sss_time'][month] = time        
        out_id.variables['SSS'][month,:,:] = sss_interp

    out_id.close()


def interp_woa2roms_sfc (A, lon_woa, lat_woa, lon_roms, lat_roms):

    j,i = mgrid[0:A.shape[0], 0:A.shape[1]]
    jigood = array((j[~A.mask], i[~A.mask])).T
    jibad = array((j[A.mask], i[A.mask])).T
    A[A.mask] = A[~A.mask][KDTree(jigood).query(jibad)[1]]
    interp_function = RegularGridInterpolator((lon_woa, lat_woa), A)
    B = interp_function((lon_roms, lat_roms))

    return B


# Command-line interface
if __name__ == "__main__":

    romscice_sss_nudging()
    

    
