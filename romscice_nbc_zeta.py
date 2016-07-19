from netCDF4 import Dataset
from numpy import *
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import KDTree

# Interpolate AVISO sea surface height climatology data (annual average) 
# to the northern boundary of the ROMS grid. Save to the given lateral boundary
# condition file, created using romscice_nbc.py.

# NB for raijin users: RegularGridInterpolator needs python/2.7.6 but the
# default is 2.7.3. Before running this script, switch them as follows:
# module unload python/2.7.3
# module unload python/2.7.3-matplotlib
# module load python/2.7.6
# module load python/2.7.6-matplotlib

def romscice_nbc_zeta (nbc_file):

    # Paths of ROMS grid file and input AVISO file
    grid_file = '../ROMS-CICE-MCT/apps/common/grid/circ30S_quarterdegree_10m.nc'
    aviso_file = '../ROMS-CICE-MCT/data/aviso_climatology.nc'

    # Read ROMS grid at northern boundary
    id = Dataset(grid_file, 'r')
    lon_rho = id.variables['lon_rho'][-1,:]
    lat_rho = id.variables['lat_rho'][-1,:]
    id.close()

    # Read AVISO data
    id = Dataset(aviso_file, 'r')
    lon_aviso_raw = id.variables['lon'][:]
    lat_aviso = id.variables['lat'][59:61]
    ssh_aviso_raw = transpose(id.variables['zos'][0,59:61,:])
    id.close()

    # The AVISO longitude axis doesn't wrap around; there is a gap between
    # almost-180W and almost-180E, and the ROMS grid has points in this gap.
    # So copy the last longitude value (mod 360) to the beginning, and the
    # first longitude value (mod 360) to the end.
    lon_aviso = zeros(size(lon_aviso_raw)+2)
    lon_aviso[0] = lon_aviso_raw[-1]-360
    lon_aviso[1:-1] = lon_aviso_raw
    lon_aviso[-1] = lon_aviso_raw[0]+360

    # Copy the SSH data to the new indices, making sure to preserve the mask
    ssh_aviso = ma.empty((size(lon_aviso), size(lat_aviso)))
    ssh_aviso[1:-1,:] = ssh_aviso_raw
    ssh_aviso[0,:] = ssh_aviso_raw[-1,:]
    ssh_aviso[-1,:] = ssh_aviso_raw[0,:]

    # Fill land  mask with nearest neighbours
    j,i = mgrid[0:ssh_aviso.shape[0], 0:ssh_aviso.shape[1]]
    jigood = array((j[~ssh_aviso.mask], i[~ssh_aviso.mask])).T
    jibad = array((j[ssh_aviso.mask], i[ssh_aviso.mask])).T
    ssh_aviso[ssh_aviso.mask] = ssh_aviso[~ssh_aviso.mask][KDTree(jigood).query(jibad)[1]]

    # Build a function for linear interpolation of SSH
    interp_function = RegularGridInterpolator((lon_aviso, lat_aviso), ssh_aviso)
    # Call this function at ROMS grid points
    zeta = interp_function((lon_rho, lat_rho))

    # Save to lateral BC file
    id = Dataset(nbc_file, 'a')
    for month in range(12):
        id.variables['zeta_north'][month,:] = zeta
    id.close()


if __name__ == "__main__":
    nbc_file = raw_input("Path to existing ROMS northern boundary condition file: ")
    romscice_nbc_zeta(nbc_file)
    
        
        
    
    
    
