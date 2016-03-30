from netCDF4 import Dataset
from numpy import *
from scipy.interpolate import RegularGridInterpolator

# Interpolate one year of monthly AVISO sea surface height data to the
# northern boundary of the ROMS grid. Save to the existing lateral boundary
# condition file, created using romscice_nbc.py.

# NB for raijin users: RegularGridInterpolator needs python/2.7.6 but the
# default is 2.7.3. Before running this script, switch them as follows:
# module unload python/2.7.3
# module unload python/2.7.3-matplotlib
# module load python/2.7.6
# module load python/2.7.6-matplotlib

def romscice_nbc_zeta (year):

    # Paths of ROMS grid file, input AVISO file (without the tail mm.nc), and
    # output ROMS-CICE boundary condition file (existing; created using
    # romscice_nbc.py)
    grid_file = '../ROMS-CICE-MCT/apps/common/grid/circ38S_quarterdegree.nc'
    aviso_base = '../ROMS-CICE-MCT/data/AVISO/dt_global_allsat_msla_h_y' + str(year) + '_m'
    nbc_file = '../ROMS-CICE-MCT/data/ECCO2/ecco2_cube92_lbc_' + str(year) + '.nc'
    # Northernmost index of AVISO to read (1-based)
    nbdry_aviso = 209

    # Read AVISO grid
    aviso_fid = Dataset(aviso_base + '01.nc', 'r')
    lon_aviso_raw = aviso_fid.variables['lon'][:]
    lat_aviso = aviso_fid.variables['lat'][0:nbdry_aviso]
    aviso_fid.close()

    # The AVISO longitude axis doesn't wrap around; there is a gap between
    # almost-180W and almost-180E, and the ROMS grid has points in this gap.
    # So copy the last longitude value (mod 360) to the beginning, and the
    # first longitude value (mod 360) to the end.
    lon_aviso = zeros(size(lon_aviso_raw)+2)
    lon_aviso[0] = lon_aviso_raw[-1]-360
    lon_aviso[1:-1] = lon_aviso_raw
    lon_aviso[-1] = lon_aviso_raw[0]+360

    # Read ROMS grid
    grid_fid = Dataset(grid_file, 'r')
    lon_rho = grid_fid.variables['lon_rho'][-1,:]
    lat_rho = grid_fid.variables['lat_rho'][-1,:]
    grid_fid.close()

    # Loop through each month of this year
    for month in range(12):

        print 'Processing month ', str(month+1), ' of 12'
        # Construct the rest of the AVISO file path
        if month+1 < 10:
            tail = '0' + str(month+1) + '.nc'
        else:
            tail = str(month+1) + '.nc'

        # Read SSH data 
        aviso_fid = Dataset(aviso_base + tail, 'r')
        sla_raw = transpose(aviso_fid.variables['sla'][0,0:nbdry_aviso,:])
        aviso_fid.close()

        # Copy the data to the new longitude indices, making sure to
        # preserve the mask
        sla = ma.array(zeros((size(lon_aviso), size(lat_aviso))))
        sla[1:-1,:] = ma.copy(sla_raw)
        sla[0,:] = ma.copy(sla_raw[-1,:])
        sla[-1,:] = ma.copy(sla_raw[0,:])

        # Fill mask with zeros since RegularGridInterpolator can't handle masks
        sla[sla.mask] = 0.0
        # Build a function for linear interpolation of sla
        interp_function = RegularGridInterpolator((lon_aviso, lat_aviso), sla)
        # Call this function at ROMS grid points
        zeta = interp_function((lon_rho, lat_rho))

        # Save to lateral BC file
        nbc_fid = Dataset(nbc_file, 'a')
        nbc_fid.variables['zeta_north'][month,:] = zeta
        nbc_fid.close()
    
        
        
    
    
    
