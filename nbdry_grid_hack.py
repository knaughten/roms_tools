from netCDF4 import Dataset
from numpy import *

# Modify the ROMS grid so that the northernmost 15 rows
# (about 3 degrees) have the same bathymetry, i.e.
# dh/dy = 0, on advice of Matt Mazloff.
# Input:
# grid_file = path to ROMS grid file
# num_pts = number of rows to modify; they will all have
#           the bathymetry in the row "num_pts" from the
#           northern boundary
def nbdry_grid_hack (grid_file, num_pts):

    # Read bathymetry and masks
    id = Dataset(grid_file, 'a')
    h = id.variables['h'][:,:]
    mask_rho = id.variables['mask_rho'][:,:]
    mask_u = id.variables['mask_u'][:,:]
    mask_v = id.variables['mask_v'][:,:]
    mask_psi = id.variables['mask_psi'][:,:]

    # Loop over longitude
    for i in range(size(h,1)):
        # Only modify columns where the northern boundary is open
        # (i.e. not land)
        if mask_rho[-1,i] == 1:
            # Set bathymetry and land mask to be the same as the
            # point "num_pts" from the northern boundary
            for j in range(1, num_pts):
                h[-j,i] = h[-num_pts,i]
                mask_rho[-j,i] = mask_rho[-num_pts,i]

    # Update the other masks
    mask_u = mask_rho[:,1:]*mask_rho[:,:-1]
    mask_v = mask_rho[1:,:]*mask_rho[:-1,:]
    mask_psi = mask_rho[1:,1:]*mask_rho[1:,:-1]*mask_rho[:-1,1:]*mask_rho[:-1,:-1]

    # Save changes
    id.variables['h'][:,:] = h
    id.variables['mask_rho'][:,:] = mask_rho
    id.variables['mask_u'][:,:] = mask_u
    id.variables['mask_v'][:,:] = mask_v
    id.variables['mask_psi'][:,:] = mask_psi

    id.close()


# Command-line interface
if __name__ == "__main__":

    grid_file='../ROMS-CICE-MCT/apps/common/grid/circ30S_quarterdegree_10m.nc'
    num_pts = 15
    nbdry_grid_hack(grid_file, num_pts)
