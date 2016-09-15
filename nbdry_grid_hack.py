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
        # Find the southernmost unmasked cell within "num_pts" of the
        # northern boundary and set all the points north of it to match
        found_pt = False
        for j in range(num_pts, -1, -1):
            if mask_rho[-j,i] == 1:
                if found_pt:
                    # Already found the right point
                    h[-j,i] = val
                else:
                    # This is the first unmasked point
                    found_pt = True
                    val = h[-j,i]

    # Save changes
    id.variables['h'][:,:] = h
    id.close()


# Command-line interface
if __name__ == "__main__":

    grid_file='../ROMS-CICE-MCT/apps/common/grid/circ30S_quarterdegree_10m.nc'
    num_pts = 15
    nbdry_grid_hack(grid_file, num_pts)
