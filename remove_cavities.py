from netCDF4 import Dataset
from numpy import *

# Remove the ice shelf cavities from the given ROMS grid file and fill them
# in with land. After you run this script, I suggest you remove the variables
# "zice" and "mask_zice" from the file using NCO, to avoid any confusion later.
# Input:
# grid_path = path to ROMS grid file to edit
def remove_cavities (grid_path):

    # Read existing grid
    id = Dataset(grid_path, 'a')
    mask_zice = id.variables['mask_zice'][:,:]
    mask_rho = id.variables['mask_rho'][:,:]
    mask_u = id.variables['mask_u'][:,:]
    mask_v = id.variables['mask_v'][:,:]
    mask_psi = id.variables['mask_psi'][:,:]
    h = id.variables['h'][:,:]

    # Select ice shelf points
    index = mask_zice == 1
    # Apply land mask
    mask_rho[index] = 0
    # Apply token bathymetry
    h[index] = 50

    # Recalculate masks on u-, v-, and psi-grids
    mask_u = mask_rho[:,:-1]*mask_rho[:,1:]
    mask_v = mask_rho[:-1,:]*mask_rho[1:,:]
    mask_psi = mask_rho[1:,1:]*mask_rho[1:,:-1]*mask_rho[:-1,1:]*mask_rho[:-1,:-1]

    # Overwrite variables in file
    id.variables['mask_rho'][:,:] = mask_rho
    id.variables['mask_u'][:,:] = mask_u
    id.variables['mask_v'][:,:] = mask_v
    id.variables['mask_psi'][:,:] = mask_psi
    id.variables['h'][:,:] = h

    id.close()


# Command-line interface
if __name__ == "__main__":

    grid_path = raw_input("Path to ROMS grid file to edit: ")
    remove_cavities(grid_path)

    
