from netCDF4 import Dataset
from numpy import *

def remove_iceshelves (grid_file):

    id = Dataset(grid_file, 'a')
    h = id.variables['h'][:,:]
    mask_rho = id.variables['mask_rho'][:,:]
    mask_zice = id.variables['mask_zice'][:,:]

    index = nonzero(mask_zice==1)
    h[index] = 50
    mask_rho[index] = 0
    mask_u = mask_rho[:,1:]*mask_rho[:,:-1]
    mask_v = mask_rho[1:,:]*mask_rho[:-1,:]
    mask_psi = mask_rho[1:,1:]*mask_rho[1:,:-1]*mask_rho[:-1,1:]*mask_rho[:-1,:-1]

    id.variables['h'][:,:] = h
    id.variables['mask_rho'][:,:] = mask_rho
    id.variables['mask_u'][:,:] = mask_u
    id.variables['mask_v'][:,:] = mask_v
    id.variables['mask_psi'][:,:] = mask_psi

    id.close()


if __name__ == "__main__":

    grid_file = raw_input("Path to ROMS grid file to edit: ")
    remove_iceshelves(grid_file)
    
    
