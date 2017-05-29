from netCDF4 import Dataset
from numpy import *

def fix_isolated_pts ():

    grid_file = '../metroms_iceshelf/apps/common/grid/circ30S_quarterdegree_tmp.nc'

    id = Dataset(grid_file, 'a')
    mask_rho = id.variables['mask_rho'][:,:]
    h = id.variables['h'][:,:]  # Fill value 50
    mask_zice = id.variables['mask_zice'][:,:]
    zice = id.variables['zice'][:,:]

    # Set as ice shelf, average of i-1, j-1, j+1
    i_vals = [1162, 1179, 430, 1184]
    j_vals = [131, 177, 183, 201]
    num_pts = len(i_vals)
    for n in range(num_pts):
        i = i_vals[n]-1
        j = j_vals[n]-1
        mask_zice[j,i] = 1
        zice[j,i] = mean(array([zice[j,i-1], zice[j-1,i], zice[j+1,i]]))

    # Set as ice shelf, average of i-1, i+1, j-1
    i_vals = [1291, 1311, 64, 2, 296, 1426, 1375, 1386, 1422, 1106, 1125, 1054, 406]
    j_vals = [85, 95, 103, 116, 118, 123, 124, 125, 125, 157, 157, 171, 181]
    num_pts = len(i_vals)
    for n in range(num_pts):
        i = i_vals[n]-1
        j = j_vals[n]-1
        mask_zice[j,i] = 1
        zice[j,i] = mean(array([zice[j,i-1], zice[j,i+1], zice[j-1,i]]))

    # Set as ice shelf, average of i-1, i+1
    i_vals = [92, 12]
    j_vals = [93, 113]
    num_pts = len(i_vals)
    for n in range(num_pts):
        i = i_vals[n]-1
        j = j_vals[n]-1
        mask_zice[j,i] = 1
        zice[j,i] = mean(array([zice[j,i-1], zice[j,i+1]]))

    # Set as ice shelf, average of i+1, j-1, j+1
    i_vals = [714, 850, 1131]
    j_vals = [101, 131, 176]
    num_pts = len(i_vals)
    for n in range(num_pts):
        i = i_vals[n]-1
        j = j_vals[n]-1
        mask_zice[j,i] = 1
        zice[j,i] = mean(array([zice[j,i+1], zice[j-1,i], zice[j+1,i]]))

    # Set as ice shelf, average of i+1, j-1
    i_vals = [1110]
    j_vals = [158]
    num_pts = len(i_vals)
    for n in range(num_pts):
        i = i_vals[n]-1
        j = j_vals[n]-1
        mask_zice[j,i] = 1
        zice[j,i] = mean(array([zice[j,i+1], zice[j-1,i]]))

    # Set as ice shelf, average of i+1, j+1
    i_vals = [1023, 1129]
    j_vals = [164, 165]
    num_pts = len(i_vals)
    for n in range(num_pts):
        i = i_vals[n]-1
        j = j_vals[n]-1
        mask_zice[j,i] = 1
        zice[j,i] = mean(array([zice[j,i+1], zice[j+1,i]]))

    # Set as ice shelf, average of i-1, j+1
    i_vals = [1178]
    j_vals = [171]
    num_pts = len(i_vals)
    for n in range(num_pts):
        i = i_vals[n]-1
        j = j_vals[n]-1
        mask_zice[j,i] = 1
        zice[j,i] = mean(array([zice[j,i-1], zice[j+1,i]]))

    # Set as ice shelf, average of i+1, j-1
    i_vals = [665, 1136, 596]
    j_vals = [178, 188, 197]
    num_pts = len(i_vals)
    for n in range(num_pts):
        i = i_vals[n]-1
        j = j_vals[n]-1
        mask_zice[j,i] = 1
        zice[j,i] = mean(array([zice[j,i+1], zice[j-1,i]]))

    # Set as ice shelf, average of i-1, j-1
    i_vals = [516]
    j_vals = [190]
    num_pts = len(i_vals)
    for n in range(num_pts):
        i = i_vals[n]-1
        j = j_vals[n]-1
        mask_zice[j,i] = 1
        zice[j,i] = mean(array([zice[j,i-1], zice[j-1,i]]))

    # Set as ice shelf, same as j-1
    i_vals = [1157, 1122]
    j_vals = [192, 158]
    num_pts = len(i_vals)
    for n in range(num_pts):
        i = i_vals[n]-1
        j = j_vals[n]-1
        mask_zice[j,i] = 1
        zice[j,i] = zice[j-1,i]

    # Set as land
    i_vals = [326, 182, 241, 204, 205, 233, 1099, 1099, 1130, 508, 624, 1138, 1143, 1159, 1163, 1166, 1169, 1172, 1183, 1194, 1180, 1199, 1203, 1202, 1184, 1209, 1210]
    j_vals = [117, 127, 135, 136, 137, 139, 159, 160, 179, 187, 187, 189, 195, 197, 207, 212, 216, 217, 222, 228, 229, 229, 230, 231, 232, 235, 235]
    num_pts = len(i_vals)
    for n in range(num_pts):
        i = i_vals[n]-1
        j = j_vals[n]-1
        mask_rho[j,i] = 0
        h[j,i] = 50.0

    # Remove ice shelf
    i_vals = [841, 853, 853, 691, 395]
    j_vals = [131, 135, 138, 142, 178]
    num_pts = len(i_vals)
    for n in range(num_pts):
        i = i_vals[n]-1
        j = j_vals[n]-1
        mask_zice[j,i] = 0
        zice[j,i] = 0.0

    # Remove land, h average of i+1, j-1, j+1
    i_vals = [1157, 1158]
    j_vals = [206, 209]
    num_pts = len(i_vals)
    for n in range(num_pts):
        i = i_vals[n]-1
        j = j_vals[n]-1
        mask_rho[j,i] = 1
        h[j,i] = mean(array([h[j,i+1], h[j-1,i], h[j+1,i]]))

    # Remove land, h average of i-1, i+1, j-1, j+1
    i_vals = [1159, 1160]
    j_vals = [202, 211]
    num_pts = len(i_vals)
    for n in range(num_pts):
        i = i_vals[n]-1
        j = j_vals[n]-1
        mask_rho[j,i] = 1
        h[j,i] = mean(array([h[j,i-1], h[j,i+1], h[j-1,i], h[j+1,i]]))

    # Calculate new land mask for u, v, psi grids
    mask_u = mask_rho[:,1:]*mask_rho[:,:-1]
    mask_v = mask_rho[1:,:]*mask_rho[:-1,:]
    mask_psi = mask_rho[1:,1:]*mask_rho[:-1,1:]*mask_rho[1:,:-1]*mask_rho[:-1,:-1]

    # Save new fields
    id.variables['mask_rho'][:,:] = mask_rho
    id.variables['mask_u'][:,:] = mask_u
    id.variables['mask_v'][:,:] = mask_v
    id.variables['h'][:,:] = h
    id.variables['mask_zice'][:,:] = mask_zice
    id.variables['zice'][:,:] = zice
    id.close()

if __name__ == "__main__":
    fix_isolated_pts()

