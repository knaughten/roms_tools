from netCDF4 import Dataset
from numpy import *
from calc_z import *

def calc_rx1 (grid_path, theta_s, theta_b, hc, N, Vstretching, out_file):

    id = Dataset(grid_path, 'r')
    lon_2d = id.variables['lon_rho'][:,:]
    lat_2d = id.variables['lat_rho'][:,:]
    h = id.variables['h'][:,:]
    zice = id.variables['zice'][:,:]
    mask_rho = id.variables['mask_rho'][:,:]
    id.close()

    z, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N, None, Vstretching)
    for k in range(N):
        tmp = z[k,:,:]
        tmp[mask_rho==0] = NaN
        z[k,:,:] = tmp

    rx1_3d_i = abs((z[1:,1:,1:]-z[1:,1:,:-1]+z[:-1,1:,1:]-z[:-1,1:,:-1])/(z[1:,1:,1:]+z[1:,1:,:-1]-z[:-1,1:,1:]-z[:-1,1:,:-1]))
    rx1_3d_j = abs((z[1:,1:,1:]-z[1:,:-1,1:]+z[:-1,1:,1:]-z[:-1,:-1,1:])/(z[1:,1:,1:]+z[1:,:-1,1:]-z[:-1,1:,1:]-z[:-1,:-1,1:]))
    rx1_3d = maximum(rx1_3d_i, rx1_3d_j)
    rx1_tmp = amax(rx1_3d, axis=0)

    rx1 = zeros(shape(lon_2d))
    rx1[1:,1:] = rx1_tmp
    rx1[0,:] = rx1[1,:]
    rx1[:,0] = rx1[:,1]
    rx1 = ma.masked_where(isnan(rx1), rx1)

    id = Dataset(out_file, 'w')
    id.createDimension('xi_rho', size(lon_2d,1))
    id.createDimension('eta_rho', size(lon_2d,0))
    id.createVariable('rx1', 'f8', ('eta_rho', 'xi_rho'))
    id.variables['rx1'][:,:] = rx1
    id.close()

    
if __name__ == "__main__":

    grid_path = raw_input("Path to grid file: ")
    theta_s = float(raw_input("theta_s: "))
    theta_b = float(raw_input("theta_b: "))
    hc = float(raw_input("hc: "))
    N = int(raw_input("N: "))
    Vstretching = int(raw_input("Vstretching (2 or 4): "))
    out_file = raw_input("Path to output file: ")

    calc_rx1 (grid_path, theta_s, theta_b, hc, N, Vstretching, out_file)
