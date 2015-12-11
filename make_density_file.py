from netCDF4 import Dataset
from numpy import *
from calc_z import *
from unesco import *

# Given an ocean history or averages file with temperature and salinity data,
# calculate density fields at each timestep using the 1980 UNESCO seawater
# equation of state. Save in a new file.
# Input:
# grid_file = path to ROMS grid file
# input_file = path to ocean history/averages file
# output_file = desired path to new density file
def make_density_file (grid_file, input_file, output_file):

    # Grid parameters
    theta_s = 0.9
    theta_b = 4.0
    hc = 40
    N = 31

    # Read grid variables
    id = Dataset(grid_file, 'r')
    h = id.variables['h'][:,:]
    zice = id.variables['zice'][:,:]
    lon = id.variables['lon_rho'][:,:]
    lat = id.variables['lat_rho'][:,:]
    id.close()
    num_lon = size(lon, 1)
    num_lat = size(lon, 0)

    # Get a 3D array of z-coordinates (metres)
    z, sc_r, Cs_r = calc_z(h, zice, lon, lat, theta_s, theta_b, hc, N)
    # Pressure is approximately equal to |z|/10
    press = abs(z)/10.0

    # Set up output file
    out_id = Dataset(output_file, 'w')
    # Define dimensions
    out_id.createDimension('xi_rho', num_lon)
    out_id.createDimension('eta_rho', num_lat)
    out_id.createDimension('s_rho', N)
    out_id.createDimension('ocean_time', None)
    # Define variables
    out_id.createVariable('lon_rho', 'f8', ('eta_rho', 'xi_rho'))
    out_id.variables['lon_rho'][:,:] = lon
    out_id.createVariable('lat_rho', 'f8', ('eta_rho', 'xi_rho'))
    out_id.variables['lat_rho'][:,:] = lat
    out_id.createVariable('sc_r', 'f8', ('s_rho'))
    out_id.variables['sc_r'].long_name = 'S-coordinate at rho-points'
    out_id.variables['sc_r'][:] = sc_r
    out_id.createVariable('ocean_time', 'f8', ('ocean_time'))
    out_id.variables['ocean_time'].units = 'seconds'
    out_id.createVariable('rho', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'))
    out_id.variables['rho'].long_name = 'density'
    out_id.variables['rho'].units = 'kg/m^3'

    # Read time values from input file
    in_id = Dataset(input_file, 'r')
    time = in_id.variables['ocean_time'][:]

    # Process each timestep individually to conserve memory
    for t in range(size(time)):
        print 'Processing timestep '+str(t+1)+' of '+str(size(time))
        # Set a new time value in the output file
        out_id.variables['ocean_time'][t] = time[t]
        # Read temperature and salinity (convert to float128 to prevent
        # overflow in UNESCO calculations)
        temp = ma.asarray(in_id.variables['temp'][t,:,:,:], dtype=float128)
        salt = ma.asarray(in_id.variables['salt'][t,:,:,:], dtype=float128)
        # Magic happens here
        rho = unesco(temp, salt, press)
        # Save the results for this timestep
        out_id.variables['rho'][t,:,:,:] = rho

    in_id.close()
    out_id.close()


# Command-line interface
if __name__ == "__main__":

    grid_file = raw_input("Path to grid file: ")
    input_file = raw_input("Path to ocean history/averages file: ")
    output_file = raw_input("Desired path to new density file: ")
    make_density_file(grid_file, input_file, output_file)
