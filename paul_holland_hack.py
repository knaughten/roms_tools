from netCDF4 import Dataset
from numpy import *
from cartesian_grid_2d import *

# Paul Holland prevents spurious deep convection by applying an extra 1500 Gt/y
# of freshwater flux spread evenly south of 60S. Using the ROMS grid, calculate
# what the rate of input should be in kg/m^2/s. Print the value to the screen.
# Input: grid_file = path to ROMS grid file
def paul_holland_hack (grid_file):

    total_fw = 1500 # Gt/y
    nbdry = -60 # Apply the freshwater evenly south of here
    sec_per_year = 365.25*24*60*60

    # Read grid and masks, making sure to get rid of the overlapping periodic
    # boundary cells that are double-counted
    id = Dataset(grid_file, 'r')
    lat = id.variables['lat_rho'][:,1:-1]
    lon = id.variables['lon_rho'][:,1:-1]
    mask_zice = id.variables['mask_zice'][:,1:-1]
    mask_rho = id.variables['mask_rho'][:,1:-1]
    id.close()
    # Mask out land and ice shelves
    mask = mask_rho - mask_zice

    # Get differentials
    dx, dy = cartesian_grid_2d(lon, lat)
    # Open ocean cells
    ocn_flag = mask == 1
    # Cells south of 60S
    loc_flag = lat < nbdry
    # Total area of all open ocean cells
    total_area = sum(dx*dy*ocn_flag)
    print 'Total area = ' + str(total_area) + ' m^2'
    # Total area of open ocean cells south of 60S
    target_area = sum(dx*dy*ocn_flag*loc_flag)
    print 'Area south of 60S = ' + str(target_area) + ' m^2'
    # Multiply by 1e12 to convert from Gt/y to kg/y
    # Divide by sec_per_year to convert from kg/y to kg/s
    # Divide by target area to get kg/m^2/s
    fw_flux = total_fw*1e12/target_area/sec_per_year
    print 'Freshwater flux to add = ' + str(fw_flux) + 'kg/m^2/s'


# Command-line interface
if __name__ == "__main__":

    grid_file = raw_input("Path to ROMS grid file: ")
    paul_holland_hack(grid_file)
