from netCDF4 import Dataset
from numpy import *

# Given a ROMS atmospheric boundary condition file with one year of monthly
# data, convert it to a file which can be used repeatedly year after year.

# Main routine
def run (filename):

    id = Dataset(filename, 'a')
    # Set time with respect to day 0 of the simulation, one value per month,
    # centered in the middle of each month
    time = (arange(12)+0.5)*365.25/12.0
    id.variables['time'][:] = time
    # Set cycle length of 1 year
    id.variables['time'].cycle_length = 365.25
    id.close()

# User parameters
if __name__ == "__main__":

    run('../ROMS-CICE-MCT/data/ERA_Interim/monthly/AN_1995_monthly.nc')
    run('../ROMS-CICE-MCT/data/ERA_Interim/monthly/FC_1995_monthly.nc')
