from netCDF4 import Dataset
from numpy import *

# Given a ROMS lateral boundary condition file with one year of monthly data,
# convert it to a file which can be used repeatedly year after year.

# Main routine
def run (filename):

    id = Dataset(filename, 'a')
    # Set time with respect to day 0 of the simulation, one value per month,
    # centered in the middle of each month
    time = (arange(12)+0.5)*365.25/12.0
    id.variables['ocean_time'][:] = time
    # Set cycle length of 1 year
    id.variables['ocean_time'].cycle_length = 365.25
    id.close()

# User parameters
if __name__ == "__main__":

    filename = '../metroms_iceshelf/data/ecco2_cube92_lbc_1995_rep.nc'
    run(filename)
