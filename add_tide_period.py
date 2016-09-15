from netCDF4 import Dataset
from numpy import *

# Given a ROMS tide file with tidal potential amplitude and phase angle
# (created using Kate's potential tide scripts), add the tidal period in seconds
# for each component.
def add_tide_period ():

    # Tide file to add periods to
    # Order is q1, o1, p1, k1, n2, m2, s2, k2
    tide_file = '../ROMS-CICE-MCT/data/pot_tides.nc'
    # Periods in seconds
    # Order is m2, s2, n2, k2, k1, o1, p1, q1, mf, mm
    period = array([44714.165191868, 43200.0012869521, 86164.0770050671, 92949.6357005365, 45570.0535117177, 86637.1997716528, 43082.0503185947, 96726.0857029666, 2380715.86358729, 1180295.54554976])
    # Put them in the right order
    period_reorder = array([period[7], period[5], period[6], period[4], period[2], period[0], period[1], period[3]])

    id = Dataset(tide_file, 'a')
    id.createVariable('tide_period', 'f8', ('tide_period'))
    id.variables['tide_period'].long_name = 'tide angular period'
    id.variables['tide_period'].units = 'seconds'
    id.variables['tide_period'][:] = period_reorder
    id.close()


# Command-line interface
if __name__ == "__main__":

    add_tide_period()
    
