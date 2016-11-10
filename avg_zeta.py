from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from cartesian_grid_2d import *

# Calculate the average sea surface height (zeta) at each timestep of the given
# ocean history file.
def avg_zeta (file_path):

    # Read time and grid variables
    file = Dataset(file_path, 'r')
    time = file.variables['ocean_time'][:]
    # Convert time from seconds to years
    time = time/(365*24*60*60)
    lon = file.variables['lon_rho'][:-15,1:-1]
    lat = file.variables['lat_rho'][:-15,1:-1]
    mask = file.variables['mask_rho'][:-15,1:-1]
    avg_zeta = []

    # Calculate dx and dy in another script
    dx, dy = cartesian_grid_2d(lon, lat)

    # Calculate dA and mask with land mask
    dA = ma.masked_where(mask==0, dx*dy)

    for l in range(size(time)):
        print 'Processing timestep ' + str(l+1) + ' of ' + str(size(time))
        # Read zeta at this timestep
        zeta = file.variables['zeta'][l,:-15,1:-1]
        # Calculate area-weighted average
        avg_zeta.append(sum(zeta*dA)/sum(dA))

    file.close()

    # Plot results
    clf()
    plot(time, avg_zeta)
    xlabel('Years')
    ylabel('Average sea surface height (m)')
    show()


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to ocean history file: ")
    avg_zeta(file_path)

    
        
