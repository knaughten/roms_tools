from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from rotate_vector_roms import *

# Calculates zonal transport through each grid cell in the Drake Passage
# and plots a timeseries of the result.
# Input:
# grid_path = path to ROMS grid file
# file_path = path to ROMS ocean history or averages file
def dpt_timeseries (grid_path, file_path):

    # Radius of the Earth in m
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Northern boundary of ROMS grid
    nbdry_val = -30

    # Bounds on Drake Passage; edit for new grids
    # i-index of single north-south line to plot (representing a zonal slice);
    # it doesn't really matter which slice of the Drake Passage this is, due
    # to volume conservation
    i_DP = 1179
    # j-indices of the southern tip of South America (j_min) and the northern
    # tip of the Antarctic Peninsula (j_max); make sure these are far enough
    # north/south to be land points, but not so far that they pass through the
    # land and become ocean again (eg Weddell Sea)
    j_min = 229
    j_max = 298

    print 'Reading grid'
    # Read angle from the grid file
    grid_id = Dataset(grid_path, 'r')
    angle = grid_id.variables['angle'][:-15,:]
    grid_id.close()
    # Read other grid variables
    id = Dataset(file_path, 'r')
    h = id.variables['h'][:-15,:-3]
    zice = id.variables['zice'][:-15,:-3]
    lon = id.variables['lon_rho'][:-15,:-3]
    lat = id.variables['lat_rho'][:-15,:-3]
    mask = id.variables['mask_rho'][:-15,:-3]

    # Interpolate latitude to the edges of each cell
    s_bdry = lat[0,:]
    middle_lat = 0.5*(lat[0:-1,:] + lat[1:,:])
    n_bdry = lat[-1,:]*0 + nbdry_val
    lat_edges = ma.concatenate((s_bdry[None,:], middle_lat, n_bdry[None,:]))
    # Subtract to get the change in latitude over each cell
    dlat = lat_edges[1:,:] - lat_edges[0:-1,:]
    # Convert from spherical to Cartesian coordinates
    # dy = r*dlat where dlat is converted to radians
    dy = r*dlat*pi/180.0

    # Read time values and convert from seconds to years
    time = id.variables['ocean_time'][:]/(365*24*60*60)

    transport = []
    # Calculate transport one timestep at a time
    for t in range(size(time)):

        print 'Processing timestep ' + str(t+1) + ' of '+str(size(time))
        # Calculate water column thickness
        zeta = id.variables['zeta'][t,:-15,:-3]
        wct = h + zice + zeta
        # Calculate dy_wct and mask with land mask
        dy_wct = ma.masked_where(mask==0, dy*wct)
        # Trim to Drake Passage bounds
        dy_wct_DP = dy_wct[j_min:j_max,i_DP]
        # Rotate velocities into lat-lon space
        ubar = id.variables['ubar'][t,:-15,:]
        vbar = id.variables['vbar'][t,:-15,:]
        ubar_lonlat, vbar_lonlat = rotate_vector_roms(ubar, vbar, angle)
        # Throw away the overlapping periodic boundary
        ubar_lonlat = ubar_lonlat[:,:-3]
        # Trim to Drake Passage bounds
        ubar_DP = ubar_lonlat[j_min:j_max,i_DP]
        # Integrate transport and convert to Sv
        transport.append(sum(ubar_DP*dy_wct_DP)*1e-6)

    id.close()

    # Plot
    # Bounds are set to +/- 16 Sv, adjust as needed
    figure()
    plot(time, transport)
    xlabel('Years')
    ylabel('Drake Passage Transport (Sv)')
    grid(True)
    show(block=False)
    #savefig('dpt_timeseries.png')


# Command-line interface
if __name__ == "__main__":

    grid_path = raw_input('Enter path to ROMS grid file: ')
    file_path = raw_input('Enter path to ocean history/averages file: ')
    dpt_timeseries(grid_path, file_path)
    


