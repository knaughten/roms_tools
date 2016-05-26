from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *

# Calculates zonal transport through each grid cell in the Drake Passage
# and plots a timeseries of the result.
# Input:
# file_path = path to ROMS ocean history or averages file
def dpt_timeseries (file_path):

    # Radius of the Earth in m
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Northern boundary of ROMS grid
    nbdry_val = -30

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
    # Read grid variables
    id = Dataset(file_path, 'r')
    h = id.variables['h'][:,:]
    zice = id.variables['zice'][:,:]
    lon = id.variables['lon_rho'][:,:]
    lat = id.variables['lat_rho'][:,:]
    mask = id.variables['mask_rho'][:,:]


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
    # Calculate water column thickness
    wct = h + zice

    # Calculate dy_wct and mask with land mask
    dy_wct = ma.masked_where(mask==0, dy*wct)
    # Trim to Drake Passage bounds
    dy_wct_DP = dy_wct[j_min:j_max,i_DP]
    lat_DP = lat[j_min:j_max,i_DP]

    # Read time values and convert from seconds to years
    time = id.variables['ocean_time'][:]/(365*24*60*60)

    transport = []
    # Calculate transport one timestep at a time
    for t in range(size(time)):

        print 'Processing timestep ' + str(t+1) + ' of '+str(size(time))
        # Read ubar and interpolate onto the rho-grid
        ubar = id.variables['ubar'][t,:,:]
        w_bdry_ubar = 0.5*(ubar[:,0] + ubar[:,-1])
        middle_ubar = 0.5*(ubar[:,0:-1] + ubar[:,1:])
        e_bdry_ubar = w_bdry_ubar[:]
        ubar_rho = ma.concatenate((w_bdry_ubar[:,None], middle_ubar, e_bdry_ubar[:,None]), axis=1)
        # Trim to Drake Passage bounds
        ubar_rho_DP = ubar_rho[j_min:j_max,i_DP]
        # Integrate transport and convert to Sv
        transport.append(sum(ubar_rho_DP*dy_wct_DP)*1e-6)

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

    file_path = raw_input('Enter path to ocean history/averages file: ')
    dpt_timeseries(file_path)
    


