from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from rotate_vector_roms import *
from interp_lon_roms import interp_lon_helper

# Calculate and plot timeseries of the Drake Passage transport during a
# ROMS simulation.
# Input:
# grid_path = path to ROMS grid file
# file_path = path to ROMS history/averages file
# log_path = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
def timeseries_dpt (grid_path, file_path, log_path):

    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Longitude of Drake Passage zonal slice (convert to ROMS bounds 0-360)
    lon0 = -67 + 360
    # Latitude bounds on Drake Passage zonal slice
    lat_min = -68
    lat_max = -54.5

    time = []
    dpt = []
    # Check if the log file exists
    if exists(log_path):
        print 'Reading previously calculated values'
        f = open(log_path, 'r')
        # Skip first line (header for time array)
        f.readline()
        for line in f:
            try:
                time.append(float(line))
            except(ValueError):
                # Reached the header for the next variable
                break
        for line in f:
            dpt.append(float(line))
        f.close()

    print 'Reading grid'
    id = Dataset(grid_path, 'r')
    h = id.variables['h'][:-15,1:-1]
    zice = id.variables['zice'][:-15,1:-1]
    lon = id.variables['lon_rho'][:-15,1:-1]
    lat = id.variables['lat_rho'][:-15,1:-1]
    # Keep the overlapping periodic boundary on "angle" for now
    angle = id.variables['angle'][:-15,:]
    id.close()

    print 'Reading data'
    id = Dataset(file_path, 'r')
    # Read time values and convert from seconds to years
    new_time = id.variables['ocean_time'][:]/(60*60*24*365.25)
    num_time = size(new_time)
    # Concatenate with time values from log file
    for t in range(num_time):
        time.append(new_time[t])
    # Calculate time-dependent water column thickness: h + zice + zeta
    zeta = id.variables['zeta'][:,:-15,1:-1]
    wct = tile(h, (num_time,1,1)) + tile(zice, (num_time,1,1)) + zeta
    # Read barotropic velocities in x-y space
    ubar_xy = id.variables['ubar'][:,:-15,:]
    vbar_xy = id.variables['vbar'][:,:-15,:]
    id.close()

    print 'Rotating velocity vector'
    ubar = ma.empty([num_time, size(lon,0), size(lon,1)])
    # Rotate one time index at a time
    for t in range(num_time):        
        ubar_tmp, vbar_tmp = rotate_vector_roms(ubar_xy[t,:,:], vbar_xy[t,:,:], angle)
        # Throw away the overlapping periodic boundary before saving
        ubar[t,:,:] = ubar_tmp[:,1:-1]

    print 'Extracting zonal slice through Drake Passage'    #
    num_lat = size(lat,0)
    # Set up arrays for zonal slices of ubar, water column thickness, latitude
    ubar_DP = ma.empty([num_time, num_lat])
    wct_DP = ma.empty([num_time, num_lat])
    lat_DP = empty([num_lat])
    # Loop over longitudes
    for j in range(num_lat):        
        lon_tmp = lon[j,:]
        # Find indices and coefficients to interpolate to lon0
        ie, iw, coeffe, coeffw = interp_lon_helper(lon_tmp, lon0)
        # Use these to interpolate all 3 variables we care about
        ubar_DP[:,j] = coeffe*ubar[:,j,ie] + coeffw*ubar[:,j,iw]
        wct_DP[:,j] = coeffe*wct[:,j,ie] + coeffw*wct[:,j,iw]
        lat_DP[j] = coeffe*lat[j,ie] + coeffw*lat[j,iw]
    # Find indices for latitude bounds
    jS = nonzero(lat_DP > lat_min)[0][0]
    jN = nonzero(lat_DP > lat_max)[0][0]
    # Trim everything to these bounds
    ubar_DP = ubar_DP[:,jS:jN]
    wct_DP = wct_DP[:,jS:jN]
    lat_DP = lat_DP[jS:jN]
    # Calculate dy
    # First calculate latitude on edges of each cell
    middle_lat = 0.5*(lat_DP[:-1] + lat_DP[1:])
    s_bdry = 2*lat_DP[0] - middle_lat[0]
    n_bdry = 2*lat_DP[-1] - middle_lat[-1]
    lat_edges = zeros(size(lat_DP)+1)
    lat_edges[0] = s_bdry
    lat_edges[1:-1] = middle_lat
    lat_edges[-1] = n_bdry
    # Now calculate difference in latitude across each cell
    dlat_DP = lat_edges[1:] - lat_edges[:-1]
    # Convert to Cartesian space for dy in metres
    dy_DP = r*dlat_DP*deg2rad

    for t in range(num_time):
        # Integrate ubar*wct*dy and convert to Sv
        dpt.append(sum(ubar_DP[t,:]*wct_DP[t,:]*dy_DP)*1e-6)

    print 'Plotting'
    clf()
    plot(time, dpt)
    xlabel('Years')
    ylabel('Drake Passage Transport (Sv)')
    grid(True)
    savefig('drakepsgtrans.png')

    print 'Saving results to log file'
    f = open(log_path, 'w')
    f.write('Time (years):\n')
    for elm in time:
        f.write(str(elm) + '\n')
    f.write('Drake Passage Transport (Sv):\n')
    for elm in dpt:
        f.write(str(elm) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    grid_path = raw_input("Path to ROMS grid file: ")
    file_path = raw_input("Path to ROMS history/averages file: ")
    log_path = raw_input("Path to logfile to save values and/or read previously calculated values: ")
    timeseries_dpt(grid_path, file_path, log_path)

    
    
    
