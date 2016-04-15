from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *

# Calculate and plot timeseries of basal mass loss from major ice shelves
# during a ROMS simulation.
# Input:
# file_path = path to ocean history/averages file
# log_path = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
def massloss (file_path, log_path):

    # Titles and figure names for each ice shelf
    names = ['Amery Ice Shelf', 'Ross Ice Shelf', 'Getz Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Abbot Ice Shelf', 'George VI Ice Shelf', 'Larsen C Ice Shelf', 'Ronne-Filchner Ice Shelf', 'Brunt and Riiser-Larsen Ice Shelves', 'Fimbul and Jelbart Ice Shelves']
    fig_names = ['amery.png', 'ross.png', 'getz.png', 'pig.png', 'abbot.png', 'george_vi.png', 'larsen_c.png', 'ronne_filchner.png', 'brunt_riiser_larsen.png', 'fimbul_jelbart.png']
    # Limits on i- and j- coordinates for each ice shelf; this will vary
    # depending on the grid
    # Note there is one extra index at the end of each array; this is because
    # the Fimbul-Jelbart region crosses the periodic boundary and therefore
    # is split into two
    i_min = [250, 700, 905, 1015, 1000, 1100, 1165, 1060, 1280, 1375, 1]
    i_max = [350, 872, 975, 1030, 1090, 1155, 1190, 1240, 1369, 1443, 12]
    j_min = [1,   20,  150, 140,  160,  150,  190,  1,    65,   80,   100]
    j_max = [125, 123, 175, 160,  185,  200,  220,  135,  116,  150,  120]
    # Density of ice in kg/m^3
    rho_ice = 916

    old_time = []
    # Check if the log file exists
    if exists(log_path):
        print 'Reading previously calculated values'
        f = open(log_path, 'r')
        # Skip the first line (header for time array)
        f.readline()
        for line in f:
            try:
                old_time.append(float(line))
            except(ValueError):
                # Reached the header for the next variable
                break
        # Set up array for mass loss values at each ice shelf
        old_massloss = empty([len(names), len(old_time)])
        index = 0
        # Loop over ice shelves
        while index < len(names):
            t = 0
            for line in f:
                try:
                    old_massloss[index, t] = float(line)
                    t += 1
                except(ValueError):
                    # Reached the header for the next ice shelf
                    break
            index += 1

    # Calculate dA (masked with ice shelf mask) and i and j coordinates
    print 'Analysing grid'
    dA, i, j = calc_grid(file_path)

    # Read time data and convert from seconds to years
    id = Dataset(file_path, 'r')
    new_time = id.variables['ocean_time'][:]/(365*25*60*60)
    if exists(log_path):
        # Concatenate with time values from log file
        start_t = len(old_time)
        time = old_time
        for t in range(size(new_time)):
            time.append(new_time[t])
        time = array(time)
    else:
        start_t = 0
        time = new_time
    # Set up array of mass loss values
    massloss = empty([len(names), size(time)])
    if exists(log_path):
        # Fill first start_t timesteps with existing values
        massloss[:,0:start_t] = old_massloss[:,:]

    # Process each timestep separately to prevent memory overflow
    for t in range(start_t, size(time)):
        print 'Processing timestep ' + str(t-start_t+1) + ' of ' + str(size(time)-start_t)

        # Read ice shelf melt rate, converting to float128 to prevent overflow
        # during integration
        ismr = ma.asarray(id.variables['m'][t-start_t,:,:], dtype=float128)
        # Convert from m/s to m/y
        ismr = ismr*365.25*24*60*60

        # Loop over ice shelves
        for index in range(len(names)):

            # Mask dA for the current ice shelf (note dA is already masked
            # with the global ice shelf mask, so just restrict the i and j
            # to isolate the given ice shelf)
            if index == len(names)-1:
                # Fimbul-Jelbart region is split into two
                dA_tmp = ma.masked_where(((i < i_min[index]) + (i > i_max[index]) + (j < j_min[index]) + (j > j_max[index]))*((i < i_min[index+1]) + (i > i_max[index+1]) + (j < j_min[index+1]) + (j > j_max[index+1])), dA)
            else:
                dA_tmp = ma.masked_where((i < i_min[index]) + (i > i_max[index]) + (j < j_min[index]) + (j > j_max[index]), dA)

            # Integrate ice shelf melt rate over area to get volume loss
            volumeloss = sum(ismr*dA_tmp)
            # Convert to mass loss in Gt/y
            massloss[index, t] = 1e-12*rho_ice*volumeloss

    id.close()

    # Plot each timeseries
    print 'Plotting'
    for index in range(len(names)):
        clf()
        plot(time, massloss[index,:])
        xlabel('Years')
        ylabel('Basal Mass Loss (Gt/y)')
        title(names[index])
        savefig(fig_names[index])
        
    print 'Saving results to log file'
    f = open(log_path, 'w')
    f.write('Time (years):\n')
    for t in range(size(time)):
        f.write(str(time[t]) + '\n')
    for index in range(len(names)):
        f.write(names[index] + '\n')
        for t in range(size(time)):
            f.write(str(massloss[index, t]) + '\n')
    f.close()


# Given the path to a ROMS grid file, calculate differential of area and
# i- and j-indices.
# Input: file_path = string containing path to ROMS history/averages file
# Output:
# dA = differential of area on the 2D rho-grid, masked with zice
# i,j = i- and j- coordinates on the rho-grid, starting at 1
def calc_grid (file_path):

    # Radius of the Earth in m
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Northern boundary of ROMS grid
    nbdry_val = -30

    # Read grid variables
    id = Dataset(file_path, 'r')
    lon = id.variables['lon_rho'][:,:]
    lat = id.variables['lat_rho'][:,:]
    zice = id.variables['zice'][:,:]
    id.close()
    # Save dimensions
    num_lat = size(lon, 0)
    num_lon = size(lon, 1)

    # Add or subtract 360 from longitude values which wrap around
    # so that longitude increases monotonically from west to east
    i = tile(arange(1, num_lon+1), (num_lat, 1))
    index1 = nonzero((i > 1200)*(lon < 100))
    lon[index1] = lon[index1] + 360
    index2 = nonzero((i < 200)*(lon > 300))
    lon[index2] = lon[index2] - 360

    # Interpolate to get longitude at the edges of each cell
    w_bdry = 0.5*(lon[:,0] + lon[:,-1] - 360)
    middle_lon = 0.5*(lon[:,0:-1] + lon[:,1:])
    e_bdry = 0.5*(lon[:,0] + 360 + lon[:,-1])
    lon_edges = ma.concatenate((w_bdry[:,None], middle_lon, e_bdry[:,None]), axis=1)
    # Subtract to get the change in longitude over each cell
    dlon = abs(lon_edges[:,1:] - lon_edges[:,0:-1])

    # Similarly for latitude
    s_bdry = lat[0,:]
    middle_lat = 0.5*(lat[0:-1,:] + lat[1:,:])
    n_bdry = lat[-1,:]*0 + nbdry_val
    lat_edges = ma.concatenate((s_bdry[None,:], middle_lat, n_bdry[None,:]))
    dlat = lat_edges[1:,:] - lat_edges[0:-1,:]

    # Also calculate j-coordinates
    j = transpose(tile(arange(1, num_lat+1), (num_lon, 1)))

    # Convert from spherical to Cartesian coordinates
    # dy = r*dlat where dlat is converted to radians
    dy = r*dlat*deg2rad    
    # dx = r*cos(lat)*dlon where lat and dlon are converted to radians
    dx = r*cos(lat*deg2rad)*dlon*deg2rad

    # Calculate dA and mask with zice
    dA = ma.masked_where(zice==0, dx*dy)

    return dA, i, j


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input('Enter path to ocean history/averages file: ')
    log_path = raw_input('Enter path to log file to save values and/or read previously calculated values: ')

    massloss(file_path, log_path)

