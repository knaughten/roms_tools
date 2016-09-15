from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib import rcParams

# Make a map of unexplained percent error in annually averaged simulated basal
# mass loss from each ice shelf that is over 5,000 km^2 in Rignot et al., 2013.
# Input:
# grid_path = path to ROMS grid file
# log_path = path to log file created by timeseries_massloss.py
# save = optional boolean to save the figure to a file, rather than displaying
#        it on the screen
# fig_name = if save=True, path to the desired filename for the figure
def massloss_map (grid_path, log_path, save=False, fig_name=None):

    # Limits on longitude and latitude for each ice shelf
    # These depend on the source geometry, in this case RTopo 1.05
    # Note there is one extra index at the end of each array; this is because
    # the Ross region crosses the line 180W and therefore is split into two
    lon_min = [-62.67, -65.5, -79.17, -85, -104.17, -102.5, -108.33, -114.5, -135.67, -149.17, -155, 144, 115, 94.17, 80.83, 65, 33.83, 19, 12.9, 9.33, -10.05, -28.33, -180, 158.33]
    lon_max = [-59.33, -60, -66.67, -28.33, -88.83, -99.17, -103.33, -111.5, -114.33, -140, -145, 146.62, 123.33, 102.5, 89.17, 75, 37.67, 33.33, 16.17, 12.88, 7.6, -10.33, -146.67, 180]
    lat_min = [-73.03, -69.35, -74.17, -83.5, -73.28, -75.5, -75.5, -75.33, -74.9, -76.42, -78, -67.83, -67.17, -66.67, -67.83, -73.67, -69.83, -71.67, -70.5, -70.75, -71.83, -76.33, -85, -84.5]
    lat_max = [-69.37, -66.13, -69.5, -74.67, -71.67, -74.17, -74.67, -73.67, -73, -75.17, -76.41, -66.67, -66.5, -64.83, -66.17, -68.33, -68.67, -68.33, -69.33, -69.83, -69.33, -71.5, -77.77, -77]
    # Observed mass loss (Rignot 2013) and uncertainty for each ice shelf, in Gt/y
    obs_massloss = [1.4, 20.7, 135.4, 155.4, 51.8, 101.2, 97.5, 45.2, 144.9, 4.2, 18.2, 7.9, 90.6, 72.6, 27.2, 35.5, -2, 21.6, 6.3, 3.9, 26.8, 9.7, 47.7]
    obs_massloss_error = [14, 67, 40, 45, 19, 8, 7, 4, 14, 2, 3, 3, 8, 15, 10, 23, 3, 18, 2, 2, 14, 16, 34]

    # Degrees to radians conversion factor
    deg2rad = pi/180
    # Northern boundary 63S for plot
    nbdry = -63+90
    # Centre of missing circle in grid
    lon_c = 50
    lat_c = -83
    # Radius of missing circle (play around with this until it works)
    radius = 10.1
    # Minimum zice
    min_zice = -10

    # Read log file
    time = []
    f = open(log_path, 'r')
    # Skip the first line (header for time array)
    f.readline()
    for line in f:
        try:
            time.append(float(line))
        except(ValueError):
            # Reached the header for the next variable
            break
    # Skip the values for the entire continent
    for line in f:
        try:
            pass
        except(ValueError):
            break
    # Set up array for mass loss values at each ice shelf
    massloss_ts = empty([len(obs_massloss), len(time)])
    index = 0
    # Loop over ice shelves
    while index < len(obs_massloss):
        t = 0
        for line in f:
            try:
                massloss_ts[index, t] = float(line)
                t += 1
            except(ValueError):
                # Reached the header for the next variable
                break
        index +=1

    # Find the time indices we care about: last year of simulation
    time = array(time)
    min_t = nonzero(time >= time[-1]-1)[0][0]
    max_t = size(time)
    # Find the average mass loss for each ice shelf over this last year
    massloss = empty(len(obs_massloss))
    for index in range(len(obs_massloss)):
        massloss[index] = mean(massloss_ts[index, min_t:max_t])

    # Read the grid
    id = Dataset(grid_path, 'r')
    lon = id.variables['lon_rho'][:-15,:-2]
    lat = id.variables['lat_rho'][:-15,:-2]
    mask_rho = id.variables['mask_rho'][:-15,:-2]
    mask_zice = id.variables['mask_zice'][:-15,:-2]
    zice = id.variables['zice'][:-15,:-2]
    id.close()
    
    # Make sure longitude goes from -180 to 180, not 0 to 360
    index = lon > 180
    lon[index] = lon[index] - 360

    # Get land/zice mask
    open_ocn = copy(mask_rho)
    open_ocn[mask_zice==1] = 0
    land_zice = ma.masked_where(open_ocn==1, open_ocn)

    # Initialise a field of ice shelf mass loss unexplained percent error
    error = ma.empty(shape(lon))
    error[:,:] = ma.masked

    # Loop over ice shelves
    for index in range(len(obs_massloss)):
        # Find the range of observations
        massloss_low = obs_massloss[index] - obs_massloss_error[index]
        massloss_high = obs_massloss[index] + obs_massloss_error[index]
        # Find the unexplained percent error in mass loss
        if massloss[index] < massloss_low:
            # Simulated mass loss too low
            error_tmp = (massloss[index] - massloss_low)/massloss_low*100
        elif massloss[index] > massloss_high:
            # Simulated mass loss too high
            error_tmp = (massloss[index] - massloss_high)/massloss_high*100
        else:
            # Simulated mass loss within observational error estimates
            error_tmp = 0
        # Modify error field for this region
        if index == len(obs_massloss)-1:
            # Ross region is split into two
            region = (lon >= lon_min[index])*(lon <= lon_max[index])*(lat >= lat_min[index])*(lat <= lat_max[index])*(mask_zice == 1) + (lon >= lon_min[index+1])*(lon <= lon_max[index+1])*(lat >= lat_min[index+1])*(lat <= lat_max[index+1])*(mask_zice == 1)
        else:
            region = (lon >= lon_min[index])*(lon <= lon_max[index])*(lat >= lat_min[index])*(lat <= lat_max[index])*(mask_zice == 1)
        error[region] = error_tmp

    # Edit zice so tiny ice shelves won't be contoured
    zice[error.mask] = 0.0

    # Convert grid to spherical coordinates
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)
    # Find centre in spherical coordinates
    x_c = -(lat_c+90)*cos(lon_c*deg2rad+pi/2)
    y_c = (lat_c+90)*sin(lon_c*deg2rad+pi/2)
    # Build a regular x-y grid and select the missing circle
    x_reg, y_reg = meshgrid(linspace(-nbdry, nbdry, num=1000), linspace(-nbdry, nbdry, num=1000))
    land_circle = zeros(shape(x_reg))
    land_circle = ma.masked_where(sqrt((x_reg-x_c)**2 + (y_reg-y_c)**2) > radius, land_circle)

    # Determine bounds on colour scale
    max_val = amax(abs(error))
    lev = linspace(-max_val, max_val, num=40)
    lev = linspace(-100, 100, num=40)
    # Space ticks on colorbar 25% apart
    max_tick = floor(max_val/25)*25

    # Set up plot
    fig = figure(figsize=(16,12))
    ax = fig.add_subplot(1,1,1, aspect='equal')
    fig.patch.set_facecolor('white')
    # First shade land and zice in grey (include zice so there are no white
    # patches near the grounding line where contours meet)
    contourf(x, y, land_zice, 1, colors=(('0.6', '0.6', '0.6')))
    # Fill in the missing cicle
    contourf(x_reg, y_reg, land_circle, 1, colors=(('0.6', '0.6', '0.6')))
    # Now shade the error in mass loss
    contourf(x, y, error, lev, cmap='RdBu_r', extend='both')    
    cbar = colorbar(ticks=arange(-max_tick, max_tick+25, 25))
    cbar.ax.tick_params(labelsize=20)
    # Add a black contour for the ice shelf front
    rcParams['contour.negative_linestyle'] = 'solid'
    contour(x, y, zice, levels=[min_zice], colors=('black'))
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    title('Bias in Ice Shelf Mass Loss (%)', fontsize=30)
    axis('off')

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    grid_path = raw_input("Path to grid file: ")
    log_path = raw_input("Path to mass loss logfile: ")
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    # Make the plot
    massloss_map(grid_path, log_path, save, fig_name)
            
    
    
    
