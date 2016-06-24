from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from cartesian_grid_2d import *

# Make a map of the magnitude of the seasonal cycle (over the last year of
# simulation) in area-averaged melt rates for each ice shelf that is over 
# 5,000 km^2 in Rignot et al., 2013.
# Input:
# grid_path = path to ROMS grid file
# log_path = path to logfile created using massloss.py
# save = optional boolean to save the figure to a file, rather than display it
#        on screen
# fig_name = if save=True, filename for figure
def ismr_seasonal_cycle (grid_path, log_path, save=False, fig_name=None):

    # Limits on longitude and latitude for each ice shelf
    # These depend on the source geometry, in this case RTopo 1.05
    # Note there is one extra index at the end of each array; this is because
    # the Ross region crosses the line 180W and therefore is split into two
    lon_min = [-62.67, -65.5, -79.17, -85, -104.17, -102.5, -108.33, -114.5, -135.67, -149.17, -155, 144, 115, 94.17, 80.83, 65, 33.83, 19, 12.9, 9.33, -10.05, -28.33, -180, 158.33]
    lon_max = [-59.33, -60, -66.67, -28.33, -88.83, -99.17, -103.33, -111.5, -114.33, -140, -145, 146.62, 123.33, 102.5, 89.17, 75, 37.67, 33.33, 16.17, 12.88, 7.6, -10.33, -146.67, 180]
    lat_min = [-73.03, -69.35, -74.17, -83.5, -73.28, -75.5, -75.5, -75.33, -74.9, -76.42, -78, -67.83, -67.17, -66.67, -67.83, -73.67, -69.83, -71.67, -70.5, -70.75, -71.83, -76.33, -85, -84.5]
    lat_max = [-69.37, -66.13, -69.5, -74.67, -71.67, -74.17, -74.67, -73.67, -73, -75.17, -76.41, -66.67, -66.5, -64.83, -66.17, -68.33, -68.67, -68.33, -69.33, -69.83, -69.33, -71.5, -77.77, -77]
    num_iceshelves = len(lon_min)-1

    # Density of ice in kg/m^3
    rho_ice = 916
    # Degrees to radians conversion factor
    deg2rad = pi/180
    # Northern boundary 63S for plot
    nbdry = -63+90
    # Centre of missing circle in grid
    lon_c = 50
    lat_c = -83
    # Radius of missing circle (play around with this until it works)
    radius = 10.1

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
    # Set up array for mass loss values at each ice shelf
    massloss_ts = empty([num_iceshelves, len(time)])
    index = 0
    # Loop over ice shelves
    while index < num_iceshelves:
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
    # Find the magnitude of the seasonal cycle in mass loss for each ice shelf
    massloss_cycle = empty(num_iceshelves)
    for index in range(num_iceshelves):
        massloss_cycle[index] = amax(massloss_ts[index, min_t:max_t]) - amin(massloss_ts[index, min_t:max_t])

    # Read the grid
    id = Dataset(grid_path, 'r')
    lon = id.variables['lon_rho'][:,:]
    lat = id.variables['lat_rho'][:,:]
    mask_rho = id.variables['mask_rho'][:,:]
    mask_zice = id.variables['mask_zice'][:,:]
    id.close()

    # Calculate dA and mask with zice
    dx, dy = cartesian_grid_2d(lon, lat)
    dA = ma.masked_where(mask_zice==0, dx*dy)
    
    # Make sure longitude goes from -180 to 180, not 0 to 360
    index = lon > 180
    lon[index] = lon[index] - 360

    # Get land/zice mask
    open_ocn = copy(mask_rho)
    open_ocn[mask_zice==1] = 0
    land_zice = ma.masked_where(open_ocn==1, open_ocn)

    # Initialise a field of ice shelf melt rate seasonal cycle
    ismr_cycle = zeros(shape(lon))
    # Mask out land and open ocean points
    ismr_cycle = ma.masked_where(mask_zice==0, ismr_cycle)

    # Loop over ice shelves
    for index in range(num_iceshelves):
        # Find the grid cells within this ice shelf
        if index == num_iceshelves-1:
            # Ross region is split into two
            region = (lon >= lon_min[index])*(lon <= lon_max[index])*(lat >= lat_min[index])*(lat <= lat_max[index])*(mask_zice == 1) + (lon >= lon_min[index+1])*(lon <= lon_max[index+1])*(lat >= lat_min[index+1])*(lat <= lat_max[index+1])*(mask_zice == 1)
        else:
            region = (lon >= lon_min[index])*(lon <= lon_max[index])*(lat >= lat_min[index])*(lat <= lat_max[index])*(mask_zice == 1)
        # Calculate the conversion factor from mass loss to area-averaged
        # melt rate for this ice shelf
        factor = 1e12/(rho_ice*sum(dA[region]))
        # Get seasonal cycle in melt rates and modify the field for this region
        ismr_cycle[region] = massloss_cycle[index]*factor

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

    # Set up colour levels
    lev = linspace(0, amax(ismr_cycle), num=50)
    print amax(ismr_cycle)

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
    contourf(x, y, ismr_cycle, lev, cmap='jet')    
    cbar = colorbar()
    cbar.ax.tick_params(labelsize=20)
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    title('Seasonal cycle in area-averaged ice shelf melt rate (m/y)', fontsize=30)
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
    ismr_seasonal_cycle(grid_path, log_path, save, fig_name)


            
    
    
    
