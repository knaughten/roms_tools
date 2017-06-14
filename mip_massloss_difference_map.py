from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib import rcParams

# Make a circumpolar Antarctic figure showing the percentage change in mass
# loss for each ice shelf in the FESOM simulation with respect to the MetROMS
# simulation, for the 2003-2008 average.
# Input:
# roms_grid = path to ROMS grid file
# roms_logfile = path to ROMS logfile from timeseries_massloss.py
# fesom_logfile = path to FESOM logfile from timeseries_massloss.py in the
#                 fesomtools repository
# save = optional boolean indicating to save the figure to a file, rather than
#        display on screen
# fig_name = if save=True, filename for figure
def mip_massloss_difference_map (roms_grid, roms_logfile, fesom_logfile, save=False, fig_name=None):

    # Year simulations start
    year_start = 1992
    # Years to analyse 
    obs_start = 2003
    obs_end = 2008
    # Number of output steps per year in FESOM
    peryear = 365/5

    # Limits on longitude and latitude for each ice shelf
    # These depend on the source geometry, in this case RTopo 1.05
    # Note there is one extra index at the end of each array; this is because
    # the Ross region crosses the line 180W and therefore is split into two
    lon_min = [-62.67, -65.5, -79.17, -85, -104.17, -102.5, -108.33, -114.5, -135.67, -149.17, -155, 144, 115, 94.17, 80.83, 65, 33.83, 19, 12.9, 9.33, -10.05, -28.33, -181, 158.33]
    lon_max = [-59.33, -60, -66.67, -28.33, -88.83, -99.17, -103.33, -111.5, -114.33, -140, -145, 146.62, 123.33, 102.5, 89.17, 75, 37.67, 33.33, 16.17, 12.88, 7.6, -10.33, -146.67, 181]
    lat_min = [-73.03, -69.35, -74.17, -83.5, -73.28, -75.5, -75.5, -75.33, -74.9, -76.42, -78, -67.83, -67.17, -66.67, -67.83, -73.67, -69.83, -71.67, -70.5, -70.75, -71.83, -76.33, -85, -84.5]
    lat_max = [-69.37, -66.13, -69.5, -74.67, -71.67, -74.17, -74.67, -73.67, -73, -75.17, -76.41, -66.67, -66.5, -64.83, -66.17, -68.33, -68.67, -68.33, -69.33, -69.83, -69.33, -71.5, -77.77, -77]
    num_shelves = len(lon_min)-1

    # Degrees to radians conversion factor
    deg2rad = pi/180
    # Northern boundary 63S for plot
    nbdry = -63+90
    # Centre of missing circle in grid
    lon_c = 50
    lat_c = -83
    # Radius of missing circle (play around with this until it works)
    radius = 10.1
    # Minimum zice in ROMS grid
    min_zice = -10

    # Read ROMS logfile
    roms_time = []
    f = open(roms_logfile, 'r')
    # Skip the first line (header for time array)
    f.readline()
    for line in f:
        try:
            roms_time.append(float(line))
        except(ValueError):
            # Reached the header for the next variable
            break
    # Skip the values for the entire continent
    for line in f:
        try:
            tmp = float(line)
        except(ValueError):
            break
    # Set up array for mass loss values at each ice shelf
    roms_massloss_ts = empty([num_shelves, len(roms_time)])
    index = 0
    # Loop over ice shelves
    while index < num_shelves:
        t = 0
        for line in f:
            try:
                roms_massloss_ts[index, t] = float(line)
                t += 1
            except(ValueError):
                # Reached the header for the next ice shelf
                break
        index +=1
    # Add start year to ROMS time array
    roms_time = array(roms_time) + year_start
    # Average between observation years
    t_start = nonzero(roms_time >= obs_start)[0][0]
    t_end = nonzero(roms_time >= obs_end+1)[0][0]
    roms_massloss = mean(roms_massloss_ts[:,t_start:t_end], axis=1)

    # Read FESOM timeseries
    tmp = []
    f = open(fesom_logfile, 'r')
    # Skip the first line (header)
    f.readline()
    # Count the number of time indices for the first variable (total mass loss
    # for all ice shelves, which we don't care about)
    num_time = 0
    for line in f:
        try:
            tmp = float(line)
            num_time += 1
        except(ValueError):
            # Reached the header for the next variable
            break
    # Set up array for mass loss values at each ice shelf
    fesom_massloss_ts = empty([num_shelves, num_time])
    # Loop over ice shelves
    index = 0
    while index < num_shelves:
        t = 0
        for line in f:
            try:
                fesom_massloss_ts[index,t] = float(line)
                t += 1
            except(ValueError):
                # Reached the header for the next ice shelf
                break
        index += 1
    # Average between observation years
    fesom_massloss = mean(fesom_massloss_ts[:,peryear*(obs_start-year_start):peryear*(obs_end+1-year_start)], axis=1)

    # Read ROMS grid
    id = Dataset(roms_grid, 'r')
    lon = id.variables['lon_rho'][:-15,:-1]
    lat = id.variables['lat_rho'][:-15,:-1]
    mask_rho = id.variables['mask_rho'][:-15,:-1]
    mask_zice = id.variables['mask_zice'][:-15,:-1]
    zice = id.variables['zice'][:-15,:-1]
    id.close()
    # Make sure longitude goes from -180 to 180, not 0 to 360
    index = lon > 180
    lon[index] = lon[index] - 360
    # Get land/zice mask
    open_ocn = copy(mask_rho)
    open_ocn[mask_zice==1] = 0
    land_zice = ma.masked_where(open_ocn==1, open_ocn)

    # Initialise plotting field of ice shelf mass loss percentage change
    massloss_change = ma.empty(shape(lon))
    massloss_change[:,:] = ma.masked
    # Loop over ice shelves
    for index in range(num_shelves):
        # Calculate percentage change in massloss, FESOM wrt MetROMS
        tmp = (fesom_massloss[index] - roms_massloss[index])/roms_massloss[index]*100
        # Modify plotting field for this region
        if index == num_shelves-1:
            # Ross region is split into two
            region = (lon >= lon_min[index])*(lon <= lon_max[index])*(lat >= lat_min[index])*(lat <= lat_max[index])*(mask_zice == 1) + (lon >= lon_min[index+1])*(lon <= lon_max[index+1])*(lat >= lat_min[index+1])*(lat <= lat_max[index+1])*(mask_zice == 1)
        else:
            region = (lon >= lon_min[index])*(lon <= lon_max[index])*(lat >= lat_min[index])*(lat <= lat_max[index])*(mask_zice == 1)
        massloss_change[region] = tmp

    # Edit zice so tiny ice shelves won't be contoured
    zice[massloss_change.mask] = 0.0
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
    # Set up colour scale
    lev = linspace(-200, 200, num=50)

    # Plot
    fig = figure(figsize=(16,12))
    ax = fig.add_subplot(1,1,1, aspect='equal')
    fig.patch.set_facecolor('white')
    # First shade land and zice in grey (include zice so there are no white
    # patches near the grounding line where contours meet)
    contourf(x, y, land_zice, 1, colors=(('0.6', '0.6', '0.6')))
    # Fill in the missing cicle
    contourf(x_reg, y_reg, land_circle, 1, colors=(('0.6', '0.6', '0.6')))
    # Now shade the percentage change in mass loss
    contourf(x, y, massloss_change, lev, cmap='RdBu_r', extend='both')
    cbar = colorbar(ticks=arange(-200, 200+50, 50))
    cbar.ax.tick_params(labelsize=20)
    # Add a black contour for the ice shelf front
    rcParams['contour.negative_linestyle'] = 'solid'
    contour(x, y, zice, levels=[min_zice], colors=('black'))
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    title('% Change in Ice Shelf Mass Loss (2003-2008 average)\nFESOM with respect to MetROMS', fontsize=30)
    axis('off')

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


if __name__ == "__main__":

    roms_grid = raw_input("Path to ROMS grid file: ")
    roms_logfile = raw_input("Path to ROMS mass loss logfile: ")
    fesom_logfile = raw_input("Path to FESOM mass loss logfile: ")
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    mip_massloss_difference (roms_grid, roms_logfile, fesom_logfile, save, fig_name)
    
