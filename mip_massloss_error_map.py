from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib import rcParams

# Make a 2x1 circumpolar Antarctic map of percentage error in mass loss for
# each ice shelf, outside the bounds given by Rignot et al. (2013), in MetROMS
# (left) and FESOM (right), for the 2003-2008 average.
# Input:
# roms_grid = path to ROMS grid file
# roms_logfile = path to ROMS logfile from timeseries_massloss.py
# fesom_logfile = path to FESOM logfile from timeseries_massloss.py in the
#                 fesomtools repository
# save = optional boolean indicating to save the figure to a file, rather than
#        display on screen
# fig_name = if save=True, filename for figure
def mip_massloss_error_map (roms_grid, roms_logfile, fesom_logfile, save=False, fig_name=None):

    # Year simulations start
    year_start = 1992
    # Years observations span
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
    # Observed mass loss (Rignot 2013) and uncertainty for each ice shelf, in Gt/y
    obs_massloss = [1.4, 20.7, 135.4, 155.4, 51.8, 101.2, 97.5, 45.2, 144.9, 4.2, 18.2, 7.9, 90.6, 72.6, 27.2, 35.5, -2, 21.6, 6.3, 3.9, 26.8, 9.7, 47.7]
    obs_massloss_error = [14, 67, 40, 45, 19, 8, 7, 4, 14, 2, 3, 3, 8, 15, 10, 23, 3, 18, 2, 2, 14, 16, 34]
    num_shelves = len(obs_massloss)

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

    # Initialise fields of ice shelf mass loss unexplained percent error in each mode
    roms_error = ma.empty(shape(lon))
    fesom_error = ma.empty(shape(lon))
    roms_error[:,:] = ma.masked
    fesom_error[:,:] = ma.masked
    # Loop over ice shelves
    for index in range(num_shelves):
        # Find the range of observations
        massloss_low = obs_massloss[index] - obs_massloss_error[index]
        massloss_high = obs_massloss[index] + obs_massloss_error[index]
        # Find the unexplained percent error in mass loss
        # ROMS
        if roms_massloss[index] < massloss_low:
            # Simulated mass loss too low
            error_tmp = (roms_massloss[index] - massloss_low)/massloss_low*100
        elif roms_massloss[index] > massloss_high:
            # Simulated mass loss too high
            error_tmp = (roms_massloss[index] - massloss_high)/massloss_high*100
        else:
            # Simulated mass loss within observational error estimates
            error_tmp = 0
        # Modify error field for this region
        if index == num_shelves-1:
            # Ross region is split into two
            region = (lon >= lon_min[index])*(lon <= lon_max[index])*(lat >= lat_min[index])*(lat <= lat_max[index])*(mask_zice == 1) + (lon >= lon_min[index+1])*(lon <= lon_max[index+1])*(lat >= lat_min[index+1])*(lat <= lat_max[index+1])*(mask_zice == 1)
        else:
            region = (lon >= lon_min[index])*(lon <= lon_max[index])*(lat >= lat_min[index])*(lat <= lat_max[index])*(mask_zice == 1)
        roms_error[region] = error_tmp
        # Repeat for FESOM
        if fesom_massloss[index] < massloss_low:
            # Simulated mass loss too low
            error_tmp = (fesom_massloss[index] - massloss_low)/massloss_low*100
        elif fesom_massloss[index] > massloss_high:
            # Simulated mass loss too high
            error_tmp = (fesom_massloss[index] - massloss_high)/massloss_high*100
        else:
            # Simulated mass loss within observational error estimates
            error_tmp = 0
        # Modify error field for this region
        if index == num_shelves-1:
            # Ross region is split into two
            region = (lon >= lon_min[index])*(lon <= lon_max[index])*(lat >= lat_min[index])*(lat <= lat_max[index])*(mask_zice == 1) + (lon >= lon_min[index+1])*(lon <= lon_max[index+1])*(lat >= lat_min[index+1])*(lat <= lat_max[index+1])*(mask_zice == 1)
        else:
            region = (lon >= lon_min[index])*(lon <= lon_max[index])*(lat >= lat_min[index])*(lat <= lat_max[index])*(mask_zice == 1)
        fesom_error[region] = error_tmp

    # Edit zice so tiny ice shelves won't be contoured
    zice[roms_error.mask] = 0.0
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
    lev = linspace(-100, 100, num=50)

    # Plot
    fig = figure(figsize=(20,9))
    fig.patch.set_facecolor('white')
    # ROMS
    ax1 = fig.add_subplot(1,2,1, aspect='equal')
    # First shade land and zice in grey (include zice so there are no white
    # patches near the grounding line where contours meet)
    contourf(x, y, land_zice, 1, colors=(('0.6', '0.6', '0.6')))
    # Fill in the missing cicle
    contourf(x_reg, y_reg, land_circle, 1, colors=(('0.6', '0.6', '0.6')))
    # Now shade the error in mass loss
    contourf(x, y, roms_error, lev, cmap='RdBu_r', extend='both')    
    # Add a black contour for the ice shelf front
    rcParams['contour.negative_linestyle'] = 'solid'
    contour(x, y, zice, levels=[min_zice], colors=('black'))
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    title('MetROMS', fontsize=30)
    axis('off')
    # FESOM
    ax2 = fig.add_subplot(1,2,2, aspect='equal')
    contourf(x, y, land_zice, 1, colors=(('0.6', '0.6', '0.6')))
    contourf(x_reg, y_reg, land_circle, 1, colors=(('0.6', '0.6', '0.6')))
    img = contourf(x, y, fesom_error, lev, cmap='RdBu_r', extend='max')   
    rcParams['contour.negative_linestyle'] = 'solid'
    contour(x, y, zice, levels=[min_zice], colors=('black'))
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    axis('off')
    title('FESOM', fontsize=30)
    # Add a colourbar on the right
    cbaxes = fig.add_axes([0.92, 0.2, 0.01, 0.6])
    cbar = colorbar(img, cax=cbaxes, ticks=arange(-100, 100+25, 25))
    cbar.ax.tick_params(labelsize=20)
    # Main title
    suptitle('Bias in Ice Shelf Mass Loss (%)\n2003-2008', fontsize=24)
    subplots_adjust(wspace=0.05)

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
    mip_massloss_error_map (roms_grid, roms_logfile, fesom_logfile, save, fig_name)
    
