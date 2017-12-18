from netCDF4 import Dataset, num2date
from numpy import *
from matplotlib.pyplot import *
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from rotate_vector_roms import *

# Used ocean_his_0102.nc, timestep 1
def bugs_acc_fig (grid_path, laplacian_file, biharmonic_file, tstep):

    # Month names for title
    month_names = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Bounds on plot (in polar coordinate transformation)
    x_min = -40
    x_max = 0
    y_min = 0
    y_max = 40
    # Minimum speed to plot (mask out values below)
    threshold = 0.12
    # Maximum speed to plot
    bound = 1

    # Read angle and grid
    id = Dataset(grid_path, 'r')
    angle = id.variables['angle'][:-15,:]
    lon = id.variables['lon_rho'][:-15,:-1]
    lat = id.variables['lat_rho'][:-15,:-1]
    mask = id.variables['mask_rho'][:-15,:-1]
    id.close()

    # Set up grey shading of land
    mask = ma.masked_where(mask==1, mask)
    grey_cmap = ListedColormap([(0.6, 0.6, 0.6)])
    x_reg, y_reg = meshgrid(linspace(x_min, x_max, num=1000), linspace(y_min, y_max, num=1000))
    land_circle = zeros(shape(x_reg))
    lon_c = 50
    lat_c = -83
    radius = 10.1
    x_c = -(lat_c+90)*cos(lon_c*deg2rad+pi/2)
    y_c = (lat_c+90)*sin(lon_c*deg2rad+pi/2)
    land_circle = ma.masked_where(sqrt((x_reg-x_c)**2 + (y_reg-y_c)**2) > radius, land_circle)
    # Truncate colourmap
    min_colour = threshold/bound
    max_colour = 1
    trunc_cmap = truncate_colormap(get_cmap('jet'), min_colour, max_colour)

    # Read surface velocity
    id = Dataset(laplacian_file, 'r')
    ur_lap = id.variables['u'][tstep-1,-1,:-15,:]
    vr_lap = id.variables['v'][tstep-1,-1,:-15,:]
    # Rotate to lon-lat space
    u_lap, v_lap = rotate_vector_roms(ur_lap, vr_lap, angle)
    # Now that they're on the same grid, get the speed
    speed_laplacian = sqrt(u_lap**2 + v_lap**2)    
    # Also read time and convert to Date object
    time_id = id.variables['ocean_time']
    time = num2date(time_id[tstep-1], units=time_id.units, calendar=time_id.calendar.lower())
    # Get the date for the title
    date_string = str(time.day) + ' ' + month_names[time.month-1] + ' ' + str(time.year)
    id.close()
    # Repeat velocity for biharmonic simulation
    id = Dataset(biharmonic_file, 'r')
    ur_bih = id.variables['u'][tstep-1,-1,:-15,:]
    vr_bih = id.variables['v'][tstep-1,-1,:-15,:]
    u_bih, v_bih = rotate_vector_roms(ur_bih, vr_bih, angle)
    speed_biharmonic = sqrt(u_bih**2 + v_bih**2)
    id.close()

    # Calculate x and y coordinates for plotting circumpolar projection
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y  = (lat+90)*sin(lon*deg2rad+pi/2)
    speed_laplacian = ma.masked_where(speed_laplacian < threshold, speed_laplacian)
    speed_biharmonic = ma.masked_where(speed_biharmonic < threshold, speed_biharmonic)

    # Make the plot
    fig = figure(figsize=(20,10))
    # Laplacian
    ax = fig.add_subplot(1,2,1, aspect='equal')
    # First shade land in grey
    pcolor(x, y, mask, cmap=grey_cmap)
    pcolor(x_reg, y_reg, land_circle, cmap=grey_cmap)
    pcolor(x, y, speed_laplacian, vmin=threshold, vmax=bound, cmap=trunc_cmap)
    title('Laplacian viscosity', fontsize=24)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    ax.set_xticks([])
    ax.set_yticks([])
    # Biharmonic
    ax = fig.add_subplot(1,2,2, aspect='equal')
    pcolor(x, y, mask, cmap=grey_cmap)
    pcolor(x_reg, y_reg, land_circle, cmap=grey_cmap)
    img = pcolor(x, y, speed_biharmonic, vmin=threshold, vmax=bound, cmap=trunc_cmap)
    title('Biharmonic viscosity', fontsize=24)
    xlim([-40, 0])
    ylim([0, 40])
    ax.set_xticks([])
    ax.set_yticks([])
    # Colourbar on the right
    cbaxes = fig.add_axes([0.94, 0.3, 0.02, 0.4])
    cbar = colorbar(img, cax=cbaxes, extend='max', ticks=arange(0.2,1+0.2,0.2))
    cbar.ax.tick_params(labelsize=16)
    suptitle('Surface speed (m/s): snapshot on ' + date_string, fontsize=30)
    subplots_adjust(wspace=0.05)
    fig.show()
    fig.savefig('bugs_acc.png')


# Truncate colourmap function from https://stackoverflow.com/questions/40929467/how-to-use-and-plot-only-a-part-of-a-colorbar-in-matplotlib
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=-1):
    if n== -1:
        n = cmap.N
    new_cmap = LinearSegmentedColormap.from_list('trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval), cmap(linspace(minval, maxval, n)))
    return new_cmap


# Command-line interface
if __name__ == "__main__":

    grid_path = raw_input("Path to ROMS grid file: ")
    laplacian_file = raw_input("Path to ocean history file for Laplacian simulation: ")
    biharmonic_file = raw_input("Path to same ocean history file for baseline simulation: ")
    tstep = int(raw_input("Timestep to plot (starting at 1): "))
    bugs_acc_fig(grid_path, laplacian_file, biharmonic_file, tstep)
    
