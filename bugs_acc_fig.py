from netCDF4 import Dataset, num2date
from numpy import *
from matplotlib.pyplot import *
from rotate_vector_roms import *

def bugs_acc_fig (grid_path, laplacian_file, biharmonic_file, tstep):

    # Month names for title
    month_names = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
    # Degrees to radians conversion factor
    deg2rad = pi/180.0

    # Read angle and grid
    id = Dataset(grid_path, 'r')
    angle = id.variables['angle'][:-15,:]
    lon = id.variables['lon_rho'][:-15,:-1]
    lat = id.variables['lat_rho'][:-15,:-1]
    id.close()

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

    # Make the plot
    fig = figure(figsize=(20,10))
    # Laplacian
    ax = fig.add_subplot(1,2,1, aspect='equal')
    pcolor(x, y, speed_laplacian, vmin=0, vmax=1.25, cmap='cool')
    title('Laplacian viscosity', fontsize=24)
    xlim([-40, 0])
    ylim([0, 40])
    ax.set_xticks([])
    ax.set_yticks([])
    # Biharmonic
    ax = fig.add_subplot(1,2,2, aspect='equal')
    img = pcolor(x, y, speed_biharmonic, vmin=0, vmax=1.25, cmap='cool')
    title('Biharmonic viscosity', fontsize=24)
    xlim([-40, 0])
    ylim([0, 40])
    ax.set_xticks([])
    ax.set_yticks([])
    # Colourbar on the right
    cbaxes = fig.add_axes([0.94, 0.3, 0.02, 0.4])
    cbar = colorbar(img, cax=cbaxes, extend='max', ticks=arange(0,1.25+0.25,0.25))
    cbar.ax.tick_params(labelsize=16)
    suptitle('Surface speed (m/s): snapshot on ' + date_string, fontsize=30)
    subplots_adjust(wspace=0.05)
    fig.show()
    fig.savefig('bugs_acc.png')


# Command-line interface
if __name__ == "__main__":

    grid_path = raw_input("Path to ROMS grid file: ")
    laplacian_file = raw_input("Path to ocean history file for Laplacian simulation: ")
    biharmonic_file = raw_input("Path to same ocean history file for baseline simulation: ")
    tstep = int(raw_input("Timestep to plot (starting at 1): "))
    bugs_acc_fig(grid_path, laplacian_file, biharmonic_file, tstep)
    
