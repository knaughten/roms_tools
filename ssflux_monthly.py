from netCDF4 import Dataset, num2date
from numpy import *
from matplotlib.pyplot import *
from monthly_avg_roms import *

def ssflux_monthly (roms_file, month, bound=None, save=False, fig_name=None):

    # Month names for plot titles
    month_name = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
    deg2rad = pi/180

    # Read ROMS grid
    id = Dataset(roms_file, 'r')
    lon = id.variables['lon_rho'][:-15,:-1]
    lat = id.variables['lat_rho'][:-15,:-1]
    num_lon = id.variables['lon_rho'].shape[1]
    num_lat = id.variables['lat_rho'].shape[0]
    id.close()

    ssflux = monthly_avg_roms(roms_file, 'ssflux', [num_lat, num_lon], month)
    ssflux = ssflux[:-15,:-1]
    ssflux *= 1e6

    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)

    if bound is None:
        bound = amax(abs(ssflux))
    lev = linspace(-bound, bound, num=50)
    
    fig = figure(figsize=(16,12))
    fig.add_subplot(1,1,1,aspect='equal')
    img = contourf(x, y, ssflux, lev, cmap='RdYlBu_r', extend='both')
    cbar = colorbar()
    cbar.ax.tick_params(labelsize=20)
    title(month_name[month] + r' surface salinity flux (10$^{-6}$ kg/m$^2$/s)', fontsize=30)
    axis('off')

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


if __name__ == "__main__":

    roms_file = raw_input("Path to ROMS output file: ")
    month = int(raw_input("Month number (1-12): "))-1
    bound = None
    get_bound = raw_input("Set bounds on colour scale (y/n)? ")
    if get_bound == 'y':
        bound = float(raw_input("Maximum absolute value (1e-6 kg/m^2/s): "))
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    ssflux_monthly(roms_file, month, bound, save, fig_name)

    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                changes = raw_input("Enter a parameter to change: (1) file path, (2) month, (3) colour bounds, (4) save/display; or enter to continue: ")
                if len(changes) == 0:
                    break
                else:
                    if int(changes) == 1:
                        roms_file = raw_input("Path to ROMS output file: ")
                    elif int(changes) == 2:
                        month = int(raw_input("Month number (1-12): "))-1
                    elif int(changes) == 3:
                        get_bound = raw_input("Set bounds on colour scale (y/n)? ")
                        if get_bound == 'y':
                            bound = float(raw_input("Maximum absolute value (1e-6 kg/m^2/s): "))
                    elif int(changes) == 4:
                        save = not save
            if save:
                fig_name = raw_input("File name for figure: ")
            ssflux_monthly(roms_file, month, bound, save, fig_name)
        else:
            break
                

    
    



