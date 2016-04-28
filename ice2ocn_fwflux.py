from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *

# Creates a circumpolar Antarctic plot of net CICE-to-ROMS freshwater flux
# in cm/day. Follows the same process as circumpolar_cice_plot.py, but
# the derived variable (FW flux - salt flux, converted to cm/day)
# requires a special case.
# Input:
# file_path = path to CICE history file
# tstep = timestep to plot in file_path
# fig_name = filename for figure
def ice2ocn_fwflux (file_path, tstep, fig_name):

    deg2rad = pi/180

    # Read data
    id = Dataset(file_path, 'r')
    data = id.variables['fresh_ai'][tstep-1,:,:] - id.variables['fsalt_ai'][tstep-1,:,:]/1000.*60.*60.*24.*100.
    lon = id.variables['TLON'][:,:]
    lat = id.variables['TLAT'][:,:]
    id.close()

    # Convert to spherical coordinates
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)
    lev = linspace(amin(data),amax(data),40)

    # Plot
    fig = figure(figsize=(16,12))
    fig.add_subplot(1,1,1, aspect='equal')
    contourf(x, y, data, lev)
    cbar = colorbar()
    cbar.ax.tick_params(labelsize=20)
    title('CICE-to-ROMS net freshwater flux (cm/day)', fontsize=30)
    axis('off')

    #savefig(fig_name)
    show()


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to CICE history file: ")
    tstep = int(raw_input("Timestep number (starting at 1): "))
    fig_name = raw_input("File name for figure: ")
    ice2ocn_fwflux(file_path, tstep, fig_name)        
        
        

    
    
