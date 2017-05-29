from netCDF4 import *
from numpy import *
from matplotlib.pyplot import *

# Create a 1x2 plot of sea ice concentration near the Amery Ice Shelf
# on 23 August (the sea ice area max), showing the difference in polynya size
# between the U3_LIM and C4_LD simulations.
def adv_polynyas ():

    # Paths to simulation directories
    paths = ['/short/m68/kaa561/advection/u3_lim/', '/short/m68/kaa561/advection/c4_l/']
    # Titles for plotting
    labels = ['a) U3_LIM', 'b) C4_LD']
    # File name: daily average for 23 August
    file_tail = 'cice/rundir/history/iceh.1992-08-23.nc'
    # Longitude and latitude bounds
    lon_min = 67
    lon_max = 86
    lat_min = -70
    lat_max = -65

    # Set up figure
    fig = figure(figsize=(18,6))
    # Loop over simulations
    for sim in range(2):
        # Read sea ice concentration in the region of interest
        id = Dataset(paths[sim]+file_tail, 'r')
        data = id.variables['aice'][0,50:250,250:400]
        if sim == 0:
            # For the first simulation, also read the grid
            lon = id.variables['TLON'][50:250,250:400]
            lat = id.variables['TLAT'][50:250,250:400]
        id.close()
        ax = fig.add_subplot(1, 2, sim+1)
        # Shade the data (pcolor not contourf so we can show each individual
        # model cell)
        img = pcolor(lon, lat, data, vmin=0, vmax=1, cmap='jet')
        # Configure plot
        title(labels[sim], fontsize=20)
        xlabel('Longitude', fontsize=18)
        ylabel('Latitude', fontsize=18)
        xlim([lon_min, lon_max])
        ylim([lat_min, lat_max])
        if sim == 1:
            # Add a colorbar on the right
            cbaxes = fig.add_axes([0.93, 0.25, 0.015, 0.5])
            cbar = colorbar(ticks=arange(0, 1+0.25, 0.25), cax=cbaxes)
            cbar.ax.tick_params(labelsize=16)

        # Set ticks the way we want them            
        lon_ticks = arange(lon_min+3, lon_max, 5)
        ax.set_xticks(lon_ticks)
        lon_labels = []
        for val in lon_ticks:
            lon_labels.append(str(int(round(val))) + r'$^{\circ}$E')
        ax.set_xticklabels(lon_labels, fontsize=16)
        lat_ticks = arange(lat_min, lat_max, 2)
        ax.set_yticks(lat_ticks)
        lat_labels = []
        for val in lat_ticks:
            lat_labels.append(str(int(round(-val))) + r'$^{\circ}$S')
        ax.set_yticklabels(lat_labels, fontsize=16)
    # Add main title
    suptitle('Sea ice concentration on 23 August', fontsize=22)

    #fig.show()
    fig.savefig('adv_polynyas.png')


# Command-line interface
if __name__ == "__main__":

    adv_polynyas()
            
    
