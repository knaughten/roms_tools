from netCDF4 import *
from numpy import *
from matplotlib.pyplot import *

def adv_polynyas ():

    paths = ['/short/m68/kaa561/advection/u3_lim/', '/short/m68/kaa561/advection/c4_l/']
    labels = ['a) U3_LIM', 'b) C4_LD']
    file_tail = 'iceh.1992-08-23.nc'
    lon_min = 100
    lon_max = 140
    lat_min = -67.1
    lat_max = -64.9
    lat_min_label = -67
    lat_max_label = -65

    fig = figure(figsize=(18,6))
    for sim in range(2):
        id = Dataset(paths[sim]+file_tail, 'r')
        data = id.variables['aice'][0,70:280,350:750]
        if sim == 0:
            lon = id.variables['TLON'][70:280,350:750]
            lat = id.variables['TLAT'][70:280,350:750]
        id.close()
        ax = fig.add_subplot(1, 2, sim+1)
        img = pcolor(lon, lat, data, vmin=0, vmax=1, cmap='jet')
        title(labels[sim], fontsize=20)
        xlabel('Longitude', fontsize=16)
        ylabel('Latitude', fontsize=16)
        xlim([lon_min, lon_max])
        ylim([lat_min, lat_max])
        if sim == 1:
            cbaxes = fig.add_axes([0.93, 0.25, 0.015, 0.5])
            cbar = colorbar(ticks=arange(0, 1+0.25, 0.25), cax=cbaxes)
            cbar.ax.tick_params(labelsize=14)

        lon_ticks = arange(lon_min, lon_max+10, 10)
        ax.set_xticks(lon_ticks)
        lon_labels = []
        for val in lon_ticks:
            lon_labels.append(str(int(round(val))) + r'$^{\circ}$E')
        ax.set_xticklabels(lon_labels, fontsize=14)
        lat_ticks = arange(lat_min_label, lat_max_label+1, 1)
        ax.set_yticks(lat_ticks)
        lat_labels = []
        for val in lat_ticks:
            lat_labels.append(str(int(round(-val))) + r'$^{\circ}$S')
        ax.set_yticklabels(lat_labels, fontsize=14)
    suptitle('Sea ice concentration on 23 August', fontsize=22)

    #fig.show()
    fig.savefig('adv_polynyas.png')


if __name__ == "__main__":

    adv_polynyas()
            
    
