from numpy import *
from netCDF4 import Dataset
from matplotlib.pyplot import *
from seasonal_avg_cice import *

def aice_hs_seasonal (cice_file, save=False, fig_name=None):

    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    deg2rad = pi/180.0

    id = Dataset(cice_file, 'r')
    lon_tmp = id.variables['TLON'][:-15,:]
    lat_tmp = id.variables['TLAT'][:-15,:]
    num_lon = id.variables['TLON'].shape[1]
    num_lat = id.variables['TLAT'].shape[0]
    id.close()
    lon = ma.empty([size(lon_tmp,0), size(lon_tmp,1)+1])
    lat = ma.empty([size(lat_tmp,0), size(lat_tmp,1)+1])
    lon[:,:-1] = lon_tmp
    lon[:,-1] = lon_tmp[:,0]
    lat[:,:-1] = lat_tmp
    lat[:,-1] = lat_tmp[:,0]

    aice_tmp = seasonal_avg_cice(cice_file, 'aice', [num_lat, num_lon])
    hs_tmp = seasonal_avg_cice(cice_file, 'hs', [num_lat, num_lon])
    aice_tmp = aice_tmp[:,:-15,:]
    hs_tmp = hs_tmp[:,:-15,:]
    aice = ma.empty([size(aice_tmp,0), size(aice_tmp,1), size(aice_tmp,2)+1])
    aice[:,:,:-1] = aice_tmp
    aice[:,:,-1] = aice_tmp[:,:,0]
    hs = ma.empty([size(hs_tmp,0), size(hs_tmp,1), size(hs_tmp,2)+1])
    hs[:,:,:-1] = hs_tmp
    hs[:,:,-1] = hs_tmp[:,:,0]

    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)

    bdry1 = -35
    bdry2 = 35
    bdry3 = -33
    bdry4 = 37

    # Set consistent colour levels
    lev1 = linspace(0, 1, num=50)
    lev2 = linspace(0, 0.5, num=50)

    # Plot
    fig = figure(figsize=(20,9))
    # Loop over seasons
    for season in range(4):
        ax = fig.add_subplot(2, 4, season+1, aspect='equal')
        img = contourf(x, y, aice[season,:,:], lev1)
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
        if season == 0:
            text(-39, 0, 'aice (1)', fontsize=21, ha='right')
        title(season_names[season], fontsize=24)
        if season == 3:
            cbaxes1 = fig.add_axes([0.92, 0.55, 0.01, 0.3])
            cbar1 = colorbar(img, ticks=arange(0,1+0.25,0.25), cax=cbaxes1)
            cbar1.ax.tick_params(labelsize=16)
        ax = fig.add_subplot(2, 4, season+5, aspect='equal')
        img = contourf(x, y, hs[season,:,:], lev2, extend='both')
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
        if season == 0:
            text(-39, 0, 'hs (m)', fontsize=21, ha='right')
        if season == 3:
            cbaxes2 = fig.add_axes([0.92, 0.15, 0.01, 0.3])
            cbar2 = colorbar(img, ticks=arange(0,0.5+0.1,0.1), cax=cbaxes2)
            cbar2.ax.tick_params(labelsize=16)
    subplots_adjust(wspace=0.025,hspace=0.025)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    cice_file = raw_input("Path to CICE file, containing at least one complete Dec-Nov period: ")
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    aice_hs_seasonal(cice_file, save, fig_name)

        
        
        

        

        

        

        
        

    
        

    
        

    
    
