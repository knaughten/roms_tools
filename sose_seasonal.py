from netCDF4 import Dataset, num2date
from numpy import *
from matplotlib.pyplot import *
from calc_z import *
from interp_lon_sose import *

# Make a 4x2 plot showing lat vs. depth slices of seasonally averaged 
# temperature (top) and salinity (bottom) at the given longitude, for the
# 2005-2010 SOSE climatology.
# Input:
# lon0 = the specific longitude to plot (between -180 and 180)
# depth_bdry = deepest depth to plot (negative, in m)
# save = optional boolean flag; if True, the figure will be saved with file name
#        fig_name; if False, the figure will display on the screen
# fig_name = optional string containing filename for figure, if save=True
def sose_seasonal (lon0, depth_bdry, save=False, fig_name=None):

    # Path to SOSE seasonal climatology file
    sose_file = '../SOSE_seasonal_climatology.nc'
    # Season names for titles
    season_names = ['DJF', 'MAM', 'JJA', 'SON']

    # Bounds on colour scale
    temp_min = -2
    temp_max = 3
    temp_ticks = 1
    salt_min = 33.8
    salt_max = 34.8
    salt_ticks = 0.4

    # Choose what to write on the title about longitude
    if lon0 < 0:
        lon_string = str(int(round(-lon0))) + r'$^{\circ}$W'
    else:
        lon_string = str(int(round(lon0))) + r'$^{\circ}$E'
    # Edit longitude bounds to be from 0 to 360, to fit with ROMS convention
    if lon0 < 0:
        lon0 += 360

    print 'Processing SOSE data'

    # Read grid and 3D data (already seasonally averaged)
    id = Dataset(sose_file, 'r')
    lon_sose = id.variables['longitude'][0,:]
    lat_sose = id.variables['latitude'][:,0]
    z_sose = id.variables['depth'][:]
    temp_3d_sose = id.variables['temp'][:,:,:,:]
    salt_3d_sose = id.variables['salt'][:,:,:,:]

    # Calculate zonal slices for each season
    temp_sose = ma.empty([4, size(z_sose), size(lat_sose,0)])
    temp_sose[:,:,:] = 0.0
    salt_sose = ma.empty([4, size(z_sose), size(lat_sose,0)])
    salt_sose[:,:,:] = 0.0
    for season in range(4):
        print 'Calculating zonal slices for ' + season_names[season]
        temp_sose[season,:,:] = interp_lon_sose(temp_3d_sose[season,:,:,:], lon_sose, lon0)
        salt_sose[season,:,:] = interp_lon_sose(salt_3d_sose[season,:,:,:], lon_sose, lon0)

    # Set colour levels
    lev1 = linspace(temp_min, temp_max, num=50)
    lev2 = linspace(salt_min, salt_max, num=50)

    # Choose southern boundary based on extent of SOSE grid
    sbdry = -80
    nbdry = -30

    # Plot
    print 'Plotting'
    fig = figure(figsize=(20,9))
    # Loop over seasons
    for season in range(4):
        # Temperature
        fig.add_subplot(2, 4, season+1)
        img = contourf(lat_sose, z_sose, temp_sose[season,:,:], lev1, cmap='jet', extend='both')
        xlim([sbdry, nbdry])
        ylim([depth_bdry, 0])
        title('Temperature (' + season_names[season] + ')', fontsize=24)
        if season == 0:
            ylabel('depth (m)', fontsize=18)
        if season == 3:
            # Add colourbar
            cbaxes1 = fig.add_axes([0.92, 0.55, 0.01, 0.3])
            cbar1 = colorbar(img, cax=cbaxes1, ticks=arange(temp_min, temp_max+temp_ticks, temp_ticks))
            cbar1.ax.tick_params(labelsize=16)
        # Salinity
        fig.add_subplot(2, 4, season+5)
        img = contourf(lat_sose, z_sose, salt_sose[season,:,:], lev2, cmap='jet', extend='both')
        xlim([sbdry, nbdry])
        ylim([depth_bdry, 0])
        title('Salinity (' + season_names[season] + ')', fontsize=24)
        if season == 0:
            ylabel('depth (m)', fontsize=18)
        xlabel('Latitude', fontsize=18)
        if season == 3:
            # Add colourbar
            cbaxes2 = fig.add_axes([0.92, 0.15, 0.01, 0.3])
            cbar2 = colorbar(img, cax=cbaxes2, ticks=arange(salt_min, salt_max+salt_ticks, salt_ticks))
            cbar2.ax.tick_params(labelsize=16)
    # Add the main title
    suptitle(lon_string, fontsize=30)

    # Finished
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()    


# Command-line interface
if __name__ == "__main__":

    lon0 = float(raw_input("Enter longitude (-180 to 180): "))
    depth_bdry = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    sose_seasonal(lon0, depth_bdry, save, fig_name)
