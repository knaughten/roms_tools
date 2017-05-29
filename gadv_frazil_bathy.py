from netCDF4 import *
from numpy import *
from matplotlib.pyplot import *

# Plot a section of the coastline from 25-45E. On the left, show anomalies
# in annually averaged frazil formation between the U3 and U3_LIM simulations.
# On the right, show bathymetry. This way the relationship between spurious
# supercooling and steep bathymetry is obvious.
def gadv_frazil_bathy ():

    # Paths to simulation directories
    path_up3l = '/short/m68/kaa561/advection/u3_lim/'
    path_up3 = '/short/m68/kaa561/advection/u3/'
    # Path to ROMS grid (containing bathymetry data)
    roms_grid = '/short/m68/kaa561/metroms_iceshelf/apps/common/grid/circ30S_quarterdegree_old.nc'
    # Filename for annually averaged CICE data
    file_tail = 'cice/rundir/history/iceh_avg.nc'
    # Bounds on lat and lon
    lon_min = 25
    lon_max = 45
    lat_min = -71
    lat_max = -64
    # Bounds for colour scales
    max_frazil = 0.5
    tick_frazil = 0.25
    max_bathy = 5
    tick_bathy = 1

    # Read U3_LIM grid and frazil data
    id = Dataset(path_up3l + file_tail, 'r')
    lon_cice = id.variables['TLON'][50:250,0:200]
    lat_cice = id.variables['TLAT'][50:250,0:200]
    frazil0 = id.variables['frazil'][0,50:250,0:200]
    id.close()
    # Read U3 frazil data
    id = Dataset(path_up3 + file_tail, 'r')
    frazil1 = id.variables['frazil'][0,50:250,0:200]
    id.close()
    # Calculate anomaly
    frazil_anom = frazil1 - frazil0

    # Read bathymetry
    id = Dataset(roms_grid, 'r')
    lon_roms = id.variables['lon_rho'][51:251,1:201]
    lat_roms = id.variables['lat_rho'][51:251,1:201]
    bathy = id.variables['h'][51:251,:201]*1e-3
    # Mask out land and ice shelf points
    mask = id.variables['mask_rho'][51:251,:201]-id.variables['mask_zice'][51:251,:201]
    id.close()
    bathy = ma.masked_where(mask==0, bathy)

    # Make the plot
    fig = figure(figsize=(18,8))
    # Frazil anomalies
    ax = fig.add_subplot(1, 2, 1)
    img1 = pcolor(lon_cice, lat_cice, frazil_anom, vmin=-max_frazil, vmax=max_frazil, cmap='RdBu_r')
    title('a) Frazil ice formation (cm/day), annual average\nAnomalies: UP3 minus UP3L', fontsize=20)
    xlabel('Longitude', fontsize=18)    
    xlim([lon_min, lon_max])
    ylim([lat_min, lat_max])
    # Colour bar
    cbaxes1 = fig.add_axes([0.05, 0.25, 0.015, 0.5])
    cbar1 = colorbar(img1, ticks=arange(-max_frazil, max_frazil+tick_frazil, tick_frazil), cax=cbaxes1, extend='both')
    cbar1.ax.tick_params(labelsize=16)
    # Make sure lon and lat ticks are the way we want them
    lon_ticks = arange(lon_min, lon_max+5, 5)
    ax.set_xticks(lon_ticks)
    lon_labels = []
    for val in lon_ticks:
        lon_labels.append(str(int(round(val))) + r'$^{\circ}$E')
    ax.set_xticklabels(lon_labels, fontsize=16)
    lat_ticks = arange(lat_min+1, lat_max+1, 2)
    ax.set_yticks(lat_ticks)
    lat_labels = []
    for val in lat_ticks:
        # No latitude labels on this plot, they bump into the colourbar
        lat_labels.append('')
    ax.set_yticklabels(lat_labels, fontsize=16)
    # Bathymetry
    ax = fig.add_subplot(1, 2, 2)
    img2 = pcolor(lon_roms, lat_roms, bathy, vmin=0, vmax=max_bathy, cmap='jet')
    title('b) Bathymetry (km)', fontsize=20)
    xlabel('Longitude', fontsize=18)
    ylabel('Latitude', fontsize=18)
    xlim([lon_min, lon_max])
    ylim([lat_min, lat_max])
    cbaxes2 = fig.add_axes([0.93, 0.25, 0.015, 0.5])
    cbar2 = colorbar(img2, ticks=arange(0, max_bathy+tick_bathy, tick_bathy), cax=cbaxes2, extend='max')
    cbar2.ax.tick_params(labelsize=16)
    lon_ticks = arange(lon_min, lon_max+5, 5)
    ax.set_xticks(lon_ticks)
    lon_labels = []
    for val in lon_ticks:
        lon_labels.append(str(int(round(val))) + r'$^{\circ}$E')
    ax.set_xticklabels(lon_labels, fontsize=16)
    lat_ticks = arange(lat_min+1, lat_max+1, 2)
    ax.set_yticks(lat_ticks)
    lat_labels = []
    for val in lat_ticks:
        lat_labels.append(str(int(round(-val))) + r'$^{\circ}$S')
    ax.set_yticklabels(lat_labels, fontsize=16)

    #fig.show()
    fig.savefig('kn_fig3.png')


# Command-line interface
if __name__ == "__main__":

    gadv_frazil_bathy()
    
