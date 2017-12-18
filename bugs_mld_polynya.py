from netCDF4 import Dataset, num2date
from numpy import *
from matplotlib.pyplot import *

def bugs_mld_polynya ():

    # Files and timesteps to read
    file_ocn = '/short/m68/kaa561/metroms_iceshelf/tmproms/run/bug_chapter/no_restoring/ocean_avg_0012.nc'
    tstep_ocn = 2
    file_ice = '/short/m68/kaa561/metroms_iceshelf/tmproms/run/bug_chapter/no_restoring/cice/rundir/history/iceh.1994-09-26.nc'
    tstep_ice = 1
    # Degrees to radians conversion factor
    deg2rad = pi/180
    # Month names for titles
    month_names = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
    # Maximum latitude to plot
    nbdry = -52+90    

    # Set up figure
    fig = figure(figsize=(16,8))
    gs = GridSpec(1,2)
    gs.update(left=0.1, right=0.9, bottom=0.05, top=0.9, wspace=0.05)

    # Read mixed layer depth (actually surface boundary layer depth)
    # and ROMS grid
    id = Dataset(file_ocn, 'r')
    lon = id.variables['lon_rho'][:,:-1]
    lat = id.variables['lat_rho'][:,:-1]
    zice = id.variables['zice'][:,:-1]
    hsbl = id.variables['Hsbl'][tstep_ocn-1,:,:-1]
    # Also read time axis and convert to Date objects
    time_id = id.variables['ocean_time']
    time = num2date(time_id[tstep_ocn-1], units=time_id.units, calendar=time_id.calendar.lower())
    id.close()
    # Mask out the ice shelves and change the sign
    mld = ma.masked_where(zice!=0, -hsbl)
    # Polar coordinate transformation
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)
    # Get the date for the title
    date_string = str(time.day) + ' ' + month_names[time.month-1] + ' ' + str(time.year)
    # Plot
    ax = subplot(gs[0,0], aspect='equal')
    lev = linspace(0, amax(mld), num=50)
    contourf(x, y, mld, lev, cmap='jet')
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    ax.set_xticks([])
    ax.set_yticks([])
    title('a) Surface boundary layer depth (m)', fontsize=24)
    # Colourbar
    cbaxes = fig.add_axes([0.03, 0.3, 0.02, 0.4])
    cbar = colorbar(cax=cbaxes, ticks=arange(0,2000+500,500))
    cbar.ax.tick_params(labelsize=16)
    
    # Read sea ice concentration and CICE grid
    id = Dataset(file_ice, 'r')
    lon_tmp = id.variables['TLON'][:,:]
    lat_tmp = id.variables['TLAT'][:,:]
    aice_tmp = id.variables['aice'][tstep_ice-1,:,:]
    id.close()
    # Wrap the periodic boundray by 1 cell
    lon = ma.empty([size(lon_tmp,0), size(lon_tmp,1)+1])
    lat = ma.empty([size(lat_tmp,0), size(lat_tmp,1)+1])
    aice = ma.empty([size(aice_tmp,0), size(aice_tmp,1)+1])
    lon[:,:-1] = lon_tmp
    lon[:,-1] = lon_tmp[:,0]
    lat[:,:-1] = lat_tmp
    lat[:,-1] = lat_tmp[:,0]
    aice[:,:-1] = aice_tmp
    aice[:,-1] = aice_tmp[:,0]
    # Convert to spherical coordinates
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)
    # Plot
    ax = subplot(gs[0,1], aspect='equal')
    lev = linspace(0, 1, num=50)
    contourf(x, y, aice, lev, cmap='jet')
    xlim([-nbdry, nbdry])
    ylim([-nbdry, nbdry])
    ax.set_xticks([])
    ax.set_yticks([])
    title('b) Sea ice concentration', fontsize=24)
    # Colourbar
    cbaxes = fig.add_axes([0.92, 0.3, 0.02, 0.4])
    cbar = colorbar(cax=cbaxes, ticks=arange(0,1+0.25,0.25))
    cbar.ax.tick_params(labelsize=16)
    # Label the timestep
    text(0.5, 0.94, date_string, ha='center', transform=fig.transFigure, fontsize=28)
    fig.show()
    fig.savefig('bugs_mld_polynya.png')


# Command-line interface
if __name__ == "__main__":

    bugs_mld_polynya()
    
    
