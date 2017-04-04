from netCDF4 import Dataset, num2date
from numpy import *
from matplotlib.pyplot import *
from cartesian_grid_2d import *
from rotate_vector_roms import *

# Recreate Figure 1 of Holland et al 2014 (doi:10.1175/JCLI-D-13-00301.1) using
# ROMS output: annually averaged barotropic streamfunction, JJA mean mixed-layer
# depth (just plot the KPP boundary layer depth for now), annually averaged
# bottom water temperature and salinity.
# Input:
# grid_path = path to ROMS grid file
# file_path = path to ROMS output file containing 1 year of 5-day averages,
#             including a continuous June-August period
def holland_fig1 (grid_path, file_path):

    deg2rad = pi/180.0

    # Read grid
    id = Dataset(grid_path, 'r')
    lon = id.variables['lon_rho'][:-15,1:]
    lat = id.variables['lat_rho'][:-15,1:]
    h = id.variables['h'][:-15,1:]
    zice = id.variables['zice'][:-15,1:]
    angle = id.variables['angle'][:-15,:]
    id.close()

    # Set up figure
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)
    fig = figure(figsize=(16,12))

    # Barotropic streamfunction
    # First read bartropic velocity vector
    id = Dataset(file_path, 'r')
    ubar_xy = mean(id.variables['ubar'][:,:-15,:], axis=0)
    vbar_xy = mean(id.variables['vbar'][:,:-15,:], axis=0)
    id.close()
    # Rotate to lon-lat space
    ubar, vbar = rotate_vector_roms(ubar_xy, vbar_xy, angle)
    # Throw away the overlapping periodic boundary
    ubar = ubar[:,1:]
    # Mask ice shelves
    ubar = ma.masked_where(zice!=0, ubar)
    # Water column thickness
    wct = h+zice
    # Horizontal differentials
    dx, dy = cartesian_grid_2d(lon, lat)
    # Indefinite integral from south to north of u*dz*dy, convert to Sv
    baro_strf = cumsum(ubar*wct*dy, axis=0)*1e-6
    # Colour levels
    lev1 = arange(-50, 150+10, 10)
    # Plot
    ax1 = fig.add_subplot(2, 2, 1, aspect='equal')
    img = contourf(x, y, baro_strf, lev1, extend='both')
    # Contour 0 Sv in black
    contour(x, y, baro_strf, levels=[0], colors=('black'))
    title('Barotropic streamfunction (Sv)', fontsize=24)    
    xlim([-35, 39])
    ylim([-35, 39])
    axis('off')
    cbaxes1 = fig.add_axes([0.07, 0.6, 0.02, 0.3])
    cbar1 = colorbar(img, ticks=arange(-50, 150+50, 50), cax=cbaxes1)
    cbar1.ax.tick_params(labelsize=16)

    # JJA mixed layer depth
    start_month = 6  # Start in June
    end_month = 8  # End in August
    start_day = 1  # First day in June
    next_startday = 1  # First day in September
    end_day = 31  # Last day in August
    prev_endday = 31  # Last day in May
    ndays_season = 92  # Number of days in June+July+August
    id = Dataset(file_path, 'r')
    # Read time axis and get dates
    time_id = id.variables['ocean_time']
    time = num2date(time_id[:], units=time_id.units, calendar=time_id.calendar.lower())
    # Find the last timestep we care about
    end_t = -1  # Missing value flag
    for t in range(size(time)-1, -1, -1):
        if time[t].month == end_month and time[t].day in range(end_day-2, end_day+1):
            end_t = t
            break
        if time[t].month == end_month+1 and time[t].day in range(next_startday, next_startday+2):
            end_t = t
            break
    # Make sure we actually found it
    if end_t == -1:
        print 'Error: ' + file_path + ' does not contain a complete JJA'
        return
    # Find the first timestep we care about
    start_t = -1  # Missing value flag
    for t in range(end_t, -1, -1):
        if time[t].month == start_month-1 and time[t].day in range(prev_endday-1, prev_endday+1):
            start_t = t
            break
        if time[t].month == start_month and time[t].day in range(start_day, start_day+3):
            start_t = t
            break
    # Make sure we found it
    if start_t == -1:
        print 'Error: ' + file_path + ' does not contain a complete JJA'
        return
    # Initialise time-averaged KPP boundary layer depth
    hsbl = ma.empty(shape(lon))
    hsbl[:,:] = 0.0    
    ndays = 0
    # Figure out how many of the 5 days represented in start_t we care about
    if time[start_t].month == start_month and time[start_t].day == start_day+2:
        start_days = 5
    elif time[start_t].month == start_month and time[start_t].day == start_day+1:
        start_days = 4
    elif time[start_t].month == start_month and time[start_t].day == start_day:
        start_days = 3
    elif time[start_t].month == start_month-1 and time[start_t].day == prev_endday:
        start_days = 2
    elif time[start_t].month == start_month-1 and time[start_t].day == prev_endday-1:
        start_days = 1
    else:
        print 'Error: starting index is month ' + str(time[start_t].month) + ', day ' + str(time[start_t].day)
        return
    # Integrate Hsbl weighted by start_days
    hsbl += id.variables['Hsbl'][start_t,:-15,1:]*start_days
    ndays += start_days
    # Between start_t and end_t, we care about all the days
    for t in range(start_t+1, end_t):
        hsbl += id.variables['Hsbl'][t,:-15,1:]*5
        ndays += 5
    # Figure out how many of the 5 days represented in end_t we care about
    if time[end_t].month == end_month+1 and time[end_t].day == next_startday+1:
        end_days = 1
    elif time[end_t].month == end_month+1 and time[end_t].day == next_startday:
        end_days = 2
    elif time[end_t].month == end_month and time[end_t].day == end_day:
        end_days = 3
    elif time[end_t].month == end_month and time[end_t].day == end_day-1:
        end_days = 4
    elif time[end_t].month == end_month and time[end_t].day == end_day-2:
        end_days = 5
    else:
        print 'Error: ending index is month ' + str(time[end_t].month) + ', day ' + str(time[end_t].day)
        return
    # Integrate weighted by end_days
    hsbl += id.variables['Hsbl'][end_t,:-15,1:]*end_days
    ndays += end_days
    if ndays != ndays_season:
        print 'Error: found ' + str(ndays) + ' days instead of ' + str(ndays_season)
        return
    id.close()
    # Convert from integral to average
    hsbl[:,:] = hsbl[:,:]/ndays
    # Mask out ice shelves, change sign, and call it mixed layer depth
    mld = ma.masked_where(zice!=0, -hsbl)
    # Colour levels
    lev2 = arange(0, 300+25, 25)
    # Plot
    ax2 = fig.add_subplot(2, 2, 2, aspect='equal')
    img = contourf(x, y, mld, lev2, extend='both')
    # Contour 100 m in black
    contour(x, y, mld, levels=[100], colors=('black'))
    title('Winter mixed layer depth (m)', fontsize=24)
    xlim([-35, 39])
    ylim([-35, 39])
    axis('off')
    cbaxes2 = fig.add_axes([0.9, 0.6, 0.02, 0.3])
    cbar2 = colorbar(img, ticks=arange(0, 300+100, 100), cax=cbaxes2)
    cbar2.ax.tick_params(labelsize=16)    
    
    # Bottom water temperature
    id = Dataset(file_path, 'r')
    bwtemp = mean(id.variables['temp'][:,0,:-15,1:], axis=0)
    id.close()
    # Mask ice shelves
    bwtemp = ma.masked_where(zice!=0, bwtemp)
    # Colour levels
    lev3 = arange(-2, 2+0.2, 0.2)
    # Plot
    ax3 = fig.add_subplot(2, 2, 3, aspect='equal')
    img = contourf(x, y, bwtemp, lev3, extend='both')
    # Contour 0C in black
    contour(x, y, bwtemp, levels=[0], colors=('black'))
    title(r'Bottom temperature ($^{\circ}$C', fontsize=24)
    xlim([-35, 39])
    ylim([-35, 39])
    axis('off')
    cbaxes3 = fig.add_axes([0.07, 0.1, 0.02, 0.3])
    cbar3 = colorbar(img, ticks=arange(-2, 2+1, 1), cax=cbaxes3)
    cbar3.ax.tick_params(labelsize=16)

    # Bottom water salinity
    id = Dataset(file_path, 'r')
    bwsalt = mean(id.variables['salt'][:,0,:-15,1:], axis=0)
    bwsalt=ma.masked_where(zice!=0, bwsalt)
    id.close()
    lev4 = arange(34.5, 34.8+0.025, 0.025)
    ax4 = fig.add_subplot(2, 2, 4, aspect='equal')
    img = contourf(x, y, bwsalt, lev4, extend='both')
    # Contour 34.65 psu in black
    contour(x, y, bwsalt, levels=[34.65], colors=('black'))
    title('Bottom salinity (psu)', fontsize=24)
    xlim([-35, 39])
    ylim([-35, 39])
    axis('off')
    cbaxes4 = fig.add_axes([0.9, 0.1, 0.02, 0.3])
    cbar4 = colorbar(img, ticks=arange(34.5, 34.8+0.1, 0.1), cax=cbaxes4)
    cbar4.ax.tick_params(labelsize=16)    

    fig.show()
    #fig.savefig('fig1_year3.png')


# Command-line interface
if __name__ == "__main__":

    grid_path = raw_input("Path to ROMS grid file: ")
    file_path = raw_input("Path to ROMS averages file containing one year of 5-day averages, including a continuous June-August period: ")
    holland_fig1 (grid_path, file_path)
    
