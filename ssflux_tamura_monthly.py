from numpy import *
from netCDF4 import Dataset, num2date
from matplotlib.pyplot import *

def ssflux_tamura_monthly (cice_file, month, save=False, fig_name=None):

    tamura_head = '/short/m68/kaa561/tamura_fluxes/Tamura_ssflux_'
    tamura_tail = '_monthly.nc'
    start_day = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    end_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    month_name = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
    deg2rad = pi/180.0
    rho_fw = 1000.0
    rho_sw = 1026.0
    mps_to_cmpday = 8.64e6

    id = Dataset(cice_file, 'r')
    cice_lon_tmp = id.variables['TLON'][:-15,:]
    cice_lat_tmp = id.variables['TLAT'][:-15,:]
    cice_lon = ma.empty([size(cice_lon_tmp,0), size(cice_lon_tmp,1)+1])
    cice_lat = ma.empty([size(cice_lat_tmp,0), size(cice_lat_tmp,1)+1])
    cice_lon[:,:-1] = cice_lon_tmp
    cice_lon[:,-1] = cice_lon_tmp[:,0]
    cice_lat[:,:-1] = cice_lat_tmp
    cice_lat[:,-1] = cice_lat_tmp[:,0]
    time_id = id.variables['time']
    cice_time = num2date(time_id[:], units=time_id.units, calendar='standard')

    next_month = mod(month+1, 12)
    prev_month = mod(month-1, 12)

    end_t = -1
    for t in range(size(cice_time)-1, -1, -1):
        if cice_time[t].month-1 == next_month and cice_time[t].day in range(start_day[next_month], start_day[next_month]+5):
            end_t = t
            break
    if end_t == -1:
        print 'Error: ' + cice_file + ' does not contain a complete ' + month_name[month]
        return

    start_t = -1
    for t in range(end_t, -1, -1):
        if cice_time[t].month-1 == month and cice_time[t].day in range(start_day[month]+1, start_day[month]+6):
            start_t = t
            break
    if start_t == -1:
        print 'Error: ' + cice_file + ' does not contain a complete ' + month_name[month]
        return

    leap_year = False
    cice_year = cice_time[end_t].year
    if month == 11:
        cice_year = cice_time[start_t].year
    if mod(cice_year, 4) == 0:
        leap_year = True
        if mod(cice_year, 100) == 0:
            leap_year = False
            if mod(cice_.year, 400) == 0:
                leap_year = True
    if leap_year:
        end_day[1] = 29

    cice_data_tmp = ma.empty(shape(cice_lon_tmp))
    cice_data_tmp[:,:] = 0.0
    num_days = 0

    if cice_time[start_t].month-1 == month and cice_time[start_t].day == start_day[month]+5:
        start_days = 5
    elif cice_time[start_t].month-1 == month and cice_time[start_t].day == start_day[month]+4:
        start_days = 4
    elif cice_time[start_t].month-1 == month and cice_time[start_t].day == start_day[month]+3:
        start_days = 3
    elif cice_time[start_t].month-1 == month and cice_time[start_t].day == start_day[month]+2:
        start_days = 2
    elif cice_time[start_t].month-1 == month and cice_time[start_t].day == start_day[month]+1:
        start_days = 1
    else:
        print 'Error: starting index is month ' + str(cice_time[start_t].month) + ', day ' + str(cice_time[start_t].day)
        return

    fresh_ai = id.variables['fresh_ai'][start_t,:-15,:]
    fsalt_ai = id.variables['fsalt_ai'][start_t,:-15,:]
    sss = id.variables['sss'][start_t,:-15,:]
    cice_data_tmp += -1/rho_fw*(fresh_ai*rho_sw/mps_to_cmpday - fsalt_ai*1e3/sss)*start_days
    num_days += start_days

    for t in range(start_t+1, end_t):
        fresh_ai = id.variables['fresh_ai'][t,:-15,:]
        fsalt_ai = id.variables['fsalt_ai'][t,:-15,:]
        sss = id.variables['sss'][t,:-15,:]
        cice_data_tmp += -1/rho_fw*(fresh_ai*rho_sw/mps_to_cmpday - fsalt_ai*1e3/sss)*5
        num_days += 5

    if cice_time[end_t].month-1 == next_month and cice_time[end_t].day == start_day[next_month]+4:
        end_days = 1
    elif cice_time[end_t].month-1 == next_month and cice_time[end_t].day == start_day[next_month]+3:
        end_days = 2
    elif cice_time[end_t].month-1 == next_month and cice_time[end_t].day == start_day[next_month]+2:
        end_days = 3
    elif cice_time[end_t].month-1 == next_month and cice_time[end_t].day == start_day[next_month]+1:
        end_days = 4
    elif cice_time[end_t].month-1 == next_month and cice_time[end_t].day == start_day[next_month]:
        end_days = 5
    else:
        print 'Error: ending index is month ' + str(cice_time[end_t].month) + ', day ' + str(cice_time[end_t].day)
        return

    fresh_ai = id.variables['fresh_ai'][end_t,:-15,:]
    fsalt_ai = id.variables['fsalt_ai'][end_t,:-15,:]
    sss = id.variables['sss'][end_t,:-15,:]
    cice_data_tmp += -1/rho_fw*(fresh_ai*rho_sw/mps_to_cmpday - fsalt_ai*1e3/sss)*end_days
    num_days += end_days

    if num_days != end_day[month]:
        print 'Error: found ' + str(num_days) + ' days instead of ' + str(end_day[month])
        return

    id.close()
    cice_data_tmp /= num_days
    cice_data_tmp *= 1e6

    cice_data = ma.empty([size(cice_data_tmp,0), size(cice_data_tmp,1)+1])
    cice_data[:,:-1] = cice_data_tmp
    cice_data[:,-1] = cice_data_tmp[:,0]

    id = Dataset(tamura_head + str(cice_year) + tamura_tail, 'r')
    tamura_lon = id.variables['longitude'][:,:]
    tamura_lat = id.variables['latitude'][:,:]
    tamura_data = id.variables['ssflux'][month,:,:]
    id.close()
    tamura_data = ma.masked_where(isnan(tamura_data), tamura_data)
    tamura_data *= 1e6

    cice_x = -(cice_lat+90)*cos(cice_lon*deg2rad+pi/2)
    cice_y = (cice_lat+90)*sin(cice_lon*deg2rad+pi/2)
    tamura_x = -(tamura_lat+90)*cos(tamura_lon*deg2rad+pi/2)
    tamura_y = (tamura_lat+90)*sin(tamura_lon*deg2rad+pi/2)

    bound = 1.0 #max(amax(abs(cice_data)), amax(abs(tamura_data)))
    lev = linspace(-bound, bound, num=50)
    bdry1 = max(amin(cice_x), amin(tamura_x))
    bdry2 = min(amax(cice_x), amax(tamura_x))
    bdry3 = max(amin(cice_y), amin(tamura_y))
    bdry4 = min(amax(cice_y), amax(tamura_y))

    fig = figure(figsize=(20,9))
    fig.add_subplot(1,2,1, aspect='equal')
    contourf(tamura_x, tamura_y, tamura_data, lev, cmap='RdYlBu_r')
    title('Tamura', fontsize=24)
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    fig.add_subplot(1,2,2, aspect='equal')
    img = contourf(cice_x, cice_y, cice_data, lev, cmap='RdYlBu_r')
    title('CICE', fontsize=24)
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    axis('off')
    cbaxes = fig.add_axes([0.3, 0.04, 0.4, 0.04])
    cbar = colorbar(img, orientation='horizontal', ticks=arange(-1,1+0.5, 0.5), cax=cbaxes, extend='both')
    cbar.ax.tick_params(labelsize=18)
    suptitle(r'Ice-to-ocean salt flux (10$^{-6}$ kg/m$^2$/s, ' + month_name[month] + ' ' + str(cice_year), fontsize=30)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


if __name__ == "__main__":

    cice_file = raw_input("Path to CICE file: ")
    month = int(raw_input("Month number (1-12): "))-1
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    ssflux_tamura_monthly(cice_file, month, save, fig_name)

    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                changes = raw_input("Enter a parameter to change: (1) file path, (2) month number, (3) save/display; or enter to continue: ")
                if len(changes) == 0:
                    break
                else:
                    if int(changes) == 1:
                        cice_file = raw_input("Path to CICE file: ")
                    elif int(changes) == 2:
                        month = int(raw_input("Month number (1-12): "))-1
                    elif int(changes) == 3:
                        save = not save
            if save:
                fig_name = raw_input("File name for figure: ")
            ssflux_tamura_monthly(cice_file, month, save, fig_name)
        else:
            break
            

    

    
        
    
        
