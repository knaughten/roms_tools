from numpy import *
from netCDF4 import Dataset, num2date
from matplotlib.pyplot import *

def ssflux_tamura_seasonal (cice_file, save=False, fig_name=None):

    tamura_head = '/short/m68/kaa561/tamura_fluxes/Tamura_ssflux_'
    tamura_tail = '_monthly.nc'
    start_month = [12, 3, 6, 9]
    end_month = [2, 5, 8, 11]
    start_day = [1, 1, 1, 1]
    end_day = [28, 31, 31, 30]
    ndays_season = [90, 92, 92, 91]
    ndays_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
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
    cice_time = num2date(time_id[:], units=time_id.units, calendar=time_id.calendar.lower())

    end_t = -1
    for t in range(size(cice_time)-1, -1, -1):
        if cice_time[t].month == start_month[0] and cice_time[t].day in range(start_day[0], start_day[0]+5):
            end_t = t
            break
    if end_t == -1:
        print 'Error: ' + cice_file + ' does not contain a complete Dec-Nov period'
        return

    start_t = -1
    for t in range(end_t-60, -1, -1):
        if cice_time[t].month == start_month[0] and cice_time[t].day in range(start_day[0]+1, start_day[0]+6):
            start_t = t
            break
    if start_t == -1:
        print 'Error: ' + cice_file + ' does not contain a complete Dec-Nov period'
        return

    leap_year = False
    cice_year = cice_time[end_t].year
    if mod(cice_year, 4) == 0:
        leap_year = True
        if mod(cice_year, 100) == 0:
            leap_year = False
            if mod(cice_year, 400) == 0:
                leap_year = True
    if leap_year:
        end_day[0] += 1
        ndays_season[0] += 1

    cice_data_tmp = ma.empty([4, size(cice_lon_tmp,0), size(cice_lon_tmp,1)])
    cice_data_tmp[:,:,:] = 0.0

    for season in range(4):
        season_days = 0
        next_season = mod(season+1, 4)

        start_t_season = -1
        for t in range(start_t, end_t+1):
            if cice_time[t].month == start_month[season] and cice_time[t].day in range(start_day[season]+1, start_day[season]+6):
                start_t_season = t
                break
        if start_t_season == -1:
            print 'Error: could not find starting timestep for season ' + season_names[season]
            return

        end_t_season = -1
        for t in range(start_t_season+1, end_t+1):
            if cice_time[t].month == start_month[next_season] and cice_time[t].day in range(start_day[next_season], start_day[next_season]+5):
                end_t_season = t
                break
        if end_t_season == -1:
            print 'Error: could not find ending timestep for season ' + season_names[season]
            return

        if cice_time[start_t_season].month == start_month[season] and cice_time[start_t_season].day == start_day[season] + 5:
            start_days = 5
        elif cice_time[start_t_season].month == start_month[season] and cice_time[start_t_season].day == start_day[season] + 4:
            start_days = 4
        elif cice_time[start_t_season].month == start_month[season] and cice_time[start_t_season].day == start_day[season]+ 3:
            start_days = 3
        elif cice_time[start_t_season].month == start_month[season] and cice_time[start_t_season].day == start_day[season] + 2:
            start_days = 2
        elif cice_time[start_t_season].month == start_month[season] and cice_time[start_t_season].day == start_day[season] + 1:
            start_days = 1
        else:
            print 'Error for season ' + season_names[season] + ': starting index is month ' + str(cice_time[start_t_season].month) + ', day ' + str(cice_time[start_t_season].day)
            return

        fresh_ai = id.variables['fresh_ai'][start_t_season,:-15,:]
        rain_ai = id.variables['rain_ai'][start_t_season,:-15,:]
        fsalt_ai = id.variables['fsalt_ai'][start_t_season,:-15,:]
        sss = id.variables['sss'][start_t_season,:-15,:]
        cice_data_tmp[season,:,:] += -1/rho_fw*((fresh_ai-rain_ai)*rho_sw/mps_to_cmpday - fsalt_ai*1e3/sss)*start_days
        season_days += start_days

        for t in range(start_t_season+1, end_t_season):
            fresh_ai = id.variables['fresh_ai'][t,:-15,:]
            rain_ai = id.variables['rain_ai'][t,:-15,:]
            fsalt_ai = id.variables['fsalt_ai'][t,:-15,:]
            sss = id.variables['sss'][t,:-15,:]
            cice_data_tmp[season,:,:] += -1/rho_fw*((fresh_ai-rain_ai)*rho_sw/mps_to_cmpday - fsalt_ai*1e3/sss)*5
            season_days += 5

        if cice_time[end_t_season].month == start_month[next_season] and cice_time[end_t_season].day == start_day[next_season] + 4:
            end_days = 1
        elif cice_time[end_t_season].month == start_month[next_season] and cice_time[end_t_season].day == start_day[next_season] + 3:
            end_days = 2
        elif cice_time[end_t_season].month == start_month[next_season] and cice_time[end_t_season].day == start_day[next_season] + 2:
            end_days = 3
        elif cice_time[end_t_season].month == start_month[next_season] and cice_time[end_t_season].day == start_day[next_season] + 1:
            end_days = 4
        elif cice_time[end_t_season].month == start_month[next_season] and cice_time[end_t_season].day == start_day[next_season]:
            end_days = 5
        else:
            print 'Error for season ' + season_names[season] + ': ending index is month ' + str(cice_time[end_t_season].month) + ', day ' + str(cice_time[end_t_season].day)
            return

        fresh_ai = id.variables['fresh_ai'][end_t_season,:-15,:]
        rain_ai = id.variables['rain_ai'][end_t_season,:-15,:]
        fsalt_ai = id.variables['fsalt_ai'][end_t_season,:-15,:]
        sss = id.variables['sss'][end_t_season,:-15,:]
        cice_data_tmp[season,:,:] += -1/rho_fw*((fresh_ai-rain_ai)*rho_sw/mps_to_cmpday - fsalt_ai*1e3/sss)*end_days
        season_days += end_days

        if season_days != ndays_season[season]:
            print 'Error: found ' + str(season_days) + ' days instead of ' + str(ndays_season[season])
            return

        cice_data_tmp[season,:,:] /= season_days

    id.close()
    cice_data_tmp *= 1e6

    cice_data = ma.empty([size(cice_data_tmp,0), size(cice_data_tmp,1), size(cice_data_tmp,2)+1])
    cice_data[:,:,:-1] = cice_data_tmp
    cice_data[:,:,-1] = cice_data_tmp[:,:,0]

    id = Dataset(tamura_head + str(cice_year) + tamura_tail, 'r')
    tamura_lon = id.variables['longitude'][:,:]
    tamura_lat = id.variables['latitude'][:,:]
    id.close()

    tamura_data = ma.empty([4, size(tamura_lon,0), size(tamura_lon,1)])
    tamura_data[:,:,:] = 0.0
    for season in range(4):
        if season == 0:
            tamura_months = [11, 0, 1]
        elif season == 1:
            tamura_months = [2, 3, 4]
        elif season == 2:
            tamura_months = [5, 6, 7]
        elif season == 3:
            tamura_months = [8, 9, 10]
        season_days = 0

        for month in tamura_months:
            if season == 0 and month == 11:
                id = Dataset(tamura_head + str(cice_year-1) + tamura_tail, 'r')
            else:
                id = Dataset(tamura_head + str(cice_year) + tamura_tail, 'r')
            tamura_data_tmp = id.variables['ssflux'][month,:,:]
            id.close()
            tamura_data_tmp = ma.masked_where(isnan(tamura_data_tmp), tamura_data_tmp)
            tamura_data[season,:,:] += tamura_data_tmp*ndays_month[month]
            season_days += ndays_month[month]
        tamura_data[season,:,:] /= season_days

    tamura_data *= 1e6

    cice_x = -(cice_lat+90)*cos(cice_lon*deg2rad+pi/2)
    cice_y = (cice_lat+90)*sin(cice_lon*deg2rad+pi/2)
    tamura_x = -(tamura_lat+90)*cos(tamura_lon*deg2rad+pi/2)
    tamura_y = (tamura_lat+90)*sin(tamura_lon*deg2rad+pi/2)

    lev = linspace(-0.5, 0.5, num=50)
    bdry1 = -35 #max(amin(cice_x), amin(tamura_x))
    bdry2 = 35 #min(amax(cice_x), amax(tamura_x))
    bdry3 = -33 #max(amin(cice_y), amin(tamura_y))
    bdry4 = 37 #min(amax(cice_y), amax(tamura_y))

    fig = figure(figsize=(20,9))
    for season in range(4):
        ax = fig.add_subplot(2, 4, season+1, aspect='equal')
        contourf(tamura_x, tamura_y, tamura_data[season,:,:], lev, cmap='RdYlBu_r', extend='both')
        if season == 0:
            text(-39, 0, 'Tamura', fontsize=24, ha='right')
        title(season_names[season], fontsize=24)
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
        ax = fig.add_subplot(2, 4, season+5, aspect='equal')
        img = contourf(cice_x, cice_y, cice_data[season,:,:], lev, cmap='RdYlBu_r', extend='both')
        if season == 0:
            text(-39, 0, 'CICE', fontsize=24, ha='right')
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
    cbaxes = fig.add_axes([0.25, 0.04, 0.5, 0.02])
    cbar = colorbar(img, orientation='horizontal', ticks=arange(-0.5,0.5+0.25, 0.25), cax=cbaxes)
    cbar.ax.tick_params(labelsize=16)
    suptitle(r'Ice-to-ocean salt flux (10$^{-6}$ kg/m$^2$/s)', fontsize=30)
    subplots_adjust(wspace=0.025,hspace=0.025)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


if __name__ == "__main__":

    cice_file = raw_input("Path to CICE file, containing at least one complete Dec-Nov period: ")
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    ssflux_tamura_seasonal(cice_file, save, fig_name)        
