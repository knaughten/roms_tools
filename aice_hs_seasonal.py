from numpy import *
from netCDF4 import Dataset, num2date
from matplotlib.pyplot import *

def aice_hs_seasonal (cice_file, save=False, fig_name=None):

    start_month = [12, 3, 6, 9]
    end_month = [2, 5, 8, 11]
    start_day = [1, 1, 1, 1]
    end_day = [28, 31, 31, 30]
    ndays_season = [90, 92, 92, 91]
    season_names = ['DJF', 'MAM', 'JJA', 'SON']
    deg2rad = pi/180.0

    id = Dataset(cice_file, 'r')
    lon_tmp = id.variables['TLON'][:-15,:]
    lat_tmp = id.variables['TLAT'][:-15,:]
    lon = ma.empty([size(lon_tmp,0), size(lon_tmp,1)+1])
    lat = ma.empty([size(lat_tmp,0), size(lat_tmp,1)+1])
    lon[:,:-1] = lon_tmp
    lon[:,-1] = lon_tmp[:,0]
    lat[:,:-1] = lat_tmp
    lat[:,-1] = lat_tmp[:,0]
    time_id = id.variables['time']
    time = num2date(time_id[:], units=time_id.units, calendar=time_id.calendar.lower())

    end_t = -1
    for t in range(size(time)-1, -1, -1):
        if time[t].month == start_month[0] and time[t].day in range(start_day[0], start_day[0]+5):
            end_t = t
            break
    if end_t == -1:
        print 'Error: ' + cice_file + ' does not contain a complete Dec-Nov period'
        return

    start_t = -1
    for t in range(end_t-60, -1, -1):
        if time[t].month == start_month[0] and time[t].day in range(start_day[0]+1, start_day[0]+6):
            start_t = t
            break
    if start_t == -1:
        print 'Error: ' + cice_file + ' does not contain a complete Dec-Nov period'
        return

    leap_year = False
    if mod(time[end_t].year, 4) == 0:
        leap_year = True
        if mod(time[end_t].year, 100) == 0:
            leap_year = False
            if mod(time[end_t].year, 400) == 0:
                leap_year = True
    if leap_year:
        end_day[0] += 1
        ndays_season[0] += 1

    aice_tmp = ma.empty([4, size(lon_tmp,0), size(lon_tmp,1)])
    aice_tmp[:,:,:] = 0.0
    hs_tmp = ma.empty([4, size(lon_tmp,0), size(lon_tmp,1)])
    hs_tmp[:,:,:] = 0.0
    for season in range(4):
        season_days = 0
        next_season = mod(season+1, 4)

        start_t_season = -1
        for t in range(start_t, end_t+1):
            if time[t].month == start_month[season] and time[t].day in range(start_day[season]+1, start_day[season]+6):
                start_t_season = t
                break
        if start_t_season == -1:
            print 'Error: could not find starting timestep for season ' + season_names[season]
            return

        end_t_season = -1
        for t in range(start_t_season+1, end_t+1):
            if time[t].month == start_month[next_season] and time[t].day in range(start_day[next_season], start_day[next_season]+5):
                end_t_season = t
                break
        if end_t_season == -1:
            print 'Error: could not find ending timestep for season ' + season_names[season]
            return

        if time[start_t_season].month == start_month[season] and time[start_t_season].day == start_day[season] + 5:
            start_days = 5
        elif time[start_t_season].month == start_month[season] and time[start_t_season].day == start_day[season] + 4:
            start_days = 4
        elif time[start_t_season].month == start_month[season] and time[start_t_season].day == start_day[season]+ 3:
            start_days = 3
        elif time[start_t_season].month == start_month[season] and time[start_t_season].day == start_day[season] + 2:
            start_days = 2
        elif time[start_t_season].month == start_month[season] and time[start_t_season].day == start_day[season] + 1:
            start_days = 1
        else:
            print 'Error for season ' + season_names[season] + ': starting index is month ' + str(time[start_t_season].month) + ', day ' + str(time[start_t_season].day)
            return

        aice_tmp[season,:,:] += id.variables['aice'][start_t_season,:-15,:]*start_days
        hs_tmp[season,:,:] += id.variables['hs'][start_t_season,:-15,:]*start_days
        season_days += start_days

        for t in range(start_t_season+1, end_t_season):
            aice_tmp[season,:,:] += id.variables['aice'][t,:-15,:]*5
            hs_tmp[season,:,:] += id.variables['hs'][t,:-15,:]*5
            season_days += 5

        if time[end_t_season].month == start_month[next_season] and time[end_t_season].day == start_day[next_season] + 4:
            end_days = 1
        elif time[end_t_season].month == start_month[next_season] and time[end_t_season].day == start_day[next_season] + 3:
            end_days = 2
        elif time[end_t_season].month == start_month[next_season] and time[end_t_season].day == start_day[next_season] + 2:
            end_days = 3
        elif time[end_t_season].month == start_month[next_season] and time[end_t_season].day == start_day[next_season] + 1:
            end_days = 4
        elif time[end_t_season].month == start_month[next_season] and time[end_t_season].day == start_day[next_season]:
            end_days = 5
        else:
            print 'Error for season ' + season_names[season] + ': ending index is month ' + str(time[end_t_season].month) + ', day ' + str(time[end_t_season].day)
            return

        aice_tmp[season,:,:] += id.variables['aice'][end_t_season,:-15,:]*end_days
        hs_tmp[season,:,:] += id.variables['hs'][end_t_season,:-15,:]*end_days
        season_days += end_days

        if season_days != ndays_season[season]:
            print 'Error: found ' + str(season_days) + ' days instead of ' + str(ndays_season[season])
            return

        aice_tmp[season,:,:] /= season_days
        hs_tmp[season,:,:] /= season_days

    id.close()

    aice = ma.empty([size(aice_tmp,0), size(aice_tmp,1), size(aice_tmp,2)+1])
    aice[:,:,:-1] = aice_tmp
    aice[:,:,-1] = aice_tmp[:,:,0]
    hs = ma.empty([size(hs_tmp,0), size(hs_tmp,1), size(hs_tmp,2)+1])
    hs[:,:,:-1] = hs_tmp
    hs[:,:,-1] = hs_tmp[:,:,0]

    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)

    # Find boundaries for each side of plot based on extent of NSIDC grid
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

        
        
        

        

        

        

        
        

    
        

    
        

    
    
