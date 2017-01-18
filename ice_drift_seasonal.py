from numpy import *
from netCDF4 import Dataset, num2date
from matplotlib.pyplot import *
from rotate_vector_cice import *

def ice_drift_seasonal (cice_file, save=False, fig_name=None):

    start_month = [2, 5, 8, 11]
    end_month = [4, 7, 10, 1]
    start_day = [1, 1, 1, 1]
    end_day = [30, 31, 31, 31]
    ndays_season = [89, 92, 92, 92]
    season_names = ['FMA', 'MMJ', 'ASO', 'NDJ']
    figure_order = [1, 2, 4, 3]
    r = 6.371e6
    deg2rad = pi/180.0
    block = 15

    id = Dataset(cice_file, 'r')
    lon_tmp = id.variables['TLON'][:-15,:]
    lat_tmp = id.variables['TLAT'][:-15,:]
    angle_tmp = id.variables['ANGLET'][:-15,:]
    lon = ma.empty([size(lon_tmp,0), size(lon_tmp,1)+1])
    lat = ma.empty([size(lat_tmp,0), size(lat_tmp,1)+1])
    angle = ma.empty([size(angle_tmp,0), size(angle_tmp,1)+1])
    lon[:,:-1] = lon_tmp
    lon[:,-1] = lon_tmp[:,0]
    lat[:,:-1] = lat_tmp
    lat[:,-1] = lat_tmp[:,0]
    angle[:,:-1] = angle_tmp
    angle[:,-1] = angle_tmp[:,0]
    time_id = id.variables['time']
    time = num2date(time_id[:], units=time_id.units, calendar=time_id.calendar.lower())

    end_t = -1
    for t in range(size(time)-1, -1, -1):
        if time[t].month == start_month[0] and time[t].day in range(start_day[0], start_day[0]+5):
            end_t = t
            break
    if end_t == -1:
        print 'Error: ' + cice_file + ' does not contain a complete Feb-Jan period'
        return

    start_t = -1
    for t in range(end_t-60, -1, -1):
        if time[t].month == start_month[0] and time[t].day in range(start_day[0]+1, start_day[0]+6):
            start_t = t
            break
    if start_t == -1:
        print 'Error: ' + cice_file + ' does not contain a complete Feb-Jan period'
        return

    leap_year = False
    if mod(time[start_t].year, 4) == 0:
        leap_year = True
        if mod(time[start_t].year, 100) == 0:
            leap_year = False
            if mod(time[start_t].year, 400) == 0:
                leap_year = True
    if leap_year:
        ndays_season[0] += 1

    aice_tmp = ma.empty([4, size(lon_tmp,0), size(lon_tmp,1)])
    aice_tmp[:,:,:] = 0.0
    uxy_tmp = ma.empty([4, size(lon_tmp,0), size(lon_tmp,1)])
    uxy_tmp[:,:,:] = 0.0
    vxy_tmp = ma.empty([4, size(lon_tmp,0), size(lon_tmp,1)])
    vxy_tmp[:,:,:] = 0.0

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
        uxy_tmp[season,:,:] += id.variables['uvel'][start_t_season,:-15,:]*start_days
        vxy_tmp[season,:,:] += id.variables['vvel'][start_t_season,:-15,:]*start_days
        season_days += start_days

        for t in range(start_t_season+1, end_t_season):
            aice_tmp[season,:,:] += id.variables['aice'][t,:-15,:]*5
            uxy_tmp[season,:,:] += id.variables['uvel'][t,:-15,:]*5
            vxy_tmp[season,:,:] += id.variables['vvel'][t,:-15,:]*5
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
        uxy_tmp[season,:,:] += id.variables['uvel'][end_t_season,:-15,:]*end_days
        vxy_tmp[season,:,:] += id.variables['vvel'][end_t_season,:-15,:]*end_days

        season_days += end_days

        if season_days != ndays_season[season]:
            print 'Error: found ' + str(season_days) + ' days instead of ' + str(ndays_season[season])
            return

        aice_tmp[season,:,:] /= season_days
        uxy_tmp[season,:,:] /= season_days
        vxy_tmp[season,:,:] /= season_days

    id.close()

    aice = ma.empty([size(aice_tmp,0), size(aice_tmp,1), size(aice_tmp,2)+1])
    aice[:,:,:-1] = aice_tmp
    aice[:,:,-1] = aice_tmp[:,:,0]
    u_xy = ma.empty([size(uxy_tmp,0), size(uxy_tmp,1), size(uxy_tmp,2)+1])
    u_xy[:,:,:-1] = uxy_tmp
    u_xy[:,:,-1] = uxy_tmp[:,:,0]
    v_xy = ma.empty([size(vxy_tmp,0), size(vxy_tmp,1), size(vxy_tmp,2)+1])
    v_xy[:,:,:-1] = vxy_tmp
    v_xy[:,:,-1] = vxy_tmp[:,:,0]

    u, v = rotate_vector_cice(u_xy, v_xy, angle)
    speed = sqrt(u**2 + v**2)
    theta = arctan2(v, u)
    theta_circ = theta - lon*deg2rad
    u_circ = speed*cos(theta_circ)
    v_circ = speed*sin(theta_circ)

    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)

    size0 = int(ceil(size(x,0)/float(block)))
    size1 = int(ceil((size(x,1)-1)/float(block)))
    x_block = ma.empty([size0, size1])
    y_block = ma.empty([size0, size1])
    u_circ_block = ma.empty([4, size0, size1])
    v_circ_block = ma.empty([4, size0, size1])
    for season in range(4):
        posn0 = range(0, size(x,0), block)
        posn0.append(size(x,0))
        posn1 = range(0, size(x,1), block)
        posn1.append(size(x,1))
        for j in range(size0):
            for i in range(size1):
                start0 = posn0[j]
                end0 = posn0[j+1]
                start1 = posn1[i]
                end1 = posn1[i+1]
                if season == 0:
                    x_block[j,i] = mean(x[start0:end0, start1:end1])
                    y_block[j,i] = mean(y[start0:end0, start1:end1])
                u_circ_block[season,j,i] = mean(u_circ[season, start0:end0, start1:end1])
                v_circ_block[season,j,i] = mean(v_circ[season, start0:end0, start1:end1])

    lev = linspace(0, 1, num=50)
    bdry1 = -35
    bdry2 = 35
    bdry3 = -33
    bdry4 = 37

    fig = figure(figsize=(16,12))
    for season in range(4):
        ax = fig.add_subplot(2, 2, figure_order[season], aspect='equal')
        img = contourf(x, y, aice[season,:,:], lev, cmap='jet')
        q = quiver(x_block, y_block, u_circ_block[season,:,:], v_circ_block[season,:,:], color='black')
        xlim([bdry1, bdry2])
        ylim([bdry3, bdry4])
        axis('off')
        title(season_names[season], fontsize=24)
        if season == 0:
            cbaxes = fig.add_axes([0.07, 0.6, 0.02, 0.3])
            cbar = colorbar(img, ticks=arange(0,1+0.25,0.25), cax=cbaxes)
            cbar.ax.tick_params(labelsize=16)
        if season == 3:
            quiverkey(q, 0.07, 0.3, 0.2, '20 cm/s', coordinates='figure', fontproperties={'size':16})
    suptitle('Sea ice concentration (1) and velocity (m/s)', fontsize=30)
    subplots_adjust(wspace=0.025,hspace=0.1)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


if __name__ == "__main__":

    cice_file = raw_input("Path to CICE file, containing at least one complete Feb-Jan period: ")
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    ice_drift_seasonal(cice_file, save, fig_name)
