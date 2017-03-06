from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.font_manager import FontProperties
from cartesian_grid_2d import *

def nic_budget_ross_weddell (loc_flag, save=False, fig_name=None):

    nic_dir_head = '/g/data/gh5/access_om_025-CORE_NYF/output'
    output_number = arange(135, 137+1)
    nic_dir_tail = '/ice/HISTORY/'
    nic_file_head = 'iceh.0'
    nic_year_number = output_number-4
    days_per_month = [31,28,31,30,31,30,31,31,30,31,30,31]
    num_time = 12*size(output_number)

    if loc_flag == 'r':
        min_i = 300
        max_i = 600
        min_j = 30
        max_j = 120
    elif loc_flag == 'w':
        min_i = 860
        max_i = 1100
        min_j = 30
        max_j = 140

    # Read CICE grid
    id = Dataset(nic_dir_head + str(output_number[0]) + nic_dir_tail + nic_file_head + str(nic_year_number[0]) + '-01.nc', 'r')
    lon = id.variables['TLON'][:,:]
    lat = id.variables['TLAT'][:,:]
    id.close()

    # Calculate elements of area
    dx, dy = cartesian_grid_2d(lon, lat)
    dA = dx*dy
    dA = dA[min_j:max_j, min_i:max_i]

    time = arange(num_time)/12.0

    fontP = FontProperties()
    fontP.set_size('small')

    aice = ma.empty([size(time), size(dA,0), size(dA,1)])
    congel = ma.empty([size(time), size(dA,0), size(dA,1)])
    frazil = ma.empty([size(time), size(dA,0), size(dA,1)])
    snoice = ma.empty([size(time), size(dA,0), size(dA,1)])
    meltt = ma.empty([size(time), size(dA,0), size(dA,1)])
    meltb = ma.empty([size(time), size(dA,0), size(dA,1)])
    meltl = ma.empty([size(time), size(dA,0), size(dA,1)])

    t = 0
    for year in range(size(output_number)):
        for month in range(12):
            if month+1 < 10:
                filename = nic_dir_head + str(output_number[year]) + nic_dir_tail + nic_file_head + str(nic_year_number[year]) + '-0' + str(month+1) + '.nc'
            else:
                filename = nic_dir_head + str(output_number[year]) + nic_dir_tail + nic_file_head + str(nic_year_number[year]) + '-' + str(month+1) + '.nc'
            id = Dataset(filename, 'r')
            aice[t,:,:] = id.variables['aice'][0,min_j:max_j,min_i:max_i]
            congel[t,:,:] = id.variables['congel'][0,min_j:max_j,min_i:max_i]
            frazil[t,:,:] = id.variables['frazil'][0,min_j:max_j,min_i:max_i]
            snoice[t,:,:] = id.variables['snoice'][0,min_j:max_j,min_i:max_i]
            meltt[t,:,:] = -1*id.variables['meltt'][0,min_j:max_j,min_i:max_i]
            meltb[t,:,:] = -1*id.variables['meltb'][0,min_j:max_j,min_i:max_i]
            meltl[t,:,:] = -1*id.variables['meltl'][0,min_j:max_j,min_i:max_i]
            id.close()
            t += 1
    
    congel_avg = []
    frazil_avg = []
    snoice_avg = []
    meltt_avg = []
    meltb_avg = []
    meltl_avg = []
    for t in range(size(time)):
        month = t % 12
        congel_avg.append(sum(congel[t,:,:]*aice[t,:,:]*dA)/sum(aice[t,:,:]*dA)*days_per_month[month])
        frazil_avg.append(sum(frazil[t,:,:]*aice[t,:,:]*dA)/sum(aice[t,:,:]*dA)*days_per_month[month])
        snoice_avg.append(sum(snoice[t,:,:]*aice[t,:,:]*dA)/sum(aice[t,:,:]*dA)*days_per_month[month])
        meltt_avg.append(sum(meltt[t,:,:]*aice[t,:,:]*dA)/sum(aice[t,:,:]*dA)*days_per_month[month])
        meltb_avg.append(sum(meltb[t,:,:]*aice[t,:,:]*dA)/sum(aice[t,:,:]*dA)*days_per_month[month])
        meltl_avg.append(sum(meltl[t,:,:]*aice[t,:,:]*dA)/sum(aice[t,:,:]*dA)*days_per_month[month])

    congel_cum = cumsum(array(congel_avg))
    frazil_cum = cumsum(array(frazil_avg))
    snoice_cum = cumsum(array(snoice_avg))
    meltt_cum = cumsum(array(meltt_avg))
    meltb_cum = cumsum(array(meltb_avg))
    meltl_cum = cumsum(array(meltl_avg))
    total_cum = congel_cum + frazil_cum + snoice_cum + meltt_cum + meltb_cum + meltl_cum

    fig, ax = subplots(figsize=(8,6))
    ax.plot(time, congel_cum, label='Congelation', color='blue', linewidth=2)
    ax.plot(time, frazil_cum, label='Frazil', color='red', linewidth=2)
    ax.plot(time, snoice_cum, label='Snow-to-ice', color='cyan', linewidth=2)
    ax.plot(time, meltt_cum, label='Top melt', color='magenta', linewidth=2)
    ax.plot(time, meltb_cum, label='Basal melt', color='green', linewidth=2)
    ax.plot(time, meltl_cum, label='Lateral melt', color='yellow', linewidth=2)
    ax.plot(time, total_cum, label='Total', color='black', linewidth=2)
    title_string = 'Cumulative volume tendency averaged over '
    if loc_flag == 'r':
        title_string += 'Ross Sea'
    elif loc_flag == 'w':
        title_string += 'Weddell Sea'
    title(title_string)
    xlabel('Years')
    ylabel('cm')
    grid(True)
    ax.legend(loc='upper left', prop=fontP)
    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


if __name__ == "__main__":

    loc_flag = raw_input("Ross Sea (r) or Weddell Sea (w)? ")
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    else:
        save = False
        fig_name = None
    nic_budget_ross_weddell(loc_flag, save, fig_name)
    

    
    
