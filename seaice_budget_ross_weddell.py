from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.font_manager import FontProperties
from cartesian_grid_2d import *

def seaice_budget_ross_weddell (cice_file, loc_flag, save=False, fig_names=None):

    if loc_flag == 'r':
        min_i = 650
        max_i = 850
        min_j = 80
        max_j = 200
    elif loc_flag == 'w':
        min_i = 1150
        max_i = 1300
        min_j = 50
        max_j = 180

    # Read CICE grid
    id = Dataset(cice_file, 'r')
    lon = id.variables['TLON'][:,:]
    lat = id.variables['TLAT'][:,:]
    # Calculate elements of area
    dx, dy = cartesian_grid_2d(lon, lat)
    dA = dx*dy
    dA = dA[min_j:max_j, min_i:max_i]
    # Read time values
    time = id.variables['time'][:]/365.25
    # Read all the fields we need for timeseries
    aice = id.variables['aice'][:,min_j:max_j,min_i:max_i]
    dvidtt = id.variables['dvidtt'][:,min_j:max_j, min_i:max_i]
    dvidtd = id.variables['dvidtd'][:,min_j:max_j, min_i:max_i]
    congel = id.variables['congel'][:,min_j:max_j, min_i:max_i]
    frazil = id.variables['frazil'][:,min_j:max_j, min_i:max_i]
    snoice = id.variables['snoice'][:,min_j:max_j, min_i:max_i]
    meltt = -1*id.variables['meltt'][:,min_j:max_j, min_i:max_i]
    meltb = -1*id.variables['meltb'][:,min_j:max_j, min_i:max_i]
    meltl = -1*id.variables['meltl'][:,min_j:max_j, min_i:max_i]
    id.close()

    fontP = FontProperties()
    fontP.set_size('small')

    dvidtt_avg = []
    dvidtd_avg = []
    for t in range(size(time)):
        dvidtt_avg.append(sum(dvidtt[t,:,:]*aice[t,:,:]*dA)/sum(aice[t,:,:]*dA))
        dvidtd_avg.append(sum(dvidtd[t,:,:]*aice[t,:,:]*dA)/sum(aice[t,:,:]*dA))
    dvidtt_cum = cumsum(array(dvidtt_avg))*5
    dvidtd_cum = cumsum(array(dvidtd_avg))*5
    dvi_cum = dvidtt_cum + dvidtd_cum

    fig1, ax1 = subplots(figsize=(8,6))
    ax1.plot(time, dvidtt_cum, label='Thermodynamics', color='blue', linewidth=2)
    ax1.plot(time, dvidtd_cum, label='Dynamics', color='green', linewidth=2)
    ax1.plot(time, dvi_cum, label='Total', color='black', linewidth=2)
    title_string = 'Cumulative volume tendency averaged over '
    if loc_flag == 'r':
        title_string += 'Ross Sea'
    elif loc_flag == 'w':
        title_string += 'Weddell Sea'
    title(title_string)
    xlabel('Years')
    ylabel('cm')
    grid(True)
    ax1.legend(loc='upper left', prop=fontP)
    if save:
        fig1.savefig(fig_names[0])
    else:
        fig1.show()

    congel_avg = []
    frazil_avg = []
    snoice_avg = []
    meltt_avg = []
    meltb_avg = []
    meltl_avg = []
    for t in range(size(time)):
        congel_avg.append(sum(congel[t,:,:]*aice[t,:,:]*dA)/sum(aice[t,:,:]*dA))
        frazil_avg.append(sum(frazil[t,:,:]*aice[t,:,:]*dA)/sum(aice[t,:,:]*dA))
        snoice_avg.append(sum(snoice[t,:,:]*aice[t,:,:]*dA)/sum(aice[t,:,:]*dA))
        meltt_avg.append(sum(meltt[t,:,:]*aice[t,:,:]*dA)/sum(aice[t,:,:]*dA))
        meltb_avg.append(sum(meltb[t,:,:]*aice[t,:,:]*dA)/sum(aice[t,:,:]*dA))
        meltl_avg.append(sum(meltl[t,:,:]*aice[t,:,:]*dA)/sum(aice[t,:,:]*dA))

    congel_cum = cumsum(array(congel_avg))*5
    frazil_cum = cumsum(array(frazil_avg))*5
    snoice_cum = cumsum(array(snoice_avg))*5
    meltt_cum = cumsum(array(meltt_avg))*5
    meltb_cum = cumsum(array(meltb_avg))*5
    meltl_cum = cumsum(array(meltl_avg))*5
    total_cum = congel_cum + frazil_cum + snoice_cum + meltt_cum + meltb_cum + meltl_cum

    fig2, ax2 = subplots(figsize=(8,6))
    ax2.plot(time, congel_cum, label='Congelation', color='blue', linewidth=2)
    ax2.plot(time, frazil_cum, label='Frazil', color='red', linewidth=2)
    ax2.plot(time, snoice_cum, label='Snow-to-ice', color='cyan', linewidth=2)
    ax2.plot(time, meltt_cum, label='Top melt', color='magenta', linewidth=2)
    ax2.plot(time, meltb_cum, label='Basal melt', color='green', linewidth=2)
    ax2.plot(time, meltl_cum, label='Lateral melt', color='yellow', linewidth=2)
    ax2.plot(time, total_cum, label='Total', color='black', linewidth=2)
    title_string = 'Cumulative volume tendency averaged over '
    if loc_flag == 'r':
        title_string += 'Ross Sea'
    elif loc_flag == 'w':
        title_string += 'Weddell Sea'
    title(title_string)
    xlabel('Years')
    ylabel('cm')
    grid(True)
    ax2.legend(loc='lower left', prop=fontP)
    if save:
        fig2.savefig(fig_names[1])
    else:
        fig2.show()


if __name__ == "__main__":

    cice_file = raw_input("Path to CICE history file: ")
    loc_flag = raw_input("Ross Sea (r) or Weddell Sea (w)? ")
    action = raw_input("Save figures (s) or display on screen (d)? ")
    if action == 's':
        save = True
        name1 = raw_input("File name for first figure (thermodynamic versus dynamic): ")
        name2 = raw_input("File name for second figure (thermodynamic terms): ")
        fig_names = [name1, name2]
    else:
        save = False
        fig_names = None
    seaice_budget_ross_weddell(cice_file, loc_flag, save, fig_names)
    

    
    
