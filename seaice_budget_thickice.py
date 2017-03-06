from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.font_manager import FontProperties
from cartesian_grid_2d import *

def seaice_budget_thickice (cice_file, save=False, fig_names=None):

    threshold = 5.0

    # Read CICE grid
    id = Dataset(cice_file, 'r')
    lon = id.variables['TLON'][:,:]
    lat = id.variables['TLAT'][:,:]
    # Calculate elements of area
    dx, dy = cartesian_grid_2d(lon, lat)
    dA = dx*dy
    # Read time values
    time = id.variables['time'][:]/365.25
    # Read final sea ice thickness
    hi_f = id.variables['hi'][-1,:,:]
    # Read all the fields we need for timeseries
    aice = id.variables['aice'][:,:,:]
    dvidtt = id.variables['dvidtt'][:,:,:]
    dvidtd = id.variables['dvidtd'][:,:,:]
    congel = id.variables['congel'][:,:,:]
    frazil = id.variables['frazil'][:,:,:]
    snoice = id.variables['snoice'][:,:,:]
    meltt = -1*id.variables['meltt'][:,:,:]
    meltb = -1*id.variables['meltb'][:,:,:]
    meltl = -1*id.variables['meltl'][:,:,:]
    id.close()

    flag = hi_f > threshold

    fontP = FontProperties()
    fontP.set_size('small')

    dvidtt_avg = []
    dvidtd_avg = []
    for t in range(size(time)):
        dvidtt_avg.append(sum(dvidtt[t,:,:]*aice[t,:,:]*dA*flag)/sum(aice[t,:,:]*dA*flag))
        dvidtd_avg.append(sum(dvidtd[t,:,:]*aice[t,:,:]*dA*flag)/sum(aice[t,:,:]*dA*flag))
    dvidtt_cum = cumsum(array(dvidtt_avg))*5
    dvidtd_cum = cumsum(array(dvidtd_avg))*5
    dvi_cum = dvidtt_cum + dvidtd_cum

    fig1, ax1 = subplots(figsize=(8,6))
    ax1.plot(time, dvidtt_cum, label='Thermodynamics', color='blue', linewidth=2)
    ax1.plot(time, dvidtd_cum, label='Dynamics', color='green', linewidth=2)
    ax1.plot(time, dvi_cum, label='Total', color='black', linewidth=2)
    title('Cumulative volume tendency averaged over thick ice regions')
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
        congel_avg.append(sum(congel[t,:,:]*aice[t,:,:]*dA*flag)/sum(aice[t,:,:]*dA*flag))
        frazil_avg.append(sum(frazil[t,:,:]*aice[t,:,:]*dA*flag)/sum(aice[t,:,:]*dA*flag))
        snoice_avg.append(sum(snoice[t,:,:]*aice[t,:,:]*dA*flag)/sum(aice[t,:,:]*dA*flag))
        meltt_avg.append(sum(meltt[t,:,:]*aice[t,:,:]*dA*flag)/sum(aice[t,:,:]*dA*flag))
        meltb_avg.append(sum(meltb[t,:,:]*aice[t,:,:]*dA*flag)/sum(aice[t,:,:]*dA*flag))
        meltl_avg.append(sum(meltl[t,:,:]*aice[t,:,:]*dA*flag)/sum(aice[t,:,:]*dA*flag))

    congel_cum = cumsum(array(congel_avg))
    frazil_cum = cumsum(array(frazil_avg))
    snoice_cum = cumsum(array(snoice_avg))
    meltt_cum = cumsum(array(meltt_avg))
    meltb_cum = cumsum(array(meltb_avg))
    meltl_cum = cumsum(array(meltl_avg))
    total_cum = congel_cum + frazil_cum + snoice_cum + meltt_cum + meltb_cum + meltl_cum

    fig2, ax2 = subplots(figsize=(8,6))
    ax2.plot(time, congel_cum, label='Congelation', color='blue', linewidth=2)
    ax2.plot(time, frazil_cum, label='Frazil', color='red', linewidth=2)
    ax2.plot(time, snoice_cum, label='Snow-to-ice', color='cyan', linewidth=2)
    ax2.plot(time, meltt_cum, label='Top melt', color='magenta', linewidth=2)
    ax2.plot(time, meltb_cum, label='Basal melt', color='green', linewidth=2)
    ax2.plot(time, meltl_cum, label='Lateral melt', color='yellow', linewidth=2)
    ax2.plot(time, total_cum, label='Total', color='black', linewidth=2)
    title('Cumulative volume tendency averaged over thick ice regions')
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
    action = raw_input("Save figures (s) or display on screen (d)? ")
    if action == 's':
        save = True
        name1 = raw_input("File name for first figure (thermodynamic versus dynamic): ")
        name2 = raw_input("File name for second figure (thermodynamic terms): ")
        fig_names = [name1, name2]
    else:
        save = False
        fig_names = None
    seaice_budget_thickice(cice_file, save, fig_names)
    

    
    
