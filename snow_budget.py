from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from cartesian_grid_2d import *

def snow_budget (cice_file):

    # Read CICE grid
    id = Dataset(cice_file, 'r')
    lon = id.variables['TLON'][:,:]
    lat = id.variables['TLAT'][:,:]
    # Calculate elements of area
    dx, dy = cartesian_grid_2d(lon, lat)
    dA = dx*dy
    # Read time values
    time = id.variables['time'][:]/365.25
    hs = id.variables['hs'][:,:,:]
    aice = id.variables['aice'][:,:,:]
    snow_ai = id.variables['snow_ai'][:,:,:]/100
    snoice = -1*id.variables['snoice'][:,:,:]/100/3.6
    melts = -1*id.variables['melts'][:,:,:]/100/3.6
    id.close()

    vs = []
    total_snow = []
    total_snoice = []
    total_melts = []
    for t in range(size(time)):
        vs.append(sum(hs[t,:,:]*aice[t,:,:]*dA))
        total_snow.append(sum(snow_ai[t,:,:]*dA))
        total_snoice.append(sum(snoice[t,:,:]*aice[t,:,:]*dA))
        total_melts.append(sum(melts[t,:,:]*aice[t,:,:]*dA))
    delta_vs = vs - vs[0]
    cum_snow = cumsum(array(total_snow))*5
    cum_snoice = cumsum(array(total_snoice))*5
    cum_melts = cumsum(array(total_melts))*5
    total = cum_snow + cum_snoice + cum_melts

    fig, ax = subplots(figsize=(12,9))
    ax.plot(time, delta_vs/1e9, label='Snow volume', color='black', linewidth=2)
    ax.plot(time, cum_snow/1e9, label='Total snowfall', color='red', linewidth=2)
    ax.plot(time, cum_snoice/1e9, label='Total snow-to-ice', color='blue', linewidth=2)
    ax.plot(time, cum_melts/1e9, label='Total snow melt', color='yellow', linewidth=2)
    ax.plot(time, total/1e9, label='Expected snow volume', color='green', linewidth=2)
    xlabel('Time (years)')
    ylabel(r'km$^3$')
    grid(True)
    ax.legend(loc='lower left')
    fig.show()


if __name__ == "__main__":

    cice_file = raw_input("Path to CICE history file: ")
    snow_budget(cice_file)

    

    

    
