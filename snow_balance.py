from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *

def snow_balance (cice_file):

    deg2rad = pi/180.0

    id = Dataset(cice_file, 'r')
    lon_tmp = id.variables['TLON'][:-15,:]
    lat_tmp = id.variables['TLAT'][:-15,:]
    hs_i = id.variables['hs'][0,:-15,:]
    hs_f = id.variables['hs'][72,:-15,:]
    total_snow = sum(id.variables['snow_ai'][0:73,:-15,:]*5, axis=0)/100.0
    total_snoice = sum(id.variables['snoice'][0:73,:-15,:]*5, axis=0)/100.0
    total_melts = sum(id.variables['melts'][0:73,:-15,:]*5, axis=0)/100.0
    id.close()

    balance_tmp = hs_i + total_snow - total_snoice - total_melts - hs_f    

    lon = ma.empty([size(lon_tmp,0), size(lon_tmp,1)+1])
    lat = ma.empty([size(lat_tmp,0), size(lat_tmp,1)+1])
    balance = ma.empty([size(balance_tmp,0), size(balance_tmp,1)+1])
    lon[:,:-1] = lon_tmp
    lon[:,-1] = lon_tmp[:,0]
    lat[:,:-1] = lat_tmp
    lat[:,-1] = lat_tmp[:,0]
    balance[:,:-1] = balance_tmp
    balance[:,-1] = balance_tmp[:,0]

    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)
    lev = linspace(-2, 0, num=50)

    fig = figure(figsize=(16,12))
    fig.add_subplot(1,1,1, aspect='equal')
    contourf(x, y, balance, lev, extend='both')
    cbar = colorbar()
    cbar.ax.tick_params(labelsize=20)
    title('Snow balance (m) over first year of simulation\nTotal snow_ai - total melts - total snoice - change in hs', fontsize=30)
    axis('off')

    fig.show()
    #fig.savefig('snow_balance.png')


if __name__ == "__main__":

    cice_file = raw_input("Path to CICE history file, containing at least one year of 5-day averages: ")
    snow_balance(cice_file)


    
