from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *

# Plot the sea surface salinity anomaly between the U3 and U3_LIM simulations
# on 23 August (sea ice area maximum).
def adv_u3_sss_anom ():

    # Paths to simulation directories
    paths = ['/short/m68/kaa561/advection/u3_lim/', '/short/m68/kaa561/advection/u3/']
    # File for 23 August daily average
    file_tail = 'iceh.1992-08-23.nc'
    # Bounds on colour scale
    max_anom = 0.2
    tick_anom = 0.1
    # Degrees to radians conversion factor
    deg2rad = pi/180.
    # Centre of missing circle in grid
    lon_c = 50
    lat_c = -83
    # Radius of missing circle
    radius = 10.5
    # Boundary of regular grid to embed circle in 
    circle_bdry = -70+90

    # Read salinity data from U3_LIM simulation; also grid and mask variables
    id = Dataset(paths[0]+file_tail, 'r')
    data0_tmp = id.variables['sss'][0,:350,:]
    lon_tmp = id.variables['TLON'][:350,:]
    lat_tmp = id.variables['TLAT'][:350,:]
    mask_tmp = id.variables['tmask'][:350,:]
    id.close()

    # Read salinity from U3 simulation
    id = Dataset(paths[1]+file_tail, 'r')
    data1_tmp = id.variables['sss'][0,:350,:]
    id.close()
    data_tmp = data1_tmp - data0_tmp

    # Wrap periodic boundary so there isn't a gap in the plot
    lon = ma.empty([size(lon_tmp,0), size(lon_tmp,1)+1])
    lat = ma.empty([size(lat_tmp,0), size(lat_tmp,1)+1])
    mask = ma.empty([size(mask_tmp,0), size(mask_tmp,1)+1])
    data = ma.empty([size(data_tmp,0), size(data_tmp,1)+1])
    lon[:,:-1] = lon_tmp
    lon[:,-1] = lon_tmp[:,0]
    lat[:,:-1] = lat_tmp
    lat[:,-1] = lat_tmp[:,0]
    mask[:,:-1] = mask_tmp
    mask[:,-1] = mask_tmp[:,0]
    data[:,:-1] = data_tmp
    data[:,-1] = data_tmp[:,0]

    # Land mask
    land = ma.masked_where(mask==1, mask)

    # Circumpolar x and y coordinates for plotting
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)
    # Coordinates of centre of missing circle
    x_c = -(lat_c+90)*cos(lon_c*deg2rad+pi/2)
    y_c = (lat_c+90)*sin(lon_c*deg2rad+pi/2)
    # Regular grid to embed missing circle in
    x_reg, y_reg = meshgrid(linspace(-circle_bdry, circle_bdry, num=100), linspace(-circle_bdry, circle_bdry, num=100))
    # Mask everything except the circle out of the regular grid
    land_circle = zeros(shape(x_reg))
    land_circle = ma.masked_where(sqrt((x_reg-x_c)**2 + (y_reg-y_c)**2) > radius, land_circle)

    fig = figure(figsize=(16,12))
    fig.add_subplot(1,1,1, aspect='equal')
    contourf(x, y, land, 1, colors=(('0.6', '0.6', '0.6')))
    contourf(x_reg, y_reg, land_circle, 1, colors=(('0.6', '0.6', '0.6')))
    img = pcolor(x, y, data, vmin=-max_anom, vmax=max_anom, cmap='RdBu_r')
    cbar = colorbar(img, ticks=arange(-max_anom, max_anom+tick_anom, tick_anom))
    cbar.ax.tick_params(labelsize=20)
    title('Sea surface salinity (psu)\nU3 - U3_LIM', fontsize=30)
    axis('off')

    #fig.show()
    fig.savefig('sss_u3_anom.png')


# Command-line interface
if __name__ == "__main__":

    adv_u3_sss_anom()


