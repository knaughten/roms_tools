from netCDF4 import *
from numpy import *
from matplotlib.pyplot import *
#import colormaps as cmaps

# Create a 3x2 plot of annually averaged frazil formation for each advection
# experiment: the absolute values for U3_LIM, and the anomalies from U3_LIM for
# the other 5 experiments.
def adv_frazil ():

    # Paths to simulation directories
    paths = ['/short/m68/kaa561/advection/u3_lim/', '/short/m68/kaa561/advection/u3/', '/short/m68/kaa561/advection/c4_l/', '/short/m68/kaa561/advection/c4_h/', '/short/m68/kaa561/advection/a4_l/', '/short/m68/kaa561/advection/a4_h/']
    # Titles for plotting
    labels = ['a) U3_LIM', 'b) U3 - U3_LIM', 'c) C4_LD - U3_LIM', 'd) C4_HD - U3_LIM', 'e) A4_LD - U3_LIM', 'f) A4_HD - U3_LIM']
    # File name: annual average
    file_tail = 'cice/rundir/history/iceh_avg.nc'

    # Bounds and ticks for colour scales
    max_abs = 0.25 #1.0
    tick_abs = 0.05 #0.25
    max_anom = 0.75 #2.0
    tick_anom = 0.25 #1.0

    # Degrees to radians conversion factor
    deg2rad = pi/180.
    # Centre of missing circle in grid
    lon_c = 50
    lat_c = -83
    # Radius of missing circle
    radius = 10.5
    # Boundary of regular grid to embed circle in
    circle_bdry = -70+90

    lon_ticks = array([-120, -60, 60, 120, 180])
    lat_ticks = array([-44, -42, -42, -44, -41])
    lon_labels = [r'120$^{\circ}$W', r'60$^{\circ}$W', r'60$^{\circ}$E', r'120$^{\circ}$E', r'180$^{\circ}$']
    lon_rot = [-60, 60, -60, 60, 0]

    # Read frazil data from U3_LIM simulation; also grid and mask variables
    id = Dataset(paths[0]+file_tail, 'r')
    data_tmp = id.variables['frazil'][0,:350,:]
    lon_tmp = id.variables['TLON'][:350,:]
    lat_tmp = id.variables['TLAT'][:350,:]
    mask_tmp = id.variables['tmask'][:350,:]
    id.close()

    # Wrap periodic boundary so there isn't a gap in the plot
    lon = ma.empty([size(lon_tmp,0), size(lon_tmp,1)+1])
    lat = ma.empty([size(lat_tmp,0), size(lat_tmp,1)+1])
    mask = ma.empty([size(mask_tmp,0), size(mask_tmp,1)+1])
    data0 = ma.empty([size(data_tmp,0), size(data_tmp,1)+1])
    lon[:,:-1] = lon_tmp
    lon[:,-1] = lon_tmp[:,0]
    lat[:,:-1] = lat_tmp
    lat[:,-1] = lat_tmp[:,0]
    mask[:,:-1] = mask_tmp
    mask[:,-1] = mask_tmp[:,0]
    data0[:,:-1] = data_tmp
    data0[:,-1] = data_tmp[:,0]

    # Land mask
    land = ma.masked_where(mask==1, mask)

    # Circumpolar x and y coordinates for plotting
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)
    # Coordinates of centre of missing circle
    x_c = -(lat_c+90)*cos(lon_c*deg2rad+pi/2)
    y_c = (lat_c+90)*sin(lon_c*deg2rad+pi/2)
    # Longitude labels
    x_ticks = -(lat_ticks+90)*cos(lon_ticks*deg2rad+pi/2)
    y_ticks = (lat_ticks+90)*sin(lon_ticks*deg2rad+pi/2)
    # Regular grid to embed missing circle in
    x_reg, y_reg = meshgrid(linspace(-circle_bdry, circle_bdry, num=100), linspace(-circle_bdry, circle_bdry, num=100))
    # Mask everything except the circle out of the regular grid
    land_circle = zeros(shape(x_reg))
    land_circle = ma.masked_where(sqrt((x_reg-x_c)**2 + (y_reg-y_c)**2) > radius, land_circle)

    # Set up figure
    fig = figure(figsize=(9,15))
    ax = fig.add_subplot(3, 2, 1, aspect='equal')
    # First shade land
    contourf(x, y, land, 1, colors=(('0.6', '0.6', '0.6')))
    # Fill in missing circle
    contourf(x_reg, y_reg, land_circle, 1, colors=(('0.6', '0.6', '0.6')))
    # Shade the frazil data (pcolor not contourf so we don't misrepresent the
    # model grid)
    img0 = pcolor(x, y, data0, vmin=0, vmax=max_abs, cmap='jet')
    # Add longitude labels
    for i in range(size(x_ticks)):
        text(x_ticks[i], y_ticks[i], lon_labels[i], ha='center', rotation=lon_rot[i])
    axis('off')
    # Add title
    title(labels[0], fontsize=20)
    # Add colorbar
    cbaxes0 = fig.add_axes([0.025, 0.7, 0.02, 0.2])
    cbar0 = colorbar(img0, ticks=arange(0, max_abs+tick_abs, tick_abs), cax=cbaxes0, extend='max')
    cbar0.ax.tick_params(labelsize=16)  

    # Loop over the other simulations
    for sim in range(1, len(paths)):
        # Read the frazil data
        id = Dataset(paths[sim]+file_tail, 'r')
        data_tmp = id.variables['frazil'][0,:350,:]
        id.close()
        # Wrap the periodic boundary
        data = ma.empty([size(data_tmp,0), size(data_tmp,1)+1])        
        data[:,:-1] = data_tmp
        data[:,-1] = data_tmp[:,0]
        # Calculate anomaly from U3_LIM
        data = data - data0
        # Add to plot, same as before
        ax = fig.add_subplot(3, 2, sim+1, aspect='equal')
        contourf(x, y, land, 1, colors=(('0.6', '0.6', '0.6')))
        contourf(x_reg, y_reg, land_circle, 1, colors=(('0.6', '0.6', '0.6')))
        img = pcolor(x, y, data, vmin=-max_anom, vmax=max_anom, cmap='RdBu_r')
        axis('off')
        title(labels[sim], fontsize=20)
        if sim == 3:
            # Only add an anomaly colorbar for one of the simulations
            cbaxes = fig.add_axes([0.025, 0.4, 0.02, 0.2])
            cbar = colorbar(img, ticks=arange(-max_anom, max_anom+tick_anom, tick_anom), cax=cbaxes, extend='both')
            cbar.ax.tick_params(labelsize=16)

    # Main title
    suptitle('Annually averaged frazil ice formation (cm/day)', fontsize=24)
    subplots_adjust(wspace=0.025,hspace=0.2)

    #fig.show()
    fig.savefig('adv_frazil.png')


# Command-line interface
if __name__ == '__main__':

    adv_frazil()
    
