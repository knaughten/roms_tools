from netCDF4 import *
from numpy import *
from matplotlib.pyplot import *

# Plot anomalies between the U3 and U3_LIM simulations in 2 circumpolar plots.
# On the left, annually averaged frazil formation. On the right, effective sea
# ice thickness on 23 August (sea ice area max).
def gadv_frazil_thickness ():

    # Paths to simulation directories
    path_up3l = '/short/m68/kaa561/advection/u3_lim/'
    path_up3 = '/short/m68/kaa561/advection/u3/'
    # Annually averaged CICE file
    file_tail_avg = 'cice/rundir/history/iceh_avg.nc'
    # Daily averaged CICE file for 23 August
    file_tail_23aug = 'cice/rundir/history/iceh.1992-08-23.nc'

    # Bounds for colour scales
    max_frazil = 0.5
    tick_frazil = 0.25
    max_thickness = 1
    tick_thickness = 0.5

    # Degrees to radians conversion factor
    deg2rad = pi/180.
    # Centre of missing circle in grid
    lon_c = 50
    lat_c = -83
    # Radius of missing circle
    radius = 10.5
    # Boundary of regular grid to embed circle in
    circle_bdry = -70+90

    # Labels for longitude around the plot
    # Lat and lon locations of labels
    lon_ticks = array([-120, -60, 60, 120, 180])
    lat_ticks = array([-58, -56, -54, -55.5, -56])
    # Label text
    lon_labels = [r'120$^{\circ}$W', r'60$^{\circ}$W', r'60$^{\circ}$E', r'120$^{\circ}$E', r'180$^{\circ}$']
    # Rotation of text
    lon_rot = [-60, 60, -60, 60, 0]

    # Read U3_LIM frazil data and grid
    id = Dataset(path_up3l + file_tail_avg, 'r')
    frazil_tmp = id.variables['frazil'][0,:275,:]
    lon_tmp = id.variables['TLON'][:275,:]
    lat_tmp = id.variables['TLAT'][:275,:]
    mask_tmp = id.variables['tmask'][:275,:]
    id.close()
    # Wrap periodic boundary so there isn't a gap in the plot
    lon = ma.empty([size(lon_tmp,0), size(lon_tmp,1)+1])
    lat = ma.empty([size(lat_tmp,0), size(lat_tmp,1)+1])
    mask = ma.empty([size(mask_tmp,0), size(mask_tmp,1)+1])
    frazil0 = ma.empty([size(frazil_tmp,0), size(frazil_tmp,1)+1])
    lon[:,:-1] = lon_tmp
    lon[:,-1] = lon_tmp[:,0]
    lat[:,:-1] = lat_tmp
    lat[:,-1] = lat_tmp[:,0]
    mask[:,:-1] = mask_tmp
    mask[:,-1] = mask_tmp[:,0]
    frazil0[:,:-1] = frazil_tmp
    frazil0[:,-1] = frazil_tmp[:,0]

    # Read U3 frazil data
    id = Dataset(path_up3 + file_tail_avg, 'r')
    frazil_tmp = id.variables['frazil'][0,:275,:]
    id.close()
    frazil1 = ma.empty([size(frazil_tmp,0), size(frazil_tmp,1)+1])
    frazil1[:,:-1] = frazil_tmp
    frazil1[:,-1] = frazil_tmp[:,0]
    # Calculate anomaly
    frazil_anom = frazil1 - frazil0

    # Repeat for thickness data
    id = Dataset(path_up3l + file_tail_23aug, 'r')
    thickness_tmp = id.variables['aice'][0,:275,:]*id.variables['hi'][0,:275,:]
    id.close()
    thickness0 = ma.empty([size(thickness_tmp,0), size(thickness_tmp,1)+1])
    thickness0[:,:-1] = thickness_tmp
    thickness0[:,-1] = thickness_tmp[:,0]
    id = Dataset(path_up3 + file_tail_23aug, 'r')
    thickness_tmp = id.variables['aice'][0,:275,:]*id.variables['hi'][0,:275,:]
    id.close()
    thickness1 = ma.empty([size(thickness_tmp,0), size(thickness_tmp,1)+1])
    thickness1[:,:-1] = thickness_tmp
    thickness1[:,-1] = thickness_tmp[:,0]
    thickness_anom = thickness1 - thickness0

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

    # Make the plot
    fig = figure(figsize=(18,8))
    # Frazil
    ax = fig.add_subplot(1, 2, 1, aspect='equal')
    # First shade land
    contourf(x, y, land, 1, colors=(('0.6', '0.6', '0.6')))
    # Fill in missing circle
    contourf(x_reg, y_reg, land_circle, 1, colors=(('0.6', '0.6', '0.6')))
    img1 = pcolor(x, y, frazil_anom, vmin=-max_frazil, vmax=max_frazil, cmap='RdBu_r')
    # Add longitude labels
    for i in range(size(x_ticks)):
        text(x_ticks[i], y_ticks[i], lon_labels[i], ha='center', rotation=lon_rot[i], fontsize=15)
    axis('off')
    title('a) Frazil ice formation (cm/day), annual average', fontsize=20)
    # Colourbar
    cbaxes1 = fig.add_axes([0.07, 0.25, 0.015, 0.5])
    cbar1 = colorbar(img1, ticks=arange(-max_frazil, max_frazil+tick_frazil, tick_frazil), cax=cbaxes1, extend='both')
    cbar1.ax.tick_params(labelsize=16)
    # Thickness
    ax = fig.add_subplot(1, 2, 2, aspect='equal')
    contourf(x, y, land, 1, colors=(('0.6', '0.6', '0.6')))
    contourf(x_reg, y_reg, land_circle, 1, colors=(('0.6', '0.6', '0.6')))
    img2 = pcolor(x, y, thickness_anom, vmin=-max_thickness, vmax=max_thickness, cmap='RdBu_r')
    for i in range(size(x_ticks)):
        text(x_ticks[i], y_ticks[i], lon_labels[i], ha='center', rotation=lon_rot[i], fontsize=15)
    axis('off')
    title('b) Effective sea ice thickness (m), 23 August', fontsize=20)
    cbaxes2 = fig.add_axes([0.93, 0.25, 0.015, 0.5])
    cbar2 = colorbar(img2, ticks=arange(-max_thickness, max_thickness+tick_thickness, tick_thickness), cax=cbaxes2, extend='both')
    cbar2.ax.tick_params(labelsize=16)

    # Main title
    suptitle('Anomalies: UP3 minus UP3L', fontsize=26)
    subplots_adjust(wspace=0.05)

    #fig.show()
    fig.savefig('kn_fig2.png')
    

# Command-line interface
if __name__ == "__main__":

    gadv_frazil_thickness()
    
