from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.colors import *
from cartesian_grid_3d import *
from unesco import *

def mip_ts_distribution_ecco2 ():

    # Beginning of ECCO2 filenames
    temp_file_head = '/short/m68/kaa561/metroms_iceshelf/data/originals/ECCO2/THETA.1440x720x50.1992'
    salt_file_head = '/short/m68/kaa561/metroms_iceshelf/data/originals/ECCO2/SALT.1440x720x50.1992'
    # Northern boundary of water masses to consider
    nbdry = -65
    # Number of temperature and salinity bins
    num_bins_temp = 1000
    num_bins_salt = 2000
    # Bounds on temperature and salinity bins (pre-computed, change if needed)
    min_salt = 32.3
    max_salt = 40.1
    min_temp = -3.1
    max_temp = 3.8
    # Bounds to actually plot
    min_salt_plot = 33.25
    max_salt_plot = 35.1
    min_temp_plot = -3
    max_temp_plot = 3.8
    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0

    print 'Setting up bins'
    # Calculate boundaries of temperature bins
    temp_bins = linspace(min_temp, max_temp, num=num_bins_temp)
    # Calculate centres of temperature bins (for plotting)
    temp_centres = 0.5*(temp_bins[:-1] + temp_bins[1:])
    # Repeat for salinity
    salt_bins = linspace(min_salt, max_salt, num=num_bins_salt)
    salt_centres = 0.5*(salt_bins[:-1] + salt_bins[1:])
    # Set up 2D array of temperature bins x salinity bins to hold average
    # depth of water masses, weighted by volume
    ts_vals = zeros([size(temp_centres), size(salt_centres)])
    # Also array to integrate volume
    volume = zeros([size(temp_centres), size(salt_centres)])
    # Calculate surface freezing point as a function of salinity as seen by
    # CICE
    freezing_pt = salt_centres/(-18.48 + 18.48/1e3*salt_centres)
    # Get 2D versions of the temperature and salinity bins
    salt_2d, temp_2d = meshgrid(salt_centres, temp_centres)
    # Calculate potential density of each combination of temperature and
    # salinity bins
    density = unesco(temp_2d, salt_2d, zeros(shape(temp_2d)))-1000
    # Density contours to plot
    density_lev = arange(26.6, 28.4, 0.2)

    print 'Reading grid'
    # Read grid from first file
    id = Dataset(temp_file_head + '01.nc', 'r')
    lon = id.variables['LONGITUDE_T'][:]
    lat = id.variables['LATITUDE_T'][:]
    z = id.variables['DEPTH_T'][:]
    id.close()
    num_lon = size(lon)
    num_lat = size(lat)
    num_depth = size(z)
    # Calculate integrands
    # Interpolate to get longitude at the edges of each cell
    lon_edges = zeros(num_lon+1)
    lon_edges[1:-1] = 0.5*(lon[:-1] + lon[1:])
    lon_edges[0] = 0.5*(lon[0] + lon[-1] - 360)
    lon_edges[-1] = 0.5*(lon[0] + 360 + lon[-1])
    dlon = lon_edges[1:] - lon_edges[:-1]
    # Similarly for latitude; linearly extrapolate for edges (which don't matter)
    lat_edges = zeros(num_lat+1)
    lat_edges[1:-1] = 0.5*(lat[:-1] + lat[1:])
    lat_edges[0] = 2*lat[0] - lat_edges[1]
    lat_edges[-1] = 2*lat[-1] - lat_edges[-2]
    dlat = lat_edges[1:] - lat_edges[:-1]
    # Make 2D versions
    lon_2d, lat_2d = meshgrid(lon, lat)
    dlon_2d, dlat_2d = meshgrid(dlon, dlat)
    # Convert to Cartesian space
    dx_2d = r*cos(lat_2d*deg2rad)*dlon_2d*deg2rad
    dy_2d = r*dlat_2d*deg2rad
    # We have z at the midpoint of each cell, now find it on the top and
    # bottom edges of each cell
    z_edges = zeros(num_depth+1)
    z_edges[1:-1] = 0.5*(z[:-1] + z[1:])
    # At the surface, z=0
    # At bottom, extrapolate
    z_edges[-1] = 2*z[-1] - z_edges[-2]
    # Now find dz
    dz_1d = z_edges[1:] - z_edges[:-1]
    # Tile each array to be 3D
    dx_3d = tile(dx_2d, (num_depth,1,1))
    dy_3d = tile(dy_2d, (num_depth,1,1))
    dz_3d = transpose(tile(dz_1d, (num_lon,num_lat,1)))
    # Get volume integrand
    dV = dx_3d*dy_3d*dz_3d

    print 'Reading data'
    # Annual average over 1992
    temp = ma.empty([num_depth, num_lat, num_lon])
    salt = ma.empty([num_depth, num_lat, num_lon])
    temp[:,:,:] = 0.0
    salt[:,:,:] = 0.0
    for month in range(12):
        if month+1 < 10:
            month_string = '0' + str(month+1)
        else:
            month_string = str(month+1)
        id = Dataset(temp_file_head + month_string + '.nc', 'r')
        temp[:,:,:] += id.variables['THETA'][0,:,:,:]
        id.close()
        id = Dataset(salt_file_head + month_string + '.nc', 'r')
        salt[:,:,:] += id.variables['SALT'][0,:,:,:]
        id.close()
    # Convert from integrals to averages
    temp /= 12.0
    salt /= 12.0

    print 'Binning temperature and salinity'
    # Loop over grid boxes
    # Find the first latitude index north of 65S; stop there
    j_max = nonzero(lat > nbdry)[0][0]
    for k in range(num_depth):
        for j in range(j_max):
            for i in range(num_lon):
                if temp[k,j,i] is ma.masked:
                    # Land
                    continue
                # Figure out which bins this falls into
                temp_index = nonzero(temp_bins > temp[k,j,i])[0][0] - 1
                salt_index = nonzero(salt_bins > salt[k,j,i])[0][0] - 1
                # Integrate depth*dV in this bin
                ts_vals[temp_index, salt_index] += z[k]*dV[k,j,i]
                volume[temp_index, salt_index] += dV[k,j,i]
    # Mask bins with zero volume
    ts_vals = ma.masked_where(volume==0, ts_vals)
    volume = ma.masked_where(volume==0, volume)
    # Convert depths from integrals to volume-averages
    ts_vals /= volume

    # Find the maximum depth for plotting
    max_depth = amax(ts_vals)
    # Make a nonlinear scale
    bounds = linspace(0, max_depth**(1.0/2.5), num=100)**2.5
    norm = BoundaryNorm(boundaries=bounds, ncolors=256)
    # Set labels for density contours
    manual_locations = [(33.4, 3.0), (33.65, 3.0), (33.9, 3.0), (34.2, 3.0), (34.45, 3.5), (34.65, 3.25), (34.9, 3.0), (35, 1.5)]

print "Plotting"
fig = figure(figsize=(9,9))
ax = fig.add_subplot(1, 1, 1)
img = pcolor(salt_centres, temp_centres, ts_vals, norm=norm, vmin=0, vmax=max_depth, cmap='jet')
# Add surface freezing point line
plot(salt_centres, freezing_pt, color='black', linestyle='dashed')
# Add density contours
cs = contour(salt_centres, temp_centres, density, density_lev, colors=(0.6,0.6,0.6), linestyles='dotted')
clabel(cs, inline=1, fontsize=14, color=(0.6,0.6,0.6), fmt='%1.1f', manual=manual_locations)
xlim([min_salt_plot, max_salt_plot])
ylim([min_temp_plot, max_temp_plot])
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)
xlabel('Salinity (psu)', fontsize=16)
ylabel(r'Temperature ($^{\circ}$C)', fontsize=16)
title('Water masses south of 65$^{\circ}$S: depth (m)\n1992 annual average, ECCO2', fontsize=20)
# Add a colourbar on the right
cbaxes = fig.add_axes([0.91, 0.3, 0.02, 0.4])
cbar = colorbar(img, cax=cbaxes, ticks=[0,50,100,200,500,1000,2000,4000])
cbar.ax.tick_params(labelsize=14)
fig.show()
fig.savefig('ts_distribution_ecco2.png')


# Command-line interface
if __name__ == "__main__":

    mip_ts_distribution_ecco2()
    
            
      
                                         

    
