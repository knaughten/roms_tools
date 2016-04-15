from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import plot,xlabel,ylabel,clf,show

# Calculate the average sea surface height (zeta) at each timestep of the given
# ocean history file.
def avg_zeta (file_path):

    # Read time and grid variables
    file = Dataset(file_path, 'r')
    time = file.variables['ocean_time'][:]
    # Convert time from seconds to years
    time = time/(365*24*60*60)
    lon = file.variables['lon_rho'][:,:]
    lat = file.variables['lat_rho'][:,:]
    mask = file.variables['mask_rho'][:,:]
    avg_zeta = []

    # Mask land points out of lat and lon
    lon = ma.masked_where(mask==0, lon)
    lat = ma.masked_where(mask==0, lat)
    num_lat = size(lon, 0)
    num_lon = size(lon, 1)
    # Add or subtract 360 from longitude values which wrap around
    # so that longitude increases monotonically from west to east
    i = tile(arange(1, num_lon+1), (num_lat, 1))
    index1 = nonzero((i > 1200)*(lon < 100))
    lon[index1] = lon[index1] + 360
    index2 = nonzero((i < 200)*(lon > 300))
    lon[index2] = lon[index2] - 360

    # Interpolate to get longitude at the edges of each cell
    w_bdry = (lon[:,0] + lon[:,num_lon-1] - 360)/2
    middle_lon = (lon[:,0:num_lon-1] + lon[:,1:num_lon])/2
    e_bdry = (lon[:,0] + 360 + lon[:,num_lon-1])/2
    lon_edges = ma.concatenate((w_bdry[:,None], middle_lon, e_bdry[:,None]), axis=1)
    # Subtract to get the change in longitude over each cell
    dlon = abs(lon_edges[:,1:num_lon+1] - lon_edges[:,0:num_lon])

    # Similarly for latitude
    s_bdry = lat[0,:]
    middle_lat = (lat[0:num_lat-1,:] + lat[1:num_lat,:])/2
    n_bdry = lat[num_lat-1,:]*0 - 50
    lat_edges = ma.concatenate((s_bdry[None,:], middle_lat, n_bdry[None,:]))
    dlat = lat_edges[1:num_lat+1,:] - lat_edges[0:num_lat,:]

    # Convert from spherical to Cartesian coordinates
    r = 6.371e6
    # dy = r*dlat where dlat is converted to radians
    dy = r*dlat*pi/180.0
    # dx = r*cos(lat)*dlon where lat and dlon are converted to radians
    dx = r*cos(pi*lat/180.0)*dlon*pi/180.0
    # Multiply to get the area of each cell
    area = nansum(dx*dy)

    for l in range(size(time)):
        print 'Processing timestep ' + str(l+1) + ' of ' + str(size(time))
        # Read zeta at this timestep
        zeta = file.variables['zeta'][l,:,:]
        # Calculate Area-weighted average
        zeta_int = nansum(zeta*dx*dy)
        avg_zeta.append(zeta_int/area)

    file.close()

    # Plot results
    clf()
    plot(time, avg_zeta)
    xlabel('Years')
    ylabel('Average sea surface height (m)')
    show()


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to ocean history file: ")
    avg_zeta(file_path)

    
        
