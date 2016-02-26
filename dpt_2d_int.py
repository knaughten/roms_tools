from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from calc_z import *

# Calculates zonal transport through each grid cell in the Drake Passage,
# vertically integrates, takes indefinite integral (cumulative sum) over
# latitude, and makes a contour plot of the 2D (latitude vs time) result.
# Input:
# grid_path = path to ROMS grid file
# file_path = path to ROMS ocean history or averages file
def dpt_2d_int (grid_path, file_path):

    # Grid parameters
    theta_s = 0.9
    theta_b = 4.0
    hc = 40
    N = 31
    # Radius of the Earth in m
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Northern boundary of ROMS grid
    nbdry_val = -38

    # Bounds on Drake Passage: edit for new grids
    # i-index of single north-south line to plot (representing a zonal slice);
    # it doesn't really matter which slice of the Drake Passage this is, due
    # to volume conservation
    i_DP = 1175
    # j-indices of the southern tip of South America (j_min) and the northern
    # tip of the Antarctic Peninsula (j_max); make sure these are far enough
    # north/south to be land points, but not so far that they pass through the
    # land and become ocean again (eg Weddell Sea)
    j_min = 210
    j_max = 300

    print 'Reading grid'
    # Read grid variables
    id = Dataset(grid_path, 'r')
    h = id.variables['h'][:,:]
    zice = id.variables['zice'][:,:]
    lon = id.variables['lon_rho'][:,:]
    lat = id.variables['lat_rho'][:,:]
    mask = id.variables['mask_rho'][:,:]
    id.close()

    # Interpolate latitude to the edges of each cell
    s_bdry = lat[0,:]
    middle_lat = 0.5*(lat[0:-1,:] + lat[1:,:])
    n_bdry = lat[-1,:]*0 + nbdry_val
    lat_edges = ma.concatenate((s_bdry[None,:], middle_lat, n_bdry[None,:]))
    # Subtract to get the change in latitude over each cell
    dlat = lat_edges[1:,:] - lat_edges[0:-1,:]

    # Convert from spherical to Cartesian coordinates
    # dy = r*dlat where dlat is converted to radians
    dy_2d = r*dlat*pi/180.0    
    # Copy into a 3D array, same at each depth level
    dy = tile(dy_2d, (N,1,1))

    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    z, sc_r, Cs_r = calc_z(h, zice, lon, lat, theta_s, theta_b, hc, N)
    # We have z at the midpoint of each cell, now find it on the top and
    # bottom edges of each cell
    z_edges = zeros((size(z,0)+1, size(z,1), size(z,2)))
    z_edges[1:-1,:,:] = 0.5*(z[0:-1,:,:] + z[1:,:,:])
    # At surface, z = 0; at bottom, set z to be the same as the midpoint of
    # the deepest cell
    z_edges[0,:,:] = z[0,:,:]
    # Now find dz
    dz = z_edges[1:,:,:] - z_edges[0:-1,:,:]

    # Get 3D land mask
    mask = tile(mask, (N,1,1))
    # Calculate dydz and mask with land mask
    dydz = ma.masked_where(mask==0, dy*dz)
    # Trim to Drake Passage bounds
    dydz_DP = dydz[:,j_min:j_max,i_DP]
    lat_DP = lat[j_min:j_max,i_DP]

    # Read time values and convert from seconds to years
    id = Dataset(file_path, 'r')
    time = id.variables['ocean_time'][:]/(365*24*60*60)

    # Set up an array of dimension time x lat to store transport values
    transport = ma.empty([size(time), j_max-j_min])
    # Calculate transport one timestep at a time
    for t in range(size(time)):

        print 'Processing timestep ' + str(t+1) + ' of '+str(size(time))
        # Read u and interpolate onto the rho-grid
        u = id.variables['u'][t,:,:,:]
        w_bdry_u = 0.5*(u[:,:,0] + u[:,:,-1])
        middle_u = 0.5*(u[:,:,0:-1] + u[:,:,1:])
        e_bdry_u = w_bdry_u[:,:]
        u_rho = ma.concatenate((w_bdry_u[:,:,None], middle_u, e_bdry_u[:,:,None]), axis=2)
        # Trim to Drake Passage bounds
        u_rho_DP = u_rho[:,j_min:j_max,i_DP]
        # Integrate transport over depth and convert to Sv
        transport[t,:] = sum(u_rho_DP*dydz_DP, axis=0)*1e-6

    # Cumulative sum; flip the latitude axis before and after so the sum
    # goes from north to south, not south to north
    transport_cs = fliplr(cumsum(fliplr(transport), axis=1))

    # Plot
    # Bounds are set to 0-150 Sv, adjust as needed
    clf()
    pcolormesh(time, lat_DP, transpose(transport_cs), vmin=0, vmax=150, cmap='jet')
    colorbar(ticks=arange(0, 150+1, 20))
    title('Drake Passage Transport (Sv), indefinite N-S integral')
    xlabel('Years')
    ylabel('Latitude')
    show()
    #savefig('dp_trans_2d_int.png')


# Command-line interface
if __name__ == "__main__":

    grid_path = raw_input('Enter path to grid file: ')
    file_path = raw_input('Enter path to ocean history/averages file: ')
    dpt_2d_int(grid_path, file_path)
    


