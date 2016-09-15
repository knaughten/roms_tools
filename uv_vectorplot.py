from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from rotate_vector_roms import *

# Make a circumpolar Antarctic plot of speed overlaid with velocity vectors at
# the given depth (surface, bottom, or vertically averaged).
# Input:
# grid_path = path to ROMS grid file
# file_path = path to ocean history/averages file
# tstep = timestep in file_path to plot (1-indexed)
# depth_key = integer flag indicating whether to plot the surface velocity (1),
#             the bottom velocity (2), or vertically averaged velocity (3)
# save = optional boolean flag indicating that the plot should be saved to a
#        file rather than displayed on the screen
# fig_name = if save=True, filename for figure
def uv_vectorplot (grid_path, file_path, tstep, depth_key, save=False, fig_name=None):

    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180
    # Side length of blocks to average vectors over (can't plot vector at every
    # single point or the plot will be way too crowded)
    block = 15

    # Read angle from grid file
    grid_id = Dataset(grid_path, 'r')
    angle = grid_id.variables['angle'][:-15,:]
    grid_id.close()
    # Read grid and velocity data
    id = Dataset(file_path, 'r')
    lon = id.variables['lon_rho'][:-15,:-2]
    lat = id.variables['lat_rho'][:-15,:-2]
    if depth_key == 1:
        # Surface u and v
        u = id.variables['u'][tstep-1,-1,:-15,:]
        v = id.variables['v'][tstep-1,-1,:-15,:]
    elif depth_key == 2:
        # Bottom u and v
        u = id.variables['u'][tstep-1,0,:-15,:]
        v = id.variables['v'][tstep-1,0,:-15,:]
    elif depth_key == 3:
        # Vertically averaged u and v
        u = id.variables['ubar'][tstep-1,:-15,:]
        v = id.variables['vbar'][tstep-1,:-15,:]
    id.close()

    # Rotate velocities to lat-lon space
    u_lonlat, v_lonlat = rotate_vector_roms(u, v, angle)
    # Throw away the overlapping periodic boundary
    u_rho = u_lonlat[:,:-2]
    v_rho = v_lonlat[:,:-2]
    # Calculate speed for the background filled contour plot
    speed = sqrt(u_rho**2 + v_rho**2)

    # Calculate X and Y coordinates for plotting circumpolar projection
    X = -(lat+90)*cos(lon*deg2rad+pi/2)
    Y  = (lat+90)*sin(lon*deg2rad+pi/2)

    # Calculate velocity components in spherical coordinate space
    # (just differentiate and rearrange spherical coordinate transformations)
    dlon_dt = u_rho/(r*cos(lat*deg2rad)*deg2rad)
    dlat_dt = v_rho/(r*deg2rad)
    # Calculate velocity components in X-Y space (just differentiate and
    # rearrange equations for X and Y above)
    dX_dt = -dlat_dt*cos(lon*deg2rad+pi/2) + (lat+90)*sin(lon*deg2rad+pi/2)*dlon_dt*deg2rad
    dY_dt = dlat_dt*sin(lon*deg2rad+pi/2) + (lat+90)*cos(lon*deg2rad+pi/2)*dlon_dt*deg2rad

    # Average X, Y, dX_dt, and dY_dt over block x block intervals
    # Calculate number of blocks
    size0 = int(ceil(size(X,0)/float(block)))
    size1 = int(ceil((size(X,1)-1)/float(block)))
    # Set up arrays for averaged fields
    X_block = ma.empty([size0, size1])
    Y_block = ma.empty([size0, size1])
    dX_dt_block = ma.empty([size0, size1])
    dY_dt_block = ma.empty([size0, size1])
    # Set up arrays containing boundary indices
    posn0 = range(0, size(X,0), block)
    posn0.append(size(X,0))
    posn1 = range(0, size(X,1), block)
    posn1.append(size(X,1))
    # Double loop to average each block (can't find a more efficient way to do
    # this)
    for j in range(size0):
        for i in range(size1):
            start0 = posn0[j]
            end0 = posn0[j+1]
            start1 = posn1[i]
            end1 = posn1[i+1]
            X_block[j,i] = mean(X[start0:end0, start1:end1])
            Y_block[j,i] = mean(Y[start0:end0, start1:end1])
            dX_dt_block[j,i] = mean(dX_dt[start0:end0, start1:end1])
            dY_dt_block[j,i] = mean(dY_dt[start0:end0, start1:end1])

    # Make the plot
    fig = figure(figsize=(16,12))
    fig.add_subplot(1,1,1, aspect='equal')
    # Contour speed values at every point
    # Use pastel colour map so overlaid vectors will show up
    contourf(X, Y, speed, 50, cmap='Paired')
    cbar = colorbar()
    cbar.ax.tick_params(labelsize=20)
    # Add vectors for each block
    quiver(X_block, Y_block, dX_dt_block, dY_dt_block, color='black')
    if depth_key == 1:
        title('Surface velocity (m/s)', fontsize=30)
    elif depth_key == 2:
        title('Bottom velocity (m/s)', fontsize=30)
    elif depth_key == 3:
        title('Vertically averaged velocity (m/s)', fontsize=30)
    axis('off')

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    grid_path = raw_input("Path to ROMS grid file: ")
    file_path = raw_input("Path to ocean history/averages file: ")
    tstep = int(raw_input("Timestep number (starting at 1): "))
    # Parse depth information
    depth_type = raw_input("Surface velocity (s), bottom velocity (b), or vertically averaged velocity (a)? ")
    if depth_type == 's':
        depth_key = 1
    elif depth_type == 'b':
        depth_key = 2
    elif depth_type == 'a':
        depth_key = 3
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    # Make the plot
    uv_vectorplot(grid_path, file_path, tstep, depth_key, save, fig_name)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the input paramters; repeat until the user
                # is finished
                changes = raw_input("Enter a parameter to change: (1) file path, (2) timestep number, (3) depth, (4) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New file path
                        file_path = raw_input("Path to ocean history/averages file: ")
                    elif int(changes) == 2:
                        # New timestep
                        tstep = int(raw_input("Timestep number (starting at 1): "))
                    elif int(changes) == 3:
                        # New depth information
                        depth_type = raw_input("Surface velocity (s), bottom velocity (b), or vertically averaged velocity (a)? ")
                        if depth_type == 's':
                            depth_key = 1
                        elif depth_type == 'b':
                            depth_key = 2
                        elif depth_type == 'a':
                            depth_key = 3
                    elif int(changes) == 4:
                        # Change from display to save, or vice versa
                        save = not save
            if save:
                # Get file name for figure
                fig_name = raw_input("File name for figure: ")
            # Make the plot
            uv_vectorplot(grid_path, file_path, tstep, depth_key, save, fig_name)
        else:
            break
                        
