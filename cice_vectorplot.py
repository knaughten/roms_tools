from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from rotate_vector_cice import *

# For any vector in the CICE output (velocities, stresses, etc.) make a
# circumpolar Antarctic plot of its magnitude, overlaid with the vectors
# themselves.
# Input:
# file_path = path to CICE history file
# tstep = timestep in file_path to plot
# xname, yname = names of the x and y vector components in file_path
# cmax = optional maximum magnitude for colour scale
# save = optional boolean flag indicating that the plot should be saved to a
#        file, rather than displayed on the screen
# fig_name = if save=True, filename for figure
def cice_vectorplot (file_path, tstep, xname, yname, cmax=None, save=False, fig_name=None):

    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180
    # Side length of blocks to average vectors over (can't plot vector at
    # every single point or the plot will be way too crowded)
    block = 15

    # Read grid (including rotation angle) and vector components
    id = Dataset(file_path, 'r')
    lon_tmp = id.variables['ULON'][:-15,:]
    lat_tmp = id.variables['ULAT'][:-15,:]
    angle_tmp = id.variables['ANGLET'][:-15,:]
    u_xy_tmp = id.variables[xname][tstep-1,:-15,:]
    v_xy_tmp = id.variables[yname][tstep-1,:-15,:]
    id.close()

    # Wrap periodic boundary by 1 cell
    lon = ma.empty([size(lon_tmp,0), size(lon_tmp,1)+1])
    lat = ma.empty([size(lat_tmp,0), size(lat_tmp,1)+1])
    angle = ma.empty([size(angle_tmp,0), size(angle_tmp,1)+1])
    u_xy = ma.empty([size(u_xy_tmp,0), size(u_xy_tmp,1)+1])
    v_xy = ma.empty([size(v_xy_tmp,0), size(v_xy_tmp,1)+1])
    lon[:,:-1] = lon_tmp
    lon[:,-1] = lon_tmp[:,0]
    lat[:,:-1] = lat_tmp
    lat[:,-1] = lat_tmp[:,0]
    angle[:,:-1] = angle_tmp
    angle[:,-1] = angle_tmp[:,0]
    u_xy[:,:-1] = u_xy_tmp
    u_xy[:,-1] = u_xy_tmp[:,0]
    v_xy[:,:-1] = v_xy_tmp
    v_xy[:,-1] = v_xy_tmp[:,0]

    # Rotate from local x-y space to lon-lat space
    u, v = rotate_vector_cice(u_xy, v_xy, angle)
    # Calculate magnitude for the background filled contour plot
    speed = sqrt(u**2 + v**2)

    # Calculate X and Y coordinates for plotting circumpolar projection
    x = -(lat+90)*cos(lon*deg2rad+pi/2)
    y = (lat+90)*sin(lon*deg2rad+pi/2)

    theta = arctan2(v, u)
    theta_circ = theta - lon*deg2rad
    u_circ = speed*cos(theta_circ)
    v_circ = speed*sin(theta_circ)

    # Average X, Y, dX_dt, and dY_dt over block x block intervals
    # Calculate number of blocks
    size0 = int(ceil(size(x,0)/float(block)))
    size1 = int(ceil((size(x,1)-1)/float(block)))
    # Set up arrays for averaged fields
    x_block = ma.empty([size0, size1])
    y_block = ma.empty([size0, size1])
    u_circ_block = ma.empty([size0, size1])
    v_circ_block = ma.empty([size0, size1])
    # Set up arrays containing boundary indices
    posn0 = range(0, size(x,0), block)
    posn0.append(size(x,0))
    posn1 = range(0, size(x,1), block)
    posn1.append(size(x,1))
    # Double loop to average each block (can't find a more efficient way to do
    # this)
    for j in range(size0):
        for i in range(size1):
            start0 = posn0[j]
            end0 = posn0[j+1]
            start1 = posn1[i]
            end1 = posn1[i+1]
            x_block[j,i] = mean(x[start0:end0, start1:end1])
            y_block[j,i] = mean(y[start0:end0, start1:end1])
            u_circ_block[j,i] = mean(u_circ[start0:end0, start1:end1])
            v_circ_block[j,i] = mean(v_circ[start0:end0, start1:end1])

    # Set up colour scale levels
    if cmax is None:
        lev = linspace(0, amax(speed), num=50)
    else:
        lev = linspace(0, cmax, num=50)

    # Make the plot
    fig = figure(figsize=(16,12))
    fig.add_subplot(1,1,1, aspect='equal')
    # Contour speed values at every point
    # Use pastel colour map so overlaid vectors will show up
    contourf(x, y, speed, lev, cmap='Paired', extend='both')
    cbar = colorbar()
    cbar.ax.tick_params(labelsize=20)
    # Add vectors for each block
    quiver(x_block, y_block, u_circ_block, v_circ_block, color='black')
    title(xname + ', ' + yname, fontsize=30)
    axis('off')

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to CICE history file: ")
    tstep = int(raw_input("Timestep number (starting at 1): "))
    xname = raw_input("Variable name of x-component: ")
    yname = raw_input("Variable name of y-component: ")
    clim = raw_input("Set upper bound on colour scale (y/n)? ")
    cmax = None
    if clim == 'y':
        cmax = float(raw_input("Upper bound: "))
    
    action = raw_input("Save figure (s) or display in window (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("File name for figure: ")
    elif action == 'd':
        save = False
        fig_name = None
    # Make the plot
    cice_vectorplot(file_path, tstep, xname, yname, cmax, save, fig_name)

    # Repeat until the user wants to exit
    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                # Ask for changes to the input parameters; repeat until the user is finished
                changes = raw_input("Enter a parameter to change: (1) file path, (2) timestep number, (3) variable names, (4) colour scale, (5) save/display; or enter to continue: ")
                if len(changes) == 0:
                    # No more changes to parameters
                    break
                else:
                    if int(changes) == 1:
                        # New file path
                        file_path = raw_input("Path to CICE history file: ")
                    elif int(changes) == 2:
                        # New timestep number
                        tstep = int(raw_input("Timestep number (starting at 1): "))
                    elif int(changes) == 3:
                        xname = raw_input("Variable name of x-component: ")
                        yname = raw_input("Variable name of y-component: ")
                    elif int(changes) == 4:
                        cmax = None
                        clim = raw_input("Set upper bound on colour scale (y/n)? ")                        
                        if clim == 'y':
                            cmax = float(raw_input("Upper bound: "))
                    elif int(changes) == 5:
                        # Change from display to save, or vice versa
                        save = not save
            if save:
                # Get file name for figure
                fig_name = raw_input("File name for figure: ")

            # Make the plot
            cice_vectorplot(file_path, tstep, xname, yname, cmax, save, fig_name)

        else:
            break
                

