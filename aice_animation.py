from numpy import *
from netCDF4 import Dataset, num2date
from matplotlib.pyplot import *
from matplotlib.animation import *

# Create an animation of sea ice concentration for the given simulation.
# Save as an mp4 file. (This doesn't run well on the CCRC desktops so use
# something like a MacBook instead.)
# In order for the mp4 saving to work, you must first type "module load ffmpeg"
# on raijin before opening ipython.
# Couldn't work out how to make this an encapsulated function so it is just
# a normal script.

# Directory containing CICE output files
directory = '/short/y99/kaa561/roms_spinup_newest/cice/'
# Number of time indices in each file
num_ts = [18, 270, 270, 270, 252]
# File number to start with for the animation (1-based)
start_file = 5
# Degrees to radians conversion factor
deg2rad = pi/180
# Names of each month for making titles
month_names = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']

# Read grid from the first file
id = Dataset(directory + 'iceh' + str(start_file) + '.nc', 'r')
lon = id.variables['TLON'][:,:]
lat = id.variables['TLAT'][:,:]
id.close()

# Calculate x and y coordinates for polar projection
x = -(lat+90)*cos(lon*deg2rad+pi/2)
y = (lat+90)*sin(lon*deg2rad+pi/2)
# Initialise plot with an array of the correct size, but entirely masked
init_data = ma.empty(shape(x))
init_data[:,:] = ma.masked

# Make the initial figure
fig = figure(figsize=(16,12))
ax = axes(xlim=(amin(x), amax(x)), ylim=(amin(y), amax(y)), aspect='equal')
lev = linspace(0, 1, num=40)
img = ax.contourf(x, y, init_data, lev, cmap='jet', extend='both')
cbar = colorbar(img)
cbar.ax.tick_params(labelsize=20)    

# Animation function to call at each time index
def animate(i):
    # Find the right file and time index to read
    file_num = start_file
    while i >= num_ts[file_num-1]:
        i -= num_ts[file_num-1]
        file_num +=1
    # Read time value and aice data
    id = Dataset(directory + 'iceh' + str(file_num) + '.nc', 'r')
    time_id = id.variables['time']
    time = num2date(time_id[i], units=time_id.units, calendar=time_id.calendar.lower()) 
    data = id.variables['aice'][i,:,:]
    id.close()
    # Clear plot to save memory
    ax.collections = []
    # Plot data
    img = contourf(x, y, data, lev, cmap='jet', extend='both')
    axis('off')
    # Add a title containing the date
    title(str(time.day) + ' ' + month_names[time.month-1], fontsize=30) # + ' ' + str(time.year), fontsize=30)
    return img

# Animate once every time index from start_file to the last file
anim = FuncAnimation(fig, func=animate, frames=range(179,252)) #sum(array(num_ts[start_file-1:])))
# Save as an mp4 with one frame per second
anim.save('aice.mp4', fps=1)
