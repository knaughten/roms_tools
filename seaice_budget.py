from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from cartesian_grid_2d import *

# Create four plots showing timeseries of the thermodynamic vs dynamic volume
# tendency, averaged over (1) the continental shelf (defined as anywhere south
# of 60S with seafloor shallower than 1500 m), and (2) the offshore region
# (everywhere else). Two plots are the absolute volume tendencies (units of 
# cm/day), the other two are cumulative over the simulation (cm).
# Input:
# cice_file = path to CICE output file containing 5-day averages for the entire
#             simulation
# roms_grid = path to ROMS grid file
# save = optional boolean flag indicating that the plots should be saved to
#        files rather than displayed on the screen
# fig_name = if save=True, an array of size 4 containing filenames for each plot
def seaice_budget (cice_file, roms_grid, save=False, fig_names=None):

    # Read bathymetry values for ROMS grid
    id = Dataset(roms_grid, 'r')
    h = id.variables['h'][1:-1,1:-1]
    id.close()

    # Read CICE grid
    id = Dataset(cice_file, 'r')
    lon = id.variables['TLON'][:,:]
    lat = id.variables['TLAT'][:,:]
    # Calculate elements of area
    dx, dy = cartesian_grid_2d(lon, lat)
    dA = dx*dy
    # Read time values
    time = id.variables['time'][:]/365.25
    # Read data (concentration and thermodynamic/dynamic volume tendencies)
    aice = id.variables['aice'][:,:,:]
    dvidtt = id.variables['dvidtt'][:,:,:]
    dvidtd = id.variables['dvidtd'][:,:,:]
    id.close()

    # Create masks for shelf and offshore region
    shelf = (lat < -60)*(h < 1500)
    offshore = invert(shelf)

    dvidtt_shelf = []
    dvidtd_shelf = []
    dvidtt_offshore = []
    dvidtd_offshore = []
    # Loop over timesteps
    for t in range(size(time)):
        # Only average over regions with at least 10% sea ice
        aice_flag = aice[t,:,:] > 0.1
        # Thermodynamic volume tendency averaged over the continental shelf
        dvidtt_shelf.append(sum(dvidtt[t,:,:]*dA*shelf*aice_flag)/sum(dA*shelf*aice_flag))
        # Dynamic volume tendency averaged over the continental shelf
        dvidtd_shelf.append(sum(dvidtd[t,:,:]*dA*shelf*aice_flag)/sum(dA*shelf*aice_flag))
        # Thermodynamic volume tendency averaged over the offshore region
        dvidtt_offshore.append(sum(dvidtt[t,:,:]*dA*offshore*aice_flag)/sum(dA*offshore*aice_flag))
        # Dynamic volume tendency averaged over the offshore region
        dvidtd_offshore.append(sum(dvidtd[t,:,:]*dA*offshore*aice_flag)/sum(dA*offshore*aice_flag))

    # Convert to arrays and sum to get total volume tendencies for each region
    dvidtt_shelf = array(dvidtt_shelf)
    dvidtd_shelf = array(dvidtd_shelf)
    dvi_shelf = dvidtt_shelf + dvidtd_shelf
    dvidtt_offshore = array(dvidtt_offshore)
    dvidtd_offshore = array(dvidtd_offshore)
    dvi_offshore = dvidtt_offshore + dvidtd_offshore

    # Set up continental shelf plot
    fig1, ax1 = subplots(figsize=(8,6))
    # Add one timeseries at a time
    ax1.plot(time, dvidtt_shelf, label='Thermodynamics', color='blue', linewidth=2)
    ax1.plot(time, dvidtd_shelf, label='Dynamics', color='green', linewidth=2)
    ax1.plot(time, dvi_shelf, label='Total', color='black', linewidth=2)
    # Configure plot
    title('Volume tendency averaged over continental shelf')
    xlabel('Time (years)')
    ylabel('cm/day')
    grid(True)
    # Add a legend
    ax1.legend(loc='upper left')
    if save:
        fig1.savefig(fig_names[0])
    else:
        fig1.show()    

    # Same for offshore plot
    fig2, ax2 = subplots(figsize=(8,6))
    ax2.plot(time, dvidtt_offshore, label='Thermodynamics', color='blue', linewidth=2)
    ax2.plot(time, dvidtd_offshore, label='Dynamics', color='green', linewidth=2)
    ax2.plot(time, dvi_offshore, label='Total', color='black', linewidth=2)
    title('Volume tendency averaged over offshore region')
    xlabel('Time (years)')
    ylabel('cm/day')
    grid(True)
    ax2.legend(loc='lower right')
    if save:
        fig2.savefig(fig_names[1])
    else:
        fig2.show()

    # Get cumulative sums of each term
    dvidtt_shelf_cum = cumsum(dvidtt_shelf)*5
    dvidtd_shelf_cum = cumsum(dvidtd_shelf)*5
    dvi_shelf_cum = cumsum(dvi_shelf)*5
    dvidtt_offshore_cum = cumsum(dvidtt_offshore)*5
    dvidtd_offshore_cum = cumsum(dvidtd_offshore)*5
    dvi_offshore_cum = cumsum(dvi_offshore)*5

    # Continental shelf cumulative plot
    fig3, ax3 = subplots(figsize=(8,6))
    ax3.plot(time, dvidtt_shelf_cum, label='Thermodynamics', color='blue', linewidth=2)
    ax3.plot(time, dvidtd_shelf_cum, label='Dynamics', color='green', linewidth=2)
    ax3.plot(time, dvi_shelf_cum, label='Total', color='black', linewidth=2)
    title('Cumulative volume tendency averaged over continental shelf')
    xlabel('Time (years)')
    ylabel('cm')
    grid(True)
    ax3.legend(loc='upper left')
    if save:
        fig3.savefig(fig_names[2])
    else:
        fig3.show()

    # Offshore cumulative plot
    fig4, ax4 = subplots(figsize=(8,6))
    ax4.plot(time, dvidtt_offshore_cum, label='Thermodynamics', color='blue', linewidth=2)
    ax4.plot(time, dvidtd_offshore_cum, label='Dynamics', color='green', linewidth=2)
    ax4.plot(time, dvi_offshore_cum, label='Total', color='black', linewidth=2)
    title('Cumulative volume tendency averaged over offshore region')
    xlabel('Time (years)')
    ylabel('cm')
    grid(True)
    ax4.legend(loc='upper right')
    if save:
        fig4.savefig(fig_names[3])
    else:
        fig4.show()


# Command-line interface
if __name__ == "__main__":

    cice_file = raw_input("Path to CICE history file: ")
    roms_grid = raw_input("Path to ROMS grid file: ")
    action = raw_input("Save figures (s) or display on screen (d)? ")
    if action == 's':
        save = True
        name1 = raw_input("File name for first figure (continental shelf): ")
        name2 = raw_input("File name for second figure (offshore): ")
        name3 = raw_input("File name for third figure (continental shelf, cumulative): ")
        name4 = raw_input("File name for fourth figure (offshore, cumulative): ")
        fig_names = [name1, name2, name3, name4]
    else:
        save = False
        fig_names = None
    seaice_budget(cice_file, roms_grid, save, fig_names)     
    

    

    

    
