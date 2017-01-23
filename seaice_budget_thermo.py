from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.font_manager import FontProperties
from cartesian_grid_2d import *

# Create four plots showing timeseries of the thermodynamic terms of sea ice
# growth and melt: congelation, frazil ice formation, snow-to-ice flooding, 
# top melt, basal melt, and lateral melt. These variables are averaged over
# (1) the continental shelf (defined as anywhere south of 60S with seafloor
# shallower than 1500 m), and (2) the offshore region (everywhere else).
# Two plots are the absolute volume tendencies (units of cm/day), the other
# two are cumulative over the simulation (cm).
# Input:
# cice_file = path to CICE output file containing 5-day averages for the entire
#             simulation
# roms_grid = path to ROMS grid file
# save = optional boolean flag indicating that the plots should be saved to
#        files rather than displayed on the screen
# fig_name = if save=True, an array of size 4 containing filenames for each plot
def seaice_budget_thermo (cice_file, roms_grid, save=False, fig_names=None):

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
    # Read all the fields we need
    aice = id.variables['aice'][:,:,:]
    congel = id.variables['congel'][:,:,:]
    frazil = id.variables['frazil'][:,:,:]
    snoice = id.variables['snoice'][:,:,:]
    meltt = -1*id.variables['meltt'][:,:,:]
    meltb = -1*id.variables['meltb'][:,:,:]
    meltl = -1*id.variables['meltl'][:,:,:]
    id.close()

    # Create masks for shelf and offshore region
    shelf = (lat < -60)*(h < 1500)
    offshore = invert(shelf)

    congel_shelf = []
    frazil_shelf = []
    snoice_shelf = []
    meltt_shelf = []
    meltb_shelf = []
    meltl_shelf = []
    congel_offshore = []
    frazil_offshore = []
    snoice_offshore = []
    meltt_offshore = []
    meltb_offshore = []
    meltl_offshore = []
    # Loop over timesteps
    for t in range(size(time)):
        # Only average over regions with at least 10% sea ice
        aice_flag = aice[t,:,:] > 0.1
        # Congelation averaged over the continental shelf
        congel_shelf.append(sum(congel[t,:,:]*dA*shelf*aice_flag)/sum(dA*shelf*aice_flag))
        # Frazil ice formation averaged over the continental shelf
        frazil_shelf.append(sum(frazil[t,:,:]*dA*shelf*aice_flag)/sum(dA*shelf*aice_flag))
        # Snow-to-ice flooding averaged over the continental shelf
        snoice_shelf.append(sum(snoice[t,:,:]*dA*shelf*aice_flag)/sum(dA*shelf*aice_flag))
        # Top melt averaged over the continental shelf
        meltt_shelf.append(sum(meltt[t,:,:]*dA*shelf*aice_flag)/sum(dA*shelf*aice_flag))
        # Basal melt averaged over the continental shelf
        meltb_shelf.append(sum(meltb[t,:,:]*dA*shelf*aice_flag)/sum(dA*shelf*aice_flag))
        # Lateral melt averaged over the continental shelf
        meltl_shelf.append(sum(meltl[t,:,:]*dA*shelf*aice_flag)/sum(dA*shelf*aice_flag))
        # Congelation averaged over the offshore region
        congel_offshore.append(sum(congel[t,:,:]*dA*offshore*aice_flag)/sum(dA*offshore*aice_flag))
        # Frazil ice formation averaged over the offshore region
        frazil_offshore.append(sum(frazil[t,:,:]*dA*offshore*aice_flag)/sum(dA*offshore*aice_flag))
        # Snow-to-ice flooding averaged over the offshore region
        snoice_offshore.append(sum(snoice[t,:,:]*dA*offshore*aice_flag)/sum(dA*offshore*aice_flag))
        # Top melt averaged over the offshore region
        meltt_offshore.append(sum(meltt[t,:,:]*dA*offshore*aice_flag)/sum(dA*offshore*aice_flag))
        # Basal melt averaged over the offshore region
        meltb_offshore.append(sum(meltb[t,:,:]*dA*offshore*aice_flag)/sum(dA*offshore*aice_flag))
        # Lateral melt averaged over the offshore region
        meltl_offshore.append(sum(meltl[t,:,:]*dA*offshore*aice_flag)/sum(dA*offshore*aice_flag))

    # Convert to arrays and sum to get total volume tendency for each region
    congel_shelf = array(congel_shelf)
    frazil_shelf = array(frazil_shelf)
    snoice_shelf = array(snoice_shelf)
    meltt_shelf = array(meltt_shelf)
    meltb_shelf = array(meltb_shelf)
    meltl_shelf = array(meltl_shelf)
    total_shelf = congel_shelf + frazil_shelf + snoice_shelf + meltt_shelf + meltb_shelf + meltl_shelf
    congel_offshore = array(congel_offshore)
    frazil_offshore = array(frazil_offshore)
    snoice_offshore = array(snoice_offshore)
    meltt_offshore = array(meltt_offshore)
    meltb_offshore = array(meltb_offshore)
    meltl_offshore = array(meltl_offshore)
    total_offshore = congel_offshore + frazil_offshore + snoice_offshore + meltt_offshore + meltb_offshore + meltl_offshore

    # Legends need small font to fit
    fontP = FontProperties()
    fontP.set_size('small')

    # Set up continental shelf plot
    fig1, ax1 = subplots(figsize=(8,6))
    # Add one timeseries at a time
    ax1.plot(time, congel_shelf, label='Congelation', color='blue', linewidth=2)
    ax1.plot(time, frazil_shelf, label='Frazil', color='red', linewidth=2)
    ax1.plot(time, snoice_shelf, label='Snow-to-ice', color='cyan', linewidth=2)
    ax1.plot(time, meltt_shelf, label='Top melt', color='magenta', linewidth=2)
    ax1.plot(time, meltb_shelf, label='Basal melt', color='green', linewidth=2)
    ax1.plot(time, meltl_shelf, label='Lateral melt', color='yellow', linewidth=2)
    ax1.plot(time, total_shelf, label='Total', color='black', linewidth=2)
    # Configure plot
    title('Volume tendency averaged over continental shelf')
    xlabel('Time (years)')
    ylabel('cm/day')
    grid(True)
    # Add a legend
    ax1.legend(loc='upper left', prop=fontP)
    if save:
        fig1.savefig(fig_names[0])
    else:
        fig1.show()

    # Same for offshore plot
    fig2, ax2 = subplots(figsize=(8,6))
    ax2.plot(time, congel_offshore, label='Congelation', color='blue', linewidth=2)
    ax2.plot(time, frazil_offshore, label='Frazil', color='red', linewidth=2)
    ax2.plot(time, snoice_offshore, label='Snow-to-ice', color='cyan', linewidth=2)
    ax2.plot(time, meltt_offshore, label='Top melt', color='magenta', linewidth=2)
    ax2.plot(time, meltb_offshore, label='Basal melt', color='green', linewidth=2)
    ax2.plot(time, meltl_offshore, label='Lateral melt', color='yellow', linewidth=2)
    ax2.plot(time, total_offshore, label='Total', color='black', linewidth=2)
    title('Volume tendency averaged over offshore region')
    xlabel('Time (years)')
    ylabel('cm/day')
    grid(True)
    ax2.legend(loc='lower right', prop=fontP)
    if save:
        fig2.savefig(fig_names[1])
    else:
        fig2.show()

    # Get cumulative sums of each term
    congel_shelf_cum = cumsum(congel_shelf)*5
    frazil_shelf_cum = cumsum(frazil_shelf)*5
    snoice_shelf_cum = cumsum(snoice_shelf)*5
    meltt_shelf_cum = cumsum(meltt_shelf)*5
    meltb_shelf_cum = cumsum(meltb_shelf)*5
    meltl_shelf_cum = cumsum(meltl_shelf)*5
    total_shelf_cum = cumsum(total_shelf)*5
    congel_offshore_cum = cumsum(congel_offshore)*5
    frazil_offshore_cum = cumsum(frazil_offshore)*5
    snoice_offshore_cum = cumsum(snoice_offshore)*5
    meltt_offshore_cum = cumsum(meltt_offshore)*5
    meltb_offshore_cum = cumsum(meltb_offshore)*5
    meltl_offshore_cum = cumsum(meltl_offshore)*5
    total_offshore_cum = cumsum(total_offshore)*5

    # Continental shelf cumulative plot
    fig3, ax3 = subplots(figsize=(8,6))
    ax3.plot(time, congel_shelf_cum, label='Congelation', color='blue', linewidth=2)
    ax3.plot(time, frazil_shelf_cum, label='Frazil', color='red', linewidth=2)
    ax3.plot(time, snoice_shelf_cum, label='Snow-to-ice', color='cyan', linewidth=2)
    ax3.plot(time, meltt_shelf_cum, label='Top melt', color='magenta', linewidth=2)
    ax3.plot(time, meltb_shelf_cum, label='Basal melt', color='green', linewidth=2)
    ax3.plot(time, meltl_shelf_cum, label='Lateral melt', color='yellow', linewidth=2)
    ax3.plot(time, total_shelf_cum, label='Total', color='black', linewidth=2)
    title('Cumulative volume tendency averaged over continental shelf')
    xlabel('Time (years)')
    ylabel('cm')
    grid(True)
    ax3.legend(loc='lower left', prop=fontP)
    if save:
        fig3.savefig(fig_names[2])
    else:
        fig3.show()

    # Offshore cumulative plot
    fig4, ax4 = subplots(figsize=(8,6))
    ax4.plot(time, congel_offshore_cum, label='Congelation', color='blue', linewidth=2)
    ax4.plot(time, frazil_offshore_cum, label='Frazil', color='red', linewidth=2)
    ax4.plot(time, snoice_offshore_cum, label='Snow-to-ice', color='cyan', linewidth=2)
    ax4.plot(time, meltt_offshore_cum, label='Top melt', color='magenta', linewidth=2)
    ax4.plot(time, meltb_offshore_cum, label='Basal melt', color='green', linewidth=2)
    ax4.plot(time, meltl_offshore_cum, label='Lateral melt', color='yellow', linewidth=2)
    ax4.plot(time, total_offshore_cum, label='Total', color='black', linewidth=2)
    title('Cumulative volume tendency averaged over offshore region')
    xlabel('Time (years)')
    ylabel('cm/day')
    grid(True)
    ax4.legend(loc='lower left', prop=fontP)
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
    seaice_budget_thermo(cice_file, roms_grid, save, fig_names)     
    

    

    

    
