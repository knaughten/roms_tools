from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *

# Saves a contour plot of ice shelf draft, with the land masked in white, and
# the non-ice-shelf-covered ocean masked in grey.
# Input:
# file_path = string containing path to grid file
# fig_name = string containing desired path to output figure (.pdf, .png, or
#            .eps recommended)
def zice (file_path, fig_name):

    # Read zice and masks
    file = Dataset(file_path, 'r')
    zice = file.variables['zice'][:-15,:-3]
    mask_rho = file.variables['mask_rho'][:-15,:-3]
    mask_zice = file.variables['mask_zice'][:-15,:-3]
    file.close()

    # Mask out ocean and land
    zice = -1*ma.masked_where(mask_zice==0, zice)
    ocn = ma.masked_where(mask_rho-mask_zice==0, mask_rho)
    land = ma.masked_where(mask_rho==1, mask_rho)

    # Configure plot
    figure(figsize=(20,8))    
    contourf(zice,50)
    cbar = colorbar(ticks=arange(500,2500,500))
    cbar.ax.tick_params(labelsize=24)
    contourf(land,1,colors=('w'))
    contourf(ocn,1,colors=((0.6,0.6,0.6)))
    title('Ice Shelf Draft (m)', fontsize=32)
    xticks([814,815], (r'0$^{\circ}$',''), fontsize=24)
    yticks([34, 324], (r'75$^{\circ}$S', r'50$^{\circ}$S'), fontsize=24)

    #savefig(fig_name, transparent=True)
    show()


# Command-line interface
if __name__ == "__main__":

    file_path = raw_input("Path to ocean grid file: ")
    fig_name = raw_input("Name of output figure: ")
    zice(file_path, fig_name)
