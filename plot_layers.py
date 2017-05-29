from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from calc_z import *
from interp_lon_roms import *

def plot_layers (grid_path, lon0, depth_min, Vstretching, theta_s, theta_b, hc, N):

    id = Dataset(grid_path, 'r')
    lon_2d = id.variables['lon_rho'][:,:]
    lat_2d = id.variables['lat_rho'][:,:]
    h = id.variables['h'][:,:]
    zice = id.variables['zice'][:,:]
    mask = id.variables['mask_rho'][:,:]
    id.close()

    z_3d, sc_r, Cs_r = calc_z(h, zice, theta_s, theta_b, hc, N, None, Vstretching)

    if lon0 < 0:
        lon_string = str(int(round(-lon0))) + r'$^{\circ}$W'
    else:
        lon_string = str(int(round(lon0))) + r'$^{\circ}$E'

    if lon0 < 0:
        lon0 += 360

    mask_3d = tile(mask, (N,1,1))

    mask_interp, z, lat = interp_lon_roms(mask_3d, z_3d, lat_2d, lon_2d, lon0)
    layer = zeros(shape(z))
    for k in range(N):        
        layer[k,:] = k+1
    layer = ma.masked_where(mask_interp==0, layer)

    # Choose latitude bounds based on land mask
    mask_sum = sum(mask_interp, axis=0)    
    # Find southernmost and northernmost unmasked j-indices
    edges = ma.flatnotmasked_edges(mask_sum)
    j_min = edges[0]
    j_max = edges[1]
    if j_min == 0:
        # There are ocean points right to the southern boundary
        # Don't do anything special
        lat_min = min(lat[:,j_min])
    else:
        # There is land everywhere at the southern boundary
        # Show the last 2 degrees of this land mask
        lat_min = min(lat[:,j_min]) - 2
    if j_max == size(mask_sum) - 1:
        # There are ocean points right to the northern boundary
        # Don't do anything special
        lat_max = max(lat[:,j_max])
    else:
        # There is land everywhere at the northern boundary
        # Show the first 2 degrees of this land mask
        lat_max = max(lat[:,j_max]) + 2

    #lat_min = -85
    #lat_max = -75

    lev = range(1, N)

    fig = figure(figsize=(18,6)) #(18,12))
    contour(lat, z, layer, lev, colors='k')
    title(lon_string) #(r"ROMS terrain-following levels through the Ross Ice Shelf cavity (180$^{\circ}$E)", fontsize=24)
    xlabel('Latitude')
    ylabel('Depth (m)')
    xlim([lat_min, lat_max])
    ylim([depth_min, 0])
    #text(-82, -200, "Ice shelf", fontsize=24)
    #text(-82, -650, "Sea floor", fontsize=24)
    #fig.savefig('ross_levels.png')
    fig.show()

    if lon0 > 180:
        lon0 -= 360


if __name__ == "__main__":

    grid_path = raw_input("Path to grid file: ")
    lon0 = float(raw_input("Enter longitude (-180 to 180): "))
    depth_min = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
    Vstretching = int(raw_input("Vstretching (2 or 4): "))
    theta_s = float(raw_input("theta_s: "))
    theta_b = float(raw_input("theta_b: "))
    hc = float(raw_input("hc: "))
    N = int(raw_input("N: "))
    plot_layers (grid_path, lon0, depth_min, Vstretching, theta_s, theta_b, hc, N)

    while True:
        repeat = raw_input("Make another plot (y/n)? ")
        if repeat == 'y':
            while True:
                changes = raw_input("Enter a parameter to change: (1) grid file, (2) longitude, (3) deepest depth, (4) Vstretching, (5) theta_s, (6) theta_b, (7) hc, (8) N; or enter to continue: ")
                if len(changes) == 0:
                    break
                else:
                    if int(changes) == 1:
                        grid_path = raw_input("Path to grid file: ")
                    elif int(changes) == 2:
                        lon0 = float(raw_input("Enter longitude (-180 to 180): "))
                    elif int(changes) == 3:
                        depth_min = -1*float(raw_input("Deepest depth to plot (positive, metres): "))
                    elif int(changes) == 4:
                        Vstretching = int(raw_input("Vstretching (2 or 4): "))
                    elif int(changes) == 5:
                        theta_s = float(raw_input("theta_s: "))
                    elif int(changes) == 6:
                        theta_b = float(raw_input("theta_b: "))
                    elif int(changes) == 7:
                        hc = float(raw_input("hc: "))
                    elif int(changes) == 8:
                        N = int(raw_input("N: "))
            plot_layers (grid_path, lon0, depth_min, Vstretching, theta_s, theta_b, hc, N)
        else:
            break
                        
    

    
