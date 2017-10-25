from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from cartesian_grid_3d import *
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from fesom_grid import *

def mip_watermass_barchart (roms_grid, roms_file, fesom_mesh_lr, fesom_mesh_hr, fesom_file_lr, fesom_file_hr):

    # Sectors to consider
    sector_names = ['a) Filchner-Ronne Ice Shelf', 'b) Eastern Weddell Region', 'c) Amery Ice Shelf', 'd) Australian Sector', 'e) Ross Sea', 'f) Amundsen Sea', 'g) Bellingshausen Sea', 'h) Larsen Ice Shelves', 'i) All Ice Shelves']
    num_sectors = len(sector_names)
    # Water masses to consider
    wm_names = ['ISW', 'MCDW', 'HSSW', 'LSSW', 'AASW']
    num_watermasses = len(wm_names)
    wm_colours = [(0.73, 0.6, 1), (1, 0.4, 0.4), (0.52, 0.88, 0.52), (0.6, 0.8, 1), (1, 1, 0)]
    # ROMS vertical grid parameters
    theta_s = 7.0
    theta_b = 2.0
    hc = 250
    N = 31
    # FESOM mesh parameters
    circumpolar = True
    cross_180 = False

    print 'Processing MetROMS'
    # Read ROMS grid variables we need
    id = Dataset(roms_grid, 'r')
    roms_lon = id.variables['lon_rho'][:,:]
    roms_lat = id.variables['lat_rho'][:,:]
    roms_h = id.variables['h'][:,:]
    roms_zice = id.variables['zice'][:,:]
    id.close()
    num_lat = size(roms_lat, 0)
    num_lon = size(roms_lon, 1)
    # Get integrands on 3D grid
    roms_dx, roms_dy, roms_dz, roms_z = cartesian_grid_3d(roms_lon, roms_lat, roms_h, roms_zice, theta_s, theta_b, hc, N)
    # Get volume integrand
    dV = roms_dx*roms_dy*roms_dz
    # Read ROMS output
    id = Dataset(roms_file, 'r')
    roms_temp = id.variables['temp'][0,:,:,:]
    roms_salt = id.variables['salt'][0,:,:,:]
    id.close()
    # Initialise volume of each water mass in each sector
    roms_vol_watermass = zeros([num_watermasses, num_sectors])
    # Calculate water mass breakdown
    for j in range(num_lat):
        for i in range(num_lon):
            # Select ice shelf points
            if roms_zice[j,i] < 0:
                # Figure out which sector this point falls into
                lon = roms_lon[j,i]
                if lon > 180:
                    lon -= 360
                lat = roms_lat[j,i]
                if lon >= -85 and lon < -30 and lat < -74:
                    # Filchner-Ronne
                    sector = 0
                elif lon >= -30 and lon < 65:
                    # Eastern Weddell region
                    sector = 1
                elif lon >= 65 and lon < 76:
                    # Amery
                    sector = 2
                elif lon >= 76 and lon < 165 and lat >= -74:
                    # Australian sector
                    sector = 3
                elif (lon >= 155 and lon < 165 and lat < -74) or (lon >= 165) or (lon < -140):
                    # Ross Sea
                    sector = 4
                elif (lon >= -140 and lon < -105) or (lon >= -105 and lon < -98 and lat < -73.1):
                    # Amundsen Sea
                    sector = 5
                elif (lon >= -104 and lon < -98 and lat >= -73.1) or (lon >= -98 and lon < -66 and lat >= -75):
                    # Bellingshausen Sea
                    sector = 6
                elif lon >= -66 and lon < -59 and lat >= -74:
                    # Larsen Ice Shelves
                    sector = 7
                else:
                    print 'No region found for lon=',str(lon),', lat=',str(lat)
                    break #return
                # Loop downward
                for k in range(N):
                    curr_temp = roms_temp[k,j,i]
                    curr_salt = roms_salt[k,j,i]
                    curr_volume = dV[k,j,i]
                    # Get surface freezing point at this salinity
                    curr_tfrz = curr_salt/(-18.48 + 18.48/1e3*curr_salt)
                    # Figure out what water mass this is
                    if curr_temp < curr_tfrz:
                        # ISW
                        wm_key = 0
                    elif curr_salt < 34:
                        # AASW
                        wm_key = 4
                    elif curr_temp > -1.5:
                        # MCDW
                        wm_key = 1
                    elif curr_salt < 34.5:
                        # LSSW
                        wm_key = 3
                    else:
                        # HSSW
                        wm_key = 2
                    # Integrate volume for the right water mass and sector
                    roms_vol_watermass[wm_key, sector] += curr_volume
                    # Also integrate total Antarctica
                    roms_vol_watermass[wm_key, -1] += curr_volume
    # Find total volume of each sector by adding up the volume of each
    # water mass
    roms_vol_sectors = sum(roms_vol_watermass, axis=0)
    # Calculate percentage of each water mass in each sector
    roms_percent_watermass = zeros([num_watermasses, num_sectors])
    for wm_key in range(num_watermasses):
        for sector in range(num_sectors):
            roms_percent_watermass[wm_key, sector] = roms_vol_watermass[wm_key, sector]/roms_vol_sectors[sector]*100                

    print 'Processing low-res FESOM'
    # Build mesh
    elements_lr = fesom_grid(fesom_mesh_lr, circumpolar, cross_180)
    id = Dataset(fesom_file_lr, 'r')
    temp_nodes_lr = id.variables['temp'][0,:]
    salt_nodes_lr = id.variables['salt'][0,:]
    id.close()
    fesom_vol_watermass_lr = zeros([num_watermasses, num_sectors])
    for i in range(len(elements_lr)):
        elm = elements_lr[i]
        if elm.cavity:
            lon = mean(elm.lon)
            lat = mean(elm.lat)
            if lon >= -85 and lon < -30 and lat < -74:
                sector = 0
            elif lon >= -30 and lon < 65:
                sector = 1
            elif lon >= 65 and lon < 76:
                sector = 2
            elif lon >= 76 and lon < 165 and lat >= -74:
                sector = 3
            elif (lon >= 155 and lon < 165 and lat < -74) or (lon >= 165) or (lon < -140):
                sector = 4
            elif (lon >= -140 and lon < -105) or (lon >= -105 and lon < -98 and lat < -73.1):
                sector = 5
            elif (lon >= -104 and lon < -98 and lat >= -73.1) or (lon >= -98 and lon < -66 and lat >= -75):
                sector = 6
            elif lon >= -66 and lon < -59 and lat >= -74:
                sector = 7
            else:
                print 'No region found for lon=',str(lon),', lat=',str(lat)
                break #return
            # Get area of 2D element
            area = elm.area()
            nodes = [elm.nodes[0], elm.nodes[1], elm.nodes[2]]
            # Loop downward
            while True:
                if nodes[0].below is None or nodes[1].below is None or nodes[2].below is None:
                    # Reached the bottom
                    break
                # Calculate average temperature, salinity, and
                # layer thickness for this 3D triangular prism
                temp_vals = []
                salt_vals = []
                dz_vals = []
                for n in range(3):
                    temp_vals.append(temp_nodes_lr[nodes[n].id])
                    salt_vals.append(salt_nodes_lr[nodes[n].id])
                    temp_vals.append(temp_nodes_lr[nodes[n].below.id])
                    salt_vals.append(salt_nodes_lr[nodes[n].below.id])
                    dz_vals.append(abs(nodes[n].depth - nodes[n].below.depth))
                    # Get ready for next iteration of loop
                    nodes[n] = nodes[n].below
                curr_temp = mean(array(temp_vals))
                curr_salt = mean(array(salt_vals))
                curr_volume = area*mean(array(dz_vals))
                curr_tfrz = -0.0575*curr_salt + 1.7105e-3*sqrt(curr_salt**3) - 2.155e-4*curr_salt**2
                if curr_temp < curr_tfrz:
                    wm_key = 0
                elif curr_salt < 34:
                    wm_key = 4
                elif curr_temp > -1.5:
                    wm_key = 1
                elif curr_salt < 34.5:
                    wm_key = 3
                else:
                    wm_key = 2
                fesom_vol_watermass_lr[wm_key, sector] += curr_volume
                fesom_vol_watermass_lr[wm_key, -1] += curr_volume
    fesom_vol_sectors_lr = sum(fesom_vol_watermass_lr, axis=0)
    fesom_percent_watermass_lr = zeros([num_watermasses, num_sectors])
    for wm_key in range(num_watermasses):
        for sector in range(num_sectors):
            fesom_percent_watermass_lr[wm_key, sector] = fesom_vol_watermass_lr[wm_key, sector]/fesom_vol_sectors_lr[sector]*100

    print 'Processing high-res FESOM'
    elements_hr = fesom_grid(fesom_mesh_hr, circumpolar, cross_180)
    fesom_vol_watermass_hr = zeros([num_watermasses, num_sectors])
    id = Dataset(fesom_file_hr, 'r')
    temp_nodes_hr = id.variables['temp'][0,:]
    salt_nodes_hr = id.variables['salt'][0,:]
    id.close()
    for i in range(len(elements_hr)):
        elm = elements_hr[i]
        if elm.cavity:
            lon = mean(elm.lon)
            lat = mean(elm.lat)
            if lon >= -85 and lon < -30 and lat < -74:
                sector = 0
            elif lon >= -30 and lon < 65:
                sector = 1
            elif lon >= 65 and lon < 76:
                sector = 2
            elif lon >= 76 and lon < 165 and lat >= -74:
                sector = 3
            elif (lon >= 155 and lon < 165 and lat < -74) or (lon >= 165) or (lon < -140):
                sector = 4
            elif (lon >= -140 and lon < -105) or (lon >= -105 and lon < -98 and lat < -73.1):
                sector = 5
            elif (lon >= -104 and lon < -98 and lat >= -73.1) or (lon >= -98 and lon < -66 and lat >= -75):
                sector = 6
            elif lon >= -66 and lon < -59 and lat >= -74:
                sector = 7
            else:
                print 'No region found for lon=',str(lon),', lat=',str(lat)
                break #return
            area = elm.area()
            nodes = [elm.nodes[0], elm.nodes[1], elm.nodes[2]]
            while True:
                if nodes[0].below is None or nodes[1].below is None or nodes[2].below is None:
                    break
                temp_vals = []
                salt_vals = []
                dz_vals = []
                for n in range(3):
                    temp_vals.append(temp_nodes_hr[nodes[n].id])
                    salt_vals.append(salt_nodes_hr[nodes[n].id])
                    temp_vals.append(temp_nodes_hr[nodes[n].below.id])
                    salt_vals.append(salt_nodes_hr[nodes[n].below.id])
                    dz_vals.append(abs(nodes[n].depth - nodes[n].below.depth))
                    nodes[n] = nodes[n].below
                curr_temp = mean(array(temp_vals))
                curr_salt = mean(array(salt_vals))
                curr_volume = area*mean(array(dz_vals))
                curr_tfrz = -0.0575*curr_salt + 1.7105e-3*sqrt(curr_salt**3) - 2.155e-4*curr_salt**2
                if curr_temp < curr_tfrz:
                    wm_key = 0
                elif curr_salt < 34:
                    wm_key = 4
                elif curr_temp > -1.5:
                    wm_key = 1
                elif curr_salt < 34.5:
                    wm_key = 3
                else:
                    wm_key = 2
                fesom_vol_watermass_hr[wm_key, sector] += curr_volume
                fesom_vol_watermass_hr[wm_key, -1] += curr_volume
    fesom_vol_sectors_hr = sum(fesom_vol_watermass_hr, axis=0)
    fesom_percent_watermass_hr = zeros([num_watermasses, num_sectors])
    for wm_key in range(num_watermasses):
        for sector in range(num_sectors):
            fesom_percent_watermass_hr[wm_key, sector] = fesom_vol_watermass_hr[wm_key, sector]/fesom_vol_sectors_hr[sector]*100

    print 'Plotting'
    fig = figure(figsize=(12,6))
    gs = GridSpec(3,3)
    gs.update(left=0.15, right=0.98, bottom=0.2, top=0.88, wspace=0.1, hspace=0.28)
    handles = []
    for sector in range(num_sectors):
        ax = subplot(gs[sector/3, sector%3])
        lefts = 0
        for wm_key in range(num_watermasses):
            ax.barh(0, roms_percent_watermass[wm_key, sector], color=wm_colours[wm_key], left=lefts, align='center')
            lefts += roms_percent_watermass[wm_key, sector]
        lefts = 0
        for wm_key in range(num_watermasses):
            ax.barh(1, fesom_percent_watermass_lr[wm_key, sector], color=wm_colours[wm_key], left=lefts, align='center')
            lefts += fesom_percent_watermass_lr[wm_key, sector]
        lefts = 0
        for wm_key in range(num_watermasses):
            tmp = ax.barh(2, fesom_percent_watermass_hr[wm_key, sector], color=wm_colours[wm_key], left=lefts, align='center')
            if sector == num_sectors-1:
                handles.append(tmp)
            lefts += fesom_percent_watermass_hr[wm_key, sector]
        xlim([0, 100])
        ax.invert_yaxis()
        ax.set_yticks(range(3))
        if sector % 3 == 0:
            ax.set_yticklabels(('MetROMS', 'FESOM (low-res)', 'FESOM (high-res)'))
        else:
            ax.set_yticklabels(('','',''))
        if sector >= num_sectors-3:
            ax.set_xlabel('% volume')
        else:
            ax.set_xticklabels([])
        ax.set_title(sector_names[sector])
    legend(handles, wm_names, ncol=num_watermasses, bbox_to_anchor=(0.35,-0.4))
    subplots_adjust(wspace=0.05, hspace=0.2)
    suptitle('Water masses in ice shelf cavities (2002-2016 average)', fontsize=20)
    fig.show()
    fig.savefig('wm_barchart.png')


# Command-line interface
if __name__ == "__main__":

    roms_grid = raw_input("Path to ROMS grid file: ")
    roms_file = raw_input("Path to ROMS time-averaged temperature and salinity file, 2002-2016: ")
    fesom_mesh_lr = raw_input("Path to FESOM low-res mesh directory: ")
    fesom_file_lr = raw_input("Path to FESOM low-res time-averaged temperature and salinity file, 2002-2016: ")
    fesom_mesh_hr = raw_input("Path to FESOM high-res mesh directory: ")
    fesom_file_hr = raw_input("Path to FESOM high-res time-averaged temperature and salinity file, 2002-2016: ")
    mip_watermass_barchart(roms_grid, roms_file, fesom_mesh_lr, fesom_mesh_hr, fesom_file_lr, fesom_file_hr)

    
    

    
    
    
    
