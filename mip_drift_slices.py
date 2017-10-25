from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.mlab import griddata
from calc_z import *
from interp_lon_roms import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from fesom_grid import *
from fesom_sidegrid import *
from triangle_area import *
from in_triangle import *

# Make a 3x2 plot of temperature (left) and salinity (right) through 0E.
# The top row is the initial conditions from ECCO2. The middle and bottom rows
# are the last January of the simulation (monthly average) from MetROMS and
# FESOM respectively.
# Input:
# roms_grid = path to ROMS grid file
# roms_file = path to file containing Jan 2016 monthly average of temperature
#             and salinity in ROMS
# fesom_mesh_path_lr, fesom_mesh_path_hr = paths to FESOM mesh directories for
#                     low-res and high-res respectively
# fesom_file_lr, fesom_file_hr = paths to files containing Jan 2016 monthly
#                averages of temperature and salinity, in low-res FESOM and
#                high-res FESOM respectively
def mip_drift_slices (roms_grid, roms_file, fesom_mesh_path_lr, fesom_file_lr, fesom_mesh_path_hr, fesom_file_hr):

    # Paths to ECCO2 files with initial conditions for temp and salt
    ecco_temp_file = '/short/m68/kaa561/metroms_iceshelf/data/originals/ECCO2/THETA.1440x720x50.199201.nc'
    ecco_salt_file = '/short/m68/kaa561/metroms_iceshelf/data/originals/ECCO2/SALT.1440x720x50.199201.nc'
    # Longitude to interpolate to (OE)
    lon0 = 0
    # Bounds on plot
    lat_min = -73
    lat_max = -30
    depth_min = -6000
    depth_max = 0
    # ROMS grid parameters
    theta_s = 7.0
    theta_b = 2.0
    hc = 250
    N = 31
    # Bounds on colour scales for temperature and salinity
    temp_min = -2
    temp_max = 6
    salt_min = 33.9
    salt_max = 34.9
    # Contours to overlay
    temp_contour = 0.75
    salt_contour = 34.5
    # Parameters for FESOM regular grid interpolation (needed for contours)
    num_lat = 500
    num_depth = 250
    r = 6.371e6
    deg2rad = pi/180.0

    # Get longitude for the title
    if lon0 < 0:
        lon_string = str(int(round(-lon0))) + r'$^{\circ}$W'
    else:
        lon_string = str(int(round(lon0))) + r'$^{\circ}$E'

    print 'Processing ECCO2'
    id = Dataset(ecco_temp_file, 'r')
    # Read grid variables
    ecco_lat = id.variables['LATITUDE_T'][:]
    ecco_depth = -1*id.variables['DEPTH_T'][:]
    if lon0 == 0:
        # Hard-coded lon0 = 0E: average between the first (0.125 E) and last
        # (359.875 E = -0.125 W) indices in the regular ECCO2 grid
        ecco_temp = 0.5*(id.variables['THETA'][0,:,:,0] + id.variables['THETA'][0,:,:,-1])
        id.close()
        id = Dataset(ecco_salt_file, 'r')
        ecco_salt = 0.5*(id.variables['SALT'][0,:,:,0] + id.variables['SALT'][0,:,:,-1])
        id.close()
    else:
        print 'lon0 is only coded for 0E at this time'
        #return

    print 'Processing ROMS'
    # Read grid variables we need
    id = Dataset(roms_grid, 'r')
    roms_lon_2d = id.variables['lon_rho'][:,:]
    roms_lat_2d = id.variables['lat_rho'][:,:]
    roms_h = id.variables['h'][:,:]
    roms_zice = id.variables['zice'][:,:]
    id.close()
    # Read temperature and salinity
    id = Dataset(roms_file, 'r')
    roms_temp_3d = id.variables['temp'][0,:,:,:]
    roms_salt_3d = id.variables['salt'][0,:,:,:]
    id.close()
    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    roms_z_3d, sc_r, Cs_r = calc_z(roms_h, roms_zice, theta_s, theta_b, hc, N)
    # Make sure we are in the range 0-360
    if lon0 < 0:
        lon0 += 360
    # Interpolate to lon0
    roms_temp, roms_z, roms_lat = interp_lon_roms(roms_temp_3d, roms_z_3d, roms_lat_2d, roms_lon_2d, lon0)
    roms_salt, roms_z, roms_lat = interp_lon_roms(roms_salt_3d, roms_z_3d, roms_lat_2d, roms_lon_2d, lon0)
    # Switch back to range -180-180
    if lon0 > 180:
        lon0 -= 360

    print 'Processing low-res FESOM'
    # Build regular elements
    elements_lr = fesom_grid(fesom_mesh_path_lr)
    # Read temperature and salinity
    id = Dataset(fesom_file_lr, 'r')
    fesom_temp_nodes_lr = id.variables['temp'][0,:]
    fesom_salt_nodes_lr = id.variables['salt'][0,:]
    id.close()
    # Make SideElements
    selements_temp_lr = fesom_sidegrid(elements_lr, fesom_temp_nodes_lr, lon0, lat_max)
    selements_salt_lr = fesom_sidegrid(elements_lr, fesom_salt_nodes_lr, lon0, lat_max)
    # Build an array of quadrilateral patches for the plot, and of data values
    # corresponding to each SideElement
    patches_lr = []
    fesom_temp_lr = []
    for selm in selements_temp_lr:
        # Make patch
        coord = transpose(vstack((selm.y, selm.z)))
        patches_lr.append(Polygon(coord, True, linewidth=0.))
        # Save data value
        fesom_temp_lr.append(selm.var)
    # Repeat for salinity
    fesom_salt_lr = []
    for selm in selements_salt_lr:
        fesom_salt_lr.append(selm.var)
    # Interpolate to regular grid so we can overlay contours
    lat_reg = linspace(lat_min, lat_max, num_lat)
    depth_reg = linspace(-depth_max, -depth_min, num_depth)
    temp_reg_lr = zeros([num_depth, num_lat])
    salt_reg_lr = zeros([num_depth, num_lat])
    temp_reg_lr[:,:] = NaN
    salt_reg_lr[:,:] = NaN
    # For each element, check if a point on the regular grid lies
    # within. If so, do barycentric interpolation to that point, at each
    # depth on the regular grid.
    for elm in elements_lr:
        # Check if this element crosses lon0
        if amin(elm.lon) < lon0 and amax(elm.lon) > lon0:
            # Check if we are within the latitude bounds
            if amax(elm.lat) > lat_min and amin(elm.lat) < lat_max:
                # Find largest regular latitude value south of Element
                tmp = nonzero(lat_reg > amin(elm.lat))[0]
                if len(tmp) == 0:
                    # Element crosses the southern boundary
                    jS = 0
                else:
                    jS = tmp[0] - 1
                # Find smallest regular latitude north of Element
                tmp = nonzero(lat_reg > amax(elm.lat))[0]
                if len(tmp) == 0:
                    # Element crosses the northern boundary
                    jN = num_lat
                else:
                    jN = tmp[0]
                for j in range(jS+1,jN):
                    # There is a chance that the regular gridpoint at j
                    # lies within this element
                    lat0 = lat_reg[j]
                    if in_triangle(elm, lon0, lat0):
                        # Yes it does
                        # Get area of entire triangle
                        area = triangle_area(elm.lon, elm.lat)
                        # Get area of each sub-triangle formed by (lon0, lat0)
                        area0 = triangle_area([lon0, elm.lon[1], elm.lon[2]], [lat0, elm.lat[1], elm.lat[2]])
                        area1 = triangle_area([lon0, elm.lon[0], elm.lon[2]], [lat0, elm.lat[0], elm.lat[2]])
                        area2 = triangle_area([lon0, elm.lon[0], elm.lon[1]], [lat0, elm.lat[0], elm.lat[1]])
                        # Find fractional area of each
                        cff = [area0/area, area1/area, area2/area]
                        # Interpolate each depth value
                        for k in range(num_depth):
                            # Linear interpolation in the vertical for the
                            # value at each corner of the triangle
                            node_vals_temp = []
                            node_vals_salt = []
                            for n in range(3):
                                id1, id2, coeff1, coeff2 = elm.nodes[n].find_depth(depth_reg[k])
                                if any(isnan(array([id1, id2, coeff1, coeff2]))):
                                    # No ocean data here (seafloor or ice shelf)
                                    node_vals_temp.append(NaN)
                                    node_vals_salt.append(NaN)
                                else:
                                    node_vals_temp.append(coeff1*fesom_temp_nodes_lr[id1] + coeff2*fesom_temp_nodes_lr[id2])
                                    node_vals_salt.append(coeff1*fesom_salt_nodes_lr[id1] + coeff2*fesom_salt_nodes_lr[id2])
                            if any(isnan(node_vals_temp)):
                                pass
                            else:
                                # Barycentric interpolation for the value at
                                # lon0, lat0
                                temp_reg_lr[k,j] = sum(array(cff)*array(node_vals_temp))
                                salt_reg_lr[k,j] = sum(array(cff)*array(node_vals_salt))
    temp_reg_lr = ma.masked_where(isnan(temp_reg_lr), temp_reg_lr)
    salt_reg_lr = ma.masked_where(isnan(salt_reg_lr), salt_reg_lr)

    print 'Processing high-res FESOM'
    elements_hr = fesom_grid(fesom_mesh_path_hr)
    id = Dataset(fesom_file_hr, 'r')
    fesom_temp_nodes_hr = id.variables['temp'][0,:]
    fesom_salt_nodes_hr = id.variables['salt'][0,:]
    id.close()
    selements_temp_hr = fesom_sidegrid(elements_hr, fesom_temp_nodes_hr, lon0, lat_max)
    selements_salt_hr = fesom_sidegrid(elements_hr, fesom_salt_nodes_hr, lon0, lat_max)
    patches_hr = []
    fesom_temp_hr = []
    for selm in selements_temp_hr:
        coord = transpose(vstack((selm.y, selm.z)))
        patches_hr.append(Polygon(coord, True, linewidth=0.))
        fesom_temp_hr.append(selm.var)
    fesom_salt_hr = []
    for selm in selements_salt_hr:
        fesom_salt_hr.append(selm.var)
    lat_reg = linspace(lat_min, lat_max, num_lat)
    temp_reg_hr = zeros([num_depth, num_lat])
    salt_reg_hr = zeros([num_depth, num_lat])
    temp_reg_hr[:,:] = NaN
    salt_reg_hr[:,:] = NaN
    for elm in elements_hr:
        if amin(elm.lon) < lon0 and amax(elm.lon) > lon0:
            if amax(elm.lat) > lat_min and amin(elm.lat) < lat_max:
                tmp = nonzero(lat_reg > amin(elm.lat))[0]
                if len(tmp) == 0:
                    jS = 0
                else:
                    jS = tmp[0] - 1
                tmp = nonzero(lat_reg > amax(elm.lat))[0]
                if len(tmp) == 0:
                    jN = num_lat
                else:
                    jN = tmp[0]
                for j in range(jS+1,jN):
                    lat0 = lat_reg[j]
                    if in_triangle(elm, lon0, lat0):
                        area = triangle_area(elm.lon, elm.lat)
                        area0 = triangle_area([lon0, elm.lon[1], elm.lon[2]], [lat0, elm.lat[1], elm.lat[2]])
                        area1 = triangle_area([lon0, elm.lon[0], elm.lon[2]], [lat0, elm.lat[0], elm.lat[2]])
                        area2 = triangle_area([lon0, elm.lon[0], elm.lon[1]], [lat0, elm.lat[0], elm.lat[1]])
                        cff = [area0/area, area1/area, area2/area]
                        for k in range(num_depth):
                            node_vals_temp = []
                            node_vals_salt = []
                            for n in range(3):
                                id1, id2, coeff1, coeff2 = elm.nodes[n].find_depth(depth_reg[k])
                                if any(isnan(array([id1, id2, coeff1, coeff2]))):
                                    node_vals_temp.append(NaN)
                                    node_vals_salt.append(NaN)
                                else:
                                    node_vals_temp.append(coeff1*fesom_temp_nodes_hr[id1] + coeff2*fesom_temp_nodes_hr[id2])
                                    node_vals_salt.append(coeff1*fesom_salt_nodes_hr[id1] + coeff2*fesom_salt_nodes_hr[id2])
                            if any(isnan(node_vals_temp)):
                                pass
                            else:
                                temp_reg_hr[k,j] = sum(array(cff)*array(node_vals_temp))
                                salt_reg_hr[k,j] = sum(array(cff)*array(node_vals_salt))
    temp_reg_hr = ma.masked_where(isnan(temp_reg_hr), temp_reg_hr)
    salt_reg_hr = ma.masked_where(isnan(salt_reg_hr), salt_reg_hr)

    depth_reg = -1*depth_reg

    # Set up axis labels the way we want them
    lat_ticks = arange(lat_min+3, lat_max+10, 10)
    lat_labels = []
    for val in lat_ticks:
        lat_labels.append(str(int(round(-val))) + r'$^{\circ}$S')
    depth_ticks = range(depth_min+1000, 0+1000, 1000)
    depth_labels = []
    for val in depth_ticks:
        depth_labels.append(str(int(round(-val))))

    print 'Plotting'
    fig = figure(figsize=(14,24))
    # ECCO2
    gs1 = GridSpec(1,2)
    gs1.update(left=0.1, right=0.95, bottom=0.7575, top=0.94, wspace=0.08)
    # Temperature
    ax = subplot(gs1[0,0])
    pcolor(ecco_lat, ecco_depth, ecco_temp, vmin=temp_min, vmax=temp_max, cmap='jet')
    # Overlay contour
    contour(ecco_lat, ecco_depth, ecco_temp, levels=[temp_contour], color='black')
    title(r'Temperature ($^{\circ}$C)', fontsize=24)
    ylabel('Depth (m)', fontsize=18)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels, fontsize=16)
    text(-64, 1000, 'a) ECCO2 initial conditions at ' + lon_string + ', January 1992', fontsize=28)
    # Salinity
    ax = subplot(gs1[0,1])
    pcolor(ecco_lat, ecco_depth, ecco_salt, vmin=salt_min, vmax=salt_max, cmap='jet')
    contour(ecco_lat, ecco_depth, ecco_salt, levels=[salt_contour], color='black')
    title('Salinity (psu)', fontsize=24)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels([])
    # MetROMS
    gs2 = GridSpec(1,2)
    gs2.update(left=0.1, right=0.95, bottom=0.525, top=0.7075, wspace=0.08)
    # Temperature
    ax = subplot(gs2[0,0])
    pcolor(roms_lat, roms_z, roms_temp, vmin=temp_min, vmax=temp_max, cmap='jet')
    contour(roms_lat, roms_z, roms_temp, levels=[temp_contour], color='black')
    ylabel('Depth (m)', fontsize=18)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels, fontsize=16)
    text(-49, 300, 'b) MetROMS, January 2016', fontsize=28)
    # Salinity
    ax = subplot(gs2[0,1])
    pcolor(roms_lat, roms_z, roms_salt, vmin=salt_min, vmax=salt_max, cmap='jet')
    contour(roms_lat, roms_z, roms_salt, levels=[salt_contour], color='black')
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels([])
    # FESOM low-res
    gs3 = GridSpec(1,2)
    gs3.update(left=0.1, right=0.95, bottom=0.2925, top=0.475, wspace=0.08)
    # Temperature
    ax = subplot(gs3[0,0])
    img = PatchCollection(patches_lr, cmap='jet')
    img.set_array(array(fesom_temp_lr))
    img.set_edgecolor('face')
    img.set_clim(vmin=temp_min, vmax=temp_max)
    ax.add_collection(img)
    # Overlay contour on regular grid
    contour(lat_reg, depth_reg, temp_reg_lr, levels=[temp_contour], color='black')
    ylabel('Depth (m)', fontsize=18)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels, fontsize=16)
    text(-53, 300, 'c) FESOM (low-res), January 2016', fontsize=28)
    # Salinity
    ax = subplot(gs3[0,1])
    img = PatchCollection(patches_lr, cmap='jet')
    img.set_array(array(fesom_salt_lr))
    img.set_edgecolor('face')
    img.set_clim(vmin=salt_min, vmax=salt_max)
    ax.add_collection(img)
    contour(lat_reg, depth_reg, salt_reg_lr, levels=[salt_contour], color='black')
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels([])
    # FESOM high-res
    gs4 = GridSpec(1,2)
    gs4.update(left=0.1, right=0.95, bottom=0.06, top=0.2425, wspace=0.08)
    # Temperature
    ax = subplot(gs4[0,0])
    img = PatchCollection(patches_hr, cmap='jet')
    img.set_array(array(fesom_temp_hr))
    img.set_edgecolor('face')
    img.set_clim(vmin=temp_min, vmax=temp_max)
    ax.add_collection(img)
    contour(lat_reg, depth_reg, temp_reg_hr, levels=[temp_contour], color='black')
    ylabel('Depth (m)', fontsize=18)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels, fontsize=16)
    text(-53, 300, 'd) FESOM (high-res), January 2016', fontsize=28)
    # Add a colorbar for temperature
    cbaxes = fig.add_axes([0.17, 0.015, 0.3, 0.015])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='both', ticks=arange(temp_min, temp_max+2, 2))
    cbar.ax.tick_params(labelsize=16)
    # Salinity
    ax = subplot(gs4[0,1])
    img = PatchCollection(patches_hr, cmap='jet')
    img.set_array(array(fesom_salt_hr))
    img.set_edgecolor('face')
    img.set_clim(vmin=salt_min, vmax=salt_max)
    ax.add_collection(img)
    contour(lat_reg, depth_reg, salt_reg_hr, levels=[salt_contour], color='black')
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels, fontsize=16)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels([])
    # Add a colorbar for salinity
    cbaxes = fig.add_axes([0.6, 0.015, 0.3, 0.02])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes, extend='both', ticks=arange(salt_min+0.1, salt_max+0.1, 0.2))
    cbar.ax.tick_params(labelsize=16)
    fig.show()
    fig.savefig('ts_drift.png')


# Command-line interface
if __name__ == "__main__":

    roms_grid = raw_input("Path to ROMS grid file: ")
    roms_file = raw_input("Path to ROMS file containing monthly averaged temperature and salinity for January 2016: ")
    fesom_mesh_path_lr = raw_input("Path to FESOM low-res mesh directory: ")
    fesom_file_lr = raw_input("Path to FESOM low-res file containing monthly averaged temperature and salinity for January 2016: ")
    fesom_mesh_path_hr = raw_input("Path to FESOM high-res mesh directory: ")
    fesom_file_hr = raw_input("Path to FESOM high-res file containing monthly averaged temperature and salinity for January 2016: ")
    mip_drift_slices(roms_grid, roms_file, fesom_mesh_path_lr, fesom_file_lr, fesom_mesh_path_hr, fesom_file_hr)
    
    
