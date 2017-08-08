from netCDF4 import Dataset
from numpy import *
from numpy.ma import MaskedArray
from matplotlib.pyplot import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.cm import *
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib import rcParams
from rotate_vector_roms import *
from cartesian_grid_3d import *
from calc_z import *
from interp_lon_roms import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from patches import *
from fesom_grid import *
from fesom_sidegrid import *
from unrotate_vector import *
from unrotate_grid import *

# This is the giant monster script to generate 8 multi-part figures showing
# ice shelf processes in the 8 regions defined in the intercomparison paper.
# Each figure is size Nx3 where N=3,4,5 is the number of variables shown
# (current options are ice shelf melt rate, ice shelf draft, bottom water
# temperature or salinity, vertically averaged velocity, zonal slices of
# temperature and salinity) and 3 is the number of models (MetROMS, low-res
# FESOM, high-res FESOM). All variables are averaged over the period 2002-2016,
# i.e. the first 10 years of the simulation are discarded as spinup, and
# seasonal cycles and interannual variability are not considered. 

# In order to avoid unnecessarily repetitive code, this script contains
# multiple functions which accept location bounds as arguments (eg to calculate
# min/max over the given region, to plot ice shelf draft over the given region
# for each model, etc.) but in order to prevent the function argument lists from
# getting out of control, most of the variables are global and so the majority
# of this script is not wrapped within a function. Be very careful if you are
# running this script in an active python session to make sure there are no 
# conflicting variable names previously defined.

# To modify for your purposes, adjust the file paths (immediately below) as
# needed, and then scroll down to "USER MODIFIED SECTION" to see examples of
# how to call the functions to create multi-part figures. You will have to
# fiddle around with the x and y bounds, GridSpec arguments, colourbar axes,
# etc. so that the positioning looks okay.

# File paths
roms_grid = '/short/m68/kaa561/metroms_iceshelf/apps/common/grid/circ30S_quarterdegree.nc'
roms_file = '/short/m68/kaa561/metroms_iceshelf/tmproms/run/intercomparison/2002_2016_avg.nc'
fesom_mesh_path_lr = '/short/y99/kaa561/FESOM/mesh/low_res/'
fesom_mesh_path_hr = '/short/y99/kaa561/FESOM/mesh/high_res/'
fesom_file_lr_o = '/short/y99/kaa561/FESOM/intercomparison_lowres/output/oce_2002_2016_avg.nc'
fesom_file_hr_o = '/short/y99/kaa561/FESOM/intercomparison_highres/output/oce_2002_2016_avg.nc'
fesom_file_lr_i = '/short/y99/kaa561/FESOM/intercomparison_lowres/output/wnet_2002_2016_avg.nc'
fesom_file_hr_i = '/short/y99/kaa561/FESOM/intercomparison_highres/output/wnet_2002_2016_avg.nc'
# Parameters for missing circle in ROMS grid
lon_c = 50
lat_c = -83
radius = 10.1
nbdry = -63+90
# Constants
deg2rad = pi/180.0
sec_per_year = 365.25*24*3600
# ROMS vertical grid parameters
theta_s = 7.0
theta_b = 2.0
hc = 250
N = 31
# Minimum zice in ROMS
min_zice = -10

print 'Reading ROMS grid'
# Read the fields we need
id = Dataset(roms_grid, 'r')
roms_lon = id.variables['lon_rho'][:,:]
roms_lat = id.variables['lat_rho'][:,:]
roms_h = id.variables['h'][:,:]
roms_mask = id.variables['mask_rho'][:,:]
roms_zice = id.variables['zice'][:,:]
roms_angle = id.variables['angle'][:,:]
id.close()
# Get land/zice mask
open_ocn = copy(roms_mask)
open_ocn[roms_zice!=0] = 0
land_zice = ma.masked_where(open_ocn==1, open_ocn)
# Convert grid to spherical coordinates
roms_x = -(roms_lat+90)*cos(roms_lon*deg2rad+pi/2)
roms_y = (roms_lat+90)*sin(roms_lon*deg2rad+pi/2)
# Find centre in spherical coordinates
x_c = -(lat_c+90)*cos(lon_c*deg2rad+pi/2)
y_c = (lat_c+90)*sin(lon_c*deg2rad+pi/2)
# Build a regular x-y grid and select the missing circle
x_reg_roms, y_reg_roms = meshgrid(linspace(-nbdry, nbdry, num=1000), linspace(-nbdry, nbdry, num=1000))
land_circle = zeros(shape(x_reg_roms))
land_circle = ma.masked_where(sqrt((x_reg_roms-x_c)**2 + (y_reg_roms-y_c)**2) > radius, land_circle)

print 'Building FESOM low-res mesh'
# Mask open ocean
elements_lr, mask_patches_lr = make_patches(fesom_mesh_path_lr, circumpolar=True, mask_cavities=True)
# Unmask ice shelves
patches_lr = iceshelf_mask(elements_lr)
# Also make a set of patches with open ocean unmasked (for bottom water T/S)
patches_all_lr = []
for elm in elements_lr:
    coord = transpose(vstack((elm.x, elm.y)))
    patches_all_lr.append(Polygon(coord, True, linewidth=0.))
# Also make non-circumpolar set of elements for zonal slices
elm2D_lr = fesom_grid(fesom_mesh_path_lr)
print 'Building FESOM high-res mesh'
elements_hr, mask_patches_hr = make_patches(fesom_mesh_path_hr, circumpolar=True, mask_cavities=True)
patches_hr = iceshelf_mask(elements_hr)
patches_all_hr = []
for elm in elements_hr:
    coord = transpose(vstack((elm.x, elm.y)))
    patches_all_hr.append(Polygon(coord, True, linewidth=0.))
elm2D_hr = fesom_grid(fesom_mesh_path_hr)

print 'Building ice shelf front contours'
# MetROMS
# Make a copy of zice field, edit so grounding line not contoured
zice_contour = ma.masked_where(roms_mask==0, roms_zice)
# FESOM low-res
contour_lines_lr = []
for elm in elements_lr:
    # Select elements where exactly 2 of the 3 nodes are in a cavity
    if count_nonzero(elm.cavity_nodes) == 2:
        # Save the coastal flags and x- and y- coordinates of these 2
        coast_tmp = []
        x_tmp = []
        y_tmp = []
        for i in range(3):
            if elm.cavity_nodes[i]:
                coast_tmp.append(elm.coast_nodes[i])
                x_tmp.append(elm.x[i])
                y_tmp.append(elm.y[i])
        # Select elements where at most 1 of these 2 nodes are coastal
        if count_nonzero(coast_tmp) < 2:
            # Draw a line between the 2 nodes
            contour_lines_lr.append([(x_tmp[0], y_tmp[0]), (x_tmp[1], y_tmp[1])])
# FESOM high-res
contour_lines_hr = []
for elm in elements_hr:
    if count_nonzero(elm.cavity_nodes) == 2:
        coast_tmp = []
        x_tmp = []
        y_tmp = []
        for i in range(3):
            if elm.cavity_nodes[i]:
                coast_tmp.append(elm.coast_nodes[i])
                x_tmp.append(elm.x[i])
                y_tmp.append(elm.y[i])
        if count_nonzero(coast_tmp) < 2:
            contour_lines_hr.append([(x_tmp[0], y_tmp[0]), (x_tmp[1], y_tmp[1])])

print 'Calculating ice shelf draft'

# ROMS
# Swap sign on existing zice field
roms_draft = -1*roms_zice
# Mask the open ocean and land
roms_draft = ma.masked_where(roms_zice==0, roms_draft)

# FESOM low-res
# Calculate draft at each element, averaged over 3 corners
# Equivalent to depth of surfce layer
fesom_draft_lr = []
for elm in elements_lr:
    if elm.cavity:
        fesom_draft_lr.append(mean([elm.nodes[0].depth, elm.nodes[1].depth, elm.nodes[2].depth]))

# FESOM high-res
fesom_draft_hr = []
for elm in elements_hr:
    if elm.cavity:
        fesom_draft_hr.append(mean([elm.nodes[0].depth, elm.nodes[1].depth, elm.nodes[2].depth]))

print 'Calculating ice shelf melt rate'

# ROMS
id = Dataset(roms_file, 'r')
# Convert from m/s to m/y
roms_melt = id.variables['m'][0,:,:]*sec_per_year
id.close()
# Mask the open ocean and land
roms_melt = ma.masked_where(roms_zice==0, roms_melt)

# FESOM low-res
# Read melt rate at 2D nodes
id = Dataset(fesom_file_lr_i, 'r')
node_melt_lr = id.variables['wnet'][0,:]*sec_per_year
id.close()
# For each element, calculate average over 3 corners
fesom_melt_lr = []
for elm in elements_lr:
    if elm.cavity:
        fesom_melt_lr.append(mean([node_melt_lr[elm.nodes[0].id], node_melt_lr[elm.nodes[1].id], node_melt_lr[elm.nodes[2].id]]))

# FESOM high-res
id = Dataset(fesom_file_hr_i, 'r')
node_melt_hr = id.variables['wnet'][0,:]*sec_per_year
id.close()
fesom_melt_hr = []
for elm in elements_hr:
    if elm.cavity:
        fesom_melt_hr.append(mean([node_melt_hr[elm.nodes[0].id], node_melt_hr[elm.nodes[1].id], node_melt_hr[elm.nodes[2].id]]))

print 'Calculating bottom water temperature'

# ROMS
id = Dataset(roms_file, 'r')
# Read bottom layer
roms_bwtemp = id.variables['temp'][0,0,:,:]
id.close()
# Mask open ocean and land
#roms_bwtemp = ma.masked_where(roms_zice==0, roms_bwtemp)

# FESOM low-res
# Read full 3D field to start
id = Dataset(fesom_file_lr_o, 'r')
node_bwtemp_lr = id.variables['temp'][0,:]
id.close()
# Calculate average over 3 corners of each bottom element
fesom_bwtemp_lr = []
for elm in elements_lr:
    #if elm.cavity:
    fesom_bwtemp_lr.append(mean([node_bwtemp_lr[elm.nodes[0].find_bottom().id], node_bwtemp_lr[elm.nodes[1].find_bottom().id], node_bwtemp_lr[elm.nodes[2].find_bottom().id]]))

# FESOM high-res
id = Dataset(fesom_file_hr_o, 'r')
node_bwtemp_hr = id.variables['temp'][0,:]
id.close()
fesom_bwtemp_hr = []
for elm in elements_hr:
    #if elm.cavity:
    fesom_bwtemp_hr.append(mean([node_bwtemp_hr[elm.nodes[0].find_bottom().id], node_bwtemp_hr[elm.nodes[1].find_bottom().id], node_bwtemp_hr[elm.nodes[2].find_bottom().id]]))    

print 'Calculating bottom water salinity'

# ROMS
id = Dataset(roms_file, 'r')
# Read bottom layer
roms_bwsalt = id.variables['salt'][0,0,:,:]
id.close()
# Mask open ocean and land
#roms_bwsalt = ma.masked_where(roms_zice==0, roms_bwsalt)

# FESOM low-res
# Read full 3D field to start
id = Dataset(fesom_file_lr_o, 'r')
node_bwsalt_lr = id.variables['salt'][0,:]
id.close()
# Calculate average over 3 corners of each bottom element
fesom_bwsalt_lr = []
for elm in elements_lr:
    #if elm.cavity:
    fesom_bwsalt_lr.append(mean([node_bwsalt_lr[elm.nodes[0].find_bottom().id], node_bwsalt_lr[elm.nodes[1].find_bottom().id], node_bwsalt_lr[elm.nodes[2].find_bottom().id]]))

# FESOM high-res
id = Dataset(fesom_file_hr_o, 'r')
node_bwsalt_hr = id.variables['salt'][0,:]
fesom_bwsalt_hr = []
id.close()
for elm in elements_hr:
    #if elm.cavity:
    fesom_bwsalt_hr.append(mean([node_bwsalt_hr[elm.nodes[0].find_bottom().id], node_bwsalt_hr[elm.nodes[1].find_bottom().id], node_bwsalt_hr[elm.nodes[2].find_bottom().id]]))    

print 'Calculating vertically averaged velocity'

# ROMS
# Read full 3D u and v
id = Dataset(roms_file, 'r')
u_3d_tmp = id.variables['u'][0,:,:,:]
v_3d_tmp = id.variables['v'][0,:,:,:]
id.close()
# Get integrands on 3D grid; we only care about dz
dx, dy, dz, z = cartesian_grid_3d(roms_lon, roms_lat, roms_h, roms_zice, theta_s, theta_b, hc, N)
# Unrotate each vertical level
u_3d = ma.empty(shape(dz))
v_3d = ma.empty(shape(dz))
num_lat_u = size(u_3d_tmp,1)
num_lon_u = size(u_3d_tmp,2)
num_lat_v = size(v_3d_tmp,1)
num_lon_v = size(v_3d_tmp,2)
for k in range(N):
    # Extend into land mask before interpolation to rho-grid so
    # the land mask doesn't change in the final plot
    for j in range(1,num_lat_u-1):
        for i in range(1,num_lon_u-1):
            # Check for masked points
            if u_3d_tmp[k,j,i] is ma.masked:
                # Look at 4 neighbours
                neighbours = ma.array([u_3d_tmp[k,j-1,i], u_3d_tmp[k,j,i-1], u_3d_tmp[k,j+1,i], u_3d_tmp[k,j,i+1]])
                # Find how many of them are unmasked
                num_unmasked = MaskedArray.count(neighbours)
                if num_unmasked > 0:
                    # There is at least one unmasked neighbour;
                    # set u_3d_tmp to their average
                    u_3d_tmp[k,j,i] = sum(neighbours)/num_unmasked
    # Repeat for v
    for j in range(1,num_lat_v-1):
        for i in range(1,num_lon_v-1):
            if v_3d_tmp[k,j,i] is ma.masked:
                neighbours = ma.array([v_3d_tmp[k,j-1,i], v_3d_tmp[k,j,i-1], v_3d_tmp[k,j+1,i], v_3d_tmp[k,j,i+1]])
                num_unmasked = MaskedArray.count(neighbours)
                if num_unmasked > 0:
                    v_3d_tmp[k,j,i] = sum(neighbours)/num_unmasked
    # Interpolate to rho grid and rotate
    u_k, v_k = rotate_vector_roms(u_3d_tmp[k,:,:], v_3d_tmp[k,:,:], roms_angle)
    u_3d[k,:,:] = u_k
    v_3d[k,:,:] = v_k
# Vertically average u and v
roms_u = sum(u_3d*dz, axis=0)/sum(dz, axis=0)
roms_v = sum(v_3d*dz, axis=0)/sum(dz, axis=0)
# Mask the open ocean and land
roms_u = ma.masked_where(roms_zice==0, roms_u)
roms_v = ma.masked_where(roms_zice==0, roms_v)
# Calculate speed
roms_speed = sqrt(roms_u**2 + roms_v**2)

# FESOM low-res
# The overlaid vectors are based on nodes not elements, so many
# of the fesom_grid data structures fail to apply and we need to
# read some of the FESOM grid files again.
# Read the cavity flag for each 2D surface node
fesom_cavity_lr = []
f = open(fesom_mesh_path_lr + 'cavity_flag_nod2d.out', 'r')
for line in f:
    tmp = int(line)
    if tmp == 1:
        fesom_cavity_lr.append(True)
    elif tmp == 0:
        fesom_cavity_lr.append(False)
    else:
        print 'Problem'
f.close()
# Save the number of 2D nodes
fesom_n2d_lr = len(fesom_cavity_lr)
# Read rotated lat and lon for each node, also depth
f = open(fesom_mesh_path_lr + 'nod3d.out', 'r')
f.readline()
rlon_lr = []
rlat_lr = []
node_depth_lr = []
for line in f:
    tmp = line.split()
    lon_tmp = float(tmp[1])
    lat_tmp = float(tmp[2])
    node_depth_tmp = -1*float(tmp[3])
    if lon_tmp < -180:
        lon_tmp += 360
    elif lon_tmp > 180:
        lon_tmp -= 360
    rlon_lr.append(lon_tmp)
    rlat_lr.append(lat_tmp)
    node_depth_lr.append(node_depth_tmp)
f.close()
# For lat and lon, only care about the 2D nodes (the first
# fesom_n2d indices)
rlon_lr = array(rlon_lr[0:fesom_n2d_lr])
rlat_lr = array(rlat_lr[0:fesom_n2d_lr])
node_depth_lr = array(node_depth_lr)
# Unrotate longitude
fesom_lon_lr, fesom_lat_lr = unrotate_grid(rlon_lr, rlat_lr)
# Calculate polar coordinates of each node
fesom_x_lr = -(fesom_lat_lr+90)*cos(fesom_lon_lr*deg2rad+pi/2)
fesom_y_lr = (fesom_lat_lr+90)*sin(fesom_lon_lr*deg2rad+pi/2)
# Read lists of which nodes are directly below which
f = open(fesom_mesh_path_lr + 'aux3d.out', 'r')
max_num_layers_lr = int(f.readline())
node_columns_lr = zeros([fesom_n2d_lr, max_num_layers_lr])
for n in range(fesom_n2d_lr):
    for k in range(max_num_layers_lr):
        node_columns_lr[n,k] = int(f.readline())
node_columns_lr = node_columns_lr.astype(int)
f.close()
# Now we can read the data
# Read full 3D field for both u and v
id = Dataset(fesom_file_lr_o, 'r')
node_ur_3d_lr = id.variables['u'][0,:]
node_vr_3d_lr = id.variables['v'][0,:]
id.close()
# Vertically average
node_ur_lr = zeros(fesom_n2d_lr)
node_vr_lr = zeros(fesom_n2d_lr)
for n in range(fesom_n2d_lr):
    # Integrate udz, vdz, and dz over this water column
    udz_col = 0
    vdz_col = 0
    dz_col = 0
    for k in range(max_num_layers_lr-1):
        if node_columns_lr[n,k+1] == -999:
            # Reached the bottom
            break
        # Trapezoidal rule
        top_id = node_columns_lr[n,k]
        bot_id = node_columns_lr[n,k+1]
        dz_tmp = node_depth_lr[bot_id-1] - node_depth_lr[top_id-1]
        udz_col += 0.5*(node_ur_3d_lr[top_id-1]+node_ur_3d_lr[bot_id-1])*dz_tmp
        vdz_col += 0.5*(node_vr_3d_lr[top_id-1]+node_vr_3d_lr[bot_id-1])*dz_tmp
        dz_col += dz_tmp
    # Convert from integrals to averages
    node_ur_lr[n] = udz_col/dz_col
    node_vr_lr[n] = vdz_col/dz_col
# Unrotate
node_u_lr, node_v_lr = unrotate_vector(rlon_lr, rlat_lr, node_ur_lr, node_vr_lr)
# Calculate speed
node_speed_lr = sqrt(node_u_lr**2 + node_v_lr**2)
# Calculate speed at each element, averaged over 3 corners
fesom_speed_lr = []
for elm in elements_lr:
    if elm.cavity:
        fesom_speed_lr.append(mean([node_speed_lr[elm.nodes[0].id], node_speed_lr[elm.nodes[1].id], node_speed_lr[elm.nodes[2].id]]))

# FESOM high-res
fesom_cavity_hr = []
f = open(fesom_mesh_path_hr + 'cavity_flag_nod2d.out', 'r')
for line in f:
    tmp = int(line)
    if tmp == 1:
        fesom_cavity_hr.append(True)
    elif tmp == 0:
        fesom_cavity_hr.append(False)
    else:
        print 'Problem'
f.close()    
fesom_n2d_hr = len(fesom_cavity_hr)
f = open(fesom_mesh_path_hr + 'nod3d.out', 'r')
f.readline()
rlon_hr = []
rlat_hr = []
node_depth_hr = []
for line in f:
    tmp = line.split()
    lon_tmp = float(tmp[1])
    lat_tmp = float(tmp[2])
    node_depth_tmp = -1*float(tmp[3])
    if lon_tmp < -180:
        lon_tmp += 360
    elif lon_tmp > 180:
        lon_tmp -= 360
    rlon_hr.append(lon_tmp)
    rlat_hr.append(lat_tmp)
    node_depth_hr.append(node_depth_tmp)
f.close()
rlon_hr = array(rlon_hr[0:fesom_n2d_hr])
rlat_hr = array(rlat_hr[0:fesom_n2d_hr])
node_depth_hr = array(node_depth_hr)
fesom_lon_hr, fesom_lat_hr = unrotate_grid(rlon_hr, rlat_hr)
fesom_x_hr = -(fesom_lat_hr+90)*cos(fesom_lon_hr*deg2rad+pi/2)
fesom_y_hr = (fesom_lat_hr+90)*sin(fesom_lon_hr*deg2rad+pi/2)
f = open(fesom_mesh_path_hr + 'aux3d.out', 'r')
max_num_layers_hr = int(f.readline())
node_columns_hr = zeros([fesom_n2d_hr, max_num_layers_hr])
for n in range(fesom_n2d_hr):
    for k in range(max_num_layers_hr):
        node_columns_hr[n,k] = int(f.readline())
node_columns_hr = node_columns_hr.astype(int)
f.close()
id = Dataset(fesom_file_hr_o, 'r')
node_ur_3d_hr = id.variables['u'][0,:]
node_vr_3d_hr = id.variables['v'][0,:]
id.close()
node_ur_hr = zeros(fesom_n2d_hr)
node_vr_hr = zeros(fesom_n2d_hr)
for n in range(fesom_n2d_hr):
    udz_col = 0
    vdz_col = 0
    dz_col = 0
    for k in range(max_num_layers_hr-1):
        if node_columns_hr[n,k+1] == -999:
            break
        top_id = node_columns_hr[n,k]
        bot_id = node_columns_hr[n,k+1]
        dz_tmp = node_depth_hr[bot_id-1] - node_depth_hr[top_id-1]
        udz_col += 0.5*(node_ur_3d_hr[top_id-1]+node_ur_3d_hr[bot_id-1])*dz_tmp
        vdz_col += 0.5*(node_vr_3d_hr[top_id-1]+node_vr_3d_hr[bot_id-1])*dz_tmp
        dz_col += dz_tmp
    node_ur_hr[n] = udz_col/dz_col
    node_vr_hr[n] = vdz_col/dz_col
node_u_hr, node_v_hr = unrotate_vector(rlon_hr, rlat_hr, node_ur_hr, node_vr_hr)
node_speed_hr = sqrt(node_u_hr**2 + node_v_hr**2)
fesom_speed_hr = []
for elm in elements_hr:
    if elm.cavity:
        fesom_speed_hr.append(mean([node_speed_hr[elm.nodes[0].id], node_speed_hr[elm.nodes[1].id], node_speed_hr[elm.nodes[2].id]])) 


# **************** USER MODIFIED SECTION ****************
# Filchner-Ronne
x_min_tmp = -14
x_max_tmp = -4.5
y_min_tmp = 1
y_max_tmp = 10
fig = figure(figsize=(8,14))
fig.patch.set_facecolor('white')
# Melt rate
gs_a = GridSpec(1,3)
gs_a.update(left=0.05, right=0.9, bottom=0.735, top=0.89, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.755, 0.025, 0.12])
cbar_ticks = arange(0, 6+3, 3)
plot_melt(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_a, cbaxes_tmp, cbar_ticks, [0.5, 3, 4.5], 'a', 1.25, [-60, -40], [-80])
# Velocity
x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr, fesom_ubin_hr, fesom_vbin_hr = make_vectors(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, 20, 20)
gs_b = GridSpec(1,3)
gs_b.update(left=0.05, right=0.9, bottom=0.555, top=0.71, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.575, 0.025, 0.12])
cbar_ticks = arange(0, 0.2+0.1, 0.1)
plot_velavg(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_b, cbaxes_tmp, cbar_ticks, x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr, fesom_ubin_hr, fesom_vbin_hr, 'b', arrow_scale=0.9, arrow_headwidth=8, arrow_headlength=9)
# Bottom water temperature
gs_c = GridSpec(1,3)
gs_c.update(left=0.05, right=0.9, bottom=0.375, top=0.53, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.395, 0.025, 0.12])
cbar_ticks = arange(-2.6, -1.8+0.4, 0.4)
plot_bwtemp(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_c, cbaxes_tmp, cbar_ticks, 'c')
# Bottom water salinity
gs_d = GridSpec(1,3)
gs_d.update(left=0.05, right=0.9, bottom=0.195, top=0.35, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.215, 0.025, 0.12])
cbar_ticks = arange(34.3, 34.7+0.2, 0.2)
plot_bwsalt(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_d, cbaxes_tmp, cbar_ticks, 'd')
# Ice shelf draft
gs_e = GridSpec(1,3)
gs_e.update(left=0.05, right=0.9, bottom=0.015, top=0.17, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.035, 0.025, 0.12])
cbar_ticks = arange(500, 1500+500, 500)
plot_draft(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_e, cbaxes_tmp, cbar_ticks, 'e')
suptitle('Filchner-Ronne Ice Shelf', fontsize=30)
fig.show()
fig.savefig('filchner_ronne.png')

# Eastern Weddell
x_min_tmp = -8
x_max_tmp = 13
y_min_tmp = 12
y_max_tmp = 21
fig = figure(figsize=(10,10))
fig.patch.set_facecolor('white')
# Melt rate
gs_a = GridSpec(1,3)
gs_a.update(left=0.11, right=0.9, bottom=0.74, top=0.89, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.77, 0.025, 0.1])
cbar_ticks = arange(0, 4+2, 2)
plot_melt(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_a, cbaxes_tmp, cbar_ticks, [1, 2, 3], 'a', 1.35, [0, 30], [-70])
# Velocity
x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr, fesom_ubin_hr, fesom_vbin_hr = make_vectors(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, 40, 20)
gs_b = GridSpec(1,3)
gs_b.update(left=0.11, right=0.9, bottom=0.57, top=0.72, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.6, 0.025, 0.1])
cbar_ticks = arange(0, 0.18+0.09, 0.09)
plot_velavg(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_b, cbaxes_tmp, cbar_ticks, x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr, fesom_ubin_hr, fesom_vbin_hr, 'b', arrow_scale=0.9, arrow_headwidth=8, arrow_headlength=9)
# Temp and salt slices through 1W: Fimbul Ice Shelf
lat_min = -71.5
lat_max = -69.5
lat_ticks = [-71.5, -71, -70.5, -70, -69.5]
lat_labels = ['', r'71$^{\circ}$S', '', r'70$^{\circ}$S', '']
depth_min = -2400
depth_max = 0
depth_ticks = [-2000, -1500, -1000, -500, 0]
depth_labels = ['2000', '1500', '1000', '500', '0']
gs_c = GridSpec(1,3)
gs_c.update(left=0.11, right=0.9, bottom=0.33, top=0.53, wspace=0.05)
cbaxes_tmp1 = fig.add_axes([0.91, 0.36, 0.025, 0.14])
cbar_ticks1 = arange(-1.8, 0.6+1.2, 1.2)
gs_d = GridSpec(1,3)
gs_d.update(left=0.11, right=0.9, bottom=0.05, top=0.25, wspace=0.05)
cbaxes_tmp2 = fig.add_axes([0.91, 0.08, 0.025, 0.14])
cbar_ticks2 = arange(34, 34.6+0.3, 0.3)
plot_zonal_ts(-1, lat_min, lat_max, lat_ticks, lat_labels, depth_min, depth_max, depth_ticks, depth_labels, gs_c, gs_d, cbaxes_tmp1, cbaxes_tmp2, cbar_ticks1, cbar_ticks2, 'c', 'd', loc_string='Fimbul Ice Shelf')
suptitle('Eastern Weddell Region', fontsize=30)
fig.show()
fig.savefig('eweddell.png')

# Amery
x_min_tmp = 15.25
x_max_tmp = 20.5
y_min_tmp = 4.75
y_max_tmp = 8
fig = figure(figsize=(10,12))
fig.patch.set_facecolor('white')
# Melt
gs_a = GridSpec(1,3)
gs_a.update(left=0.05, right=0.9, bottom=0.735, top=0.89, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.755, 0.025, 0.12])
cbar_ticks = arange(0, 40+20, 20)
plot_melt(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_a, cbaxes_tmp, cbar_ticks, [2, 5, 11], 'a', 1.25, [70], [-70])
# Velocity
x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr, fesom_ubin_hr, fesom_vbin_hr = make_vectors(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, 20, 15)
gs_b = GridSpec(1,3)
gs_b.update(left=0.05, right=0.9, bottom=0.555, top=0.71, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.575, 0.025, 0.12])
cbar_ticks = arange(0.1, 0.3+0.1, 0.1)
plot_velavg(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_b, cbaxes_tmp, cbar_ticks, x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr, fesom_ubin_hr, fesom_vbin_hr, 'b', arrow_scale=0.9, arrow_headwidth=7, arrow_headlength=8)
# Bottom water temperature
gs_c = GridSpec(1,3)
gs_c.update(left=0.05, right=0.9, bottom=0.375, top=0.53, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.395, 0.025, 0.12])
cbar_ticks = arange(-2.4, -1.4+0.5, 0.5)
plot_bwtemp(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_c, cbaxes_tmp, cbar_ticks, 'c')
# Bottom water salinity
gs_d = GridSpec(1,3)
gs_d.update(left=0.05, right=0.9, bottom=0.195, top=0.35, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.215, 0.025, 0.12])
cbar_ticks = arange(34, 34.4+0.2, 0.2)
plot_bwsalt(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_d, cbaxes_tmp, cbar_ticks, 'd')
# Ice shelf draft
gs_e = GridSpec(1,3)
gs_e.update(left=0.05, right=0.9, bottom=0.015, top=0.17, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.035, 0.025, 0.12])
cbar_ticks = arange(200, 2200+1000, 1000)
plot_draft(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_e, cbaxes_tmp, cbar_ticks, 'e')
suptitle('Amery Ice Shelf', fontsize=30)
fig.show()
fig.savefig('amery.png')

# Australian sector
x_min_tmp = 12
x_max_tmp = 25.5
y_min_tmp = -20
y_max_tmp = 4
fig = figure(figsize=(8,14))
fig.patch.set_facecolor('white')
# Melt
gs_a = GridSpec(1,3)
gs_a.update(left=0.05, right=0.9, bottom=0.68, top=0.9, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.72, 0.025, 0.15])
cbar_ticks = arange(0, 4+2, 2)
plot_melt(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_a, cbaxes_tmp, cbar_ticks, [1, 2, 3], 'a', 1.15, [90, 120], [-65])
# Bottom water temperature
gs_b = GridSpec(1,3)
gs_b.update(left=0.05, right=0.9, bottom=0.435, top=0.655, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.475, 0.025, 0.15])
cbar_ticks = arange(-1.8, -0.8+0.5, 0.5)
plot_bwtemp(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_b, cbaxes_tmp, cbar_ticks, 'b')
# Bottom water salinity
gs_c = GridSpec(1,3)
gs_c.update(left=0.05, right=0.9, bottom=0.19, top=0.41, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.23, 0.025, 0.15])
cbar_ticks = arange(34.2, 34.5+0.1, 0.1)
plot_bwsalt(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_c, cbaxes_tmp, cbar_ticks, 'c')
# Velocity for Totten
x_min_tmp = 20
x_max_tmp = 21.5
y_min_tmp = -10.9
y_max_tmp = -9.7
x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr, fesom_ubin_hr, fesom_vbin_hr = make_vectors(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, 18, 18)
gs_d = GridSpec(1,3)
gs_d.update(left=0.05, right=0.9, bottom=0.025, top=0.165, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.05, 0.025, 0.08])
cbar_ticks = arange(0.02, 0.08+0.03, 0.03)
plot_velavg(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_d, cbaxes_tmp, cbar_ticks, x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr, fesom_ubin_hr, fesom_vbin_hr, 'd', loc_string='Totten Ice Shelf', arrow_scale=0.7, arrow_headwidth=8, arrow_headlength=9)
suptitle('Australian sector', fontsize=30)
fig.show()
fig.savefig('australian.png')

# Ross
x_min_tmp = -9.5
x_max_tmp = 4
y_min_tmp = -13
y_max_tmp = -4.75
fig = figure(figsize=(10,12))
fig.patch.set_facecolor('white')
# Melt
gs_a = GridSpec(1,3)
gs_a.update(left=0.1, right=0.9, bottom=0.735, top=0.89, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.765, 0.025, 0.1])
cbar_ticks = arange(0, 8+4, 4)
plot_melt(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_a, cbaxes_tmp, cbar_ticks, [0.5, 2, 4], 'a', 1.25, [180, -160], [-80])
# Velocity
x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr, fesom_ubin_hr, fesom_vbin_hr = make_vectors(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, 20, 15)
gs_b = GridSpec(1,3)
gs_b.update(left=0.1, right=0.9, bottom=0.56, top=0.715, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.59, 0.025, 0.1])
cbar_ticks = arange(0, 0.18+0.09, 0.09)
plot_velavg(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_b, cbaxes_tmp, cbar_ticks, x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr, fesom_ubin_hr, fesom_vbin_hr, 'b', arrow_scale=0.5, arrow_headwidth=8, arrow_headlength=9)
# Draft
gs_c = GridSpec(1,3)
gs_c.update(left=0.1, right=0.9, bottom=0.385, top=0.54, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.415, 0.025, 0.1])
cbar_ticks = arange(200, 600+200, 200)
plot_draft(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_c, cbaxes_tmp, cbar_ticks, 'c')
# Temp and salt slices through 180E
lat_min = -84.5
lat_max = -77
lat_ticks = [-84, -82, -80, -78]
lat_labels = [r'84$^{\circ}$S', r'82$^{\circ}$S', r'80$^{\circ}$S', r'78$^{\circ}$S']
depth_min = -750
depth_max = 0
depth_ticks = [-750, -500, -250, 0]
depth_labels = ['750', '500', '250', '0']
gs_d = GridSpec(1,3)
gs_d.update(left=0.1, right=0.9, bottom=0.23, top=0.36, wspace=0.05)
cbaxes_tmp1 = fig.add_axes([0.91, 0.245, 0.025, 0.1])
cbar_ticks1 = arange(-2.1, -1.5+0.3, 0.3)
gs_e = GridSpec(1,3)
gs_e.update(left=0.1, right=0.9, bottom=0.05, top=0.18, wspace=0.05)
cbaxes_tmp2 = fig.add_axes([0.91, 0.065, 0.025, 0.1])
cbar_ticks2 = arange(34.4, 35+0.3, 0.3)
plot_zonal_ts(180, lat_min, lat_max, lat_ticks, lat_labels, depth_min, depth_max, depth_ticks, depth_labels, gs_d, gs_e, cbaxes_tmp1, cbaxes_tmp2, cbar_ticks1, cbar_ticks2, 'c', 'd')
suptitle('Ross Sea', fontsize=30)
fig.show()
fig.savefig('ross.png')

# Amundsen Sea
x_min_tmp = -15.5
x_max_tmp = -10.5
y_min_tmp = -11.25
y_max_tmp = -2.25
fig = figure(figsize=(8,14))
fig.patch.set_facecolor('white')
# Melt
gs_a = GridSpec(1,3)
gs_a.update(left=0.05, right=0.9, bottom=0.68, top=0.9, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.72, 0.025, 0.15])
cbar_ticks = arange(0, 14+7, 7)
plot_melt(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_a, cbaxes_tmp, cbar_ticks, [1, 3, 6], 'a', 1.15, [-130, -110], [-75])
# Bottom water temperature
gs_b = GridSpec(1,3)
gs_b.update(left=0.05, right=0.9, bottom=0.435, top=0.655, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.475, 0.025, 0.15])
cbar_ticks = arange(-1.6, 0+0.8, 0.8)
plot_bwtemp(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_b, cbaxes_tmp, cbar_ticks, 'b')
# Bottom water salinity
gs_c = GridSpec(1,3)
gs_c.update(left=0.05, right=0.9, bottom=0.19, top=0.41, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.23, 0.025, 0.15])
cbar_ticks = arange(34.1, 34.5+0.2, 0.2)
plot_bwsalt(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_c, cbaxes_tmp, cbar_ticks, 'c')
# Velocity for PIG
x_min_tmp = -15.6
x_max_tmp = -14.1
y_min_tmp = -3.5
y_max_tmp = -2.25
x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr, fesom_ubin_hr, fesom_vbin_hr = make_vectors(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, 18, 18)
gs_d = GridSpec(1,3)
gs_d.update(left=0.05, right=0.9, bottom=0.025, top=0.165, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.05, 0.025, 0.08])
cbar_ticks = arange(0.03, 0.09+0.03, 0.03)
plot_velavg(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_d, cbaxes_tmp, cbar_ticks, x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr, fesom_ubin_hr, fesom_vbin_hr, 'd', loc_string='Pine Island Glacier Ice Shelf', arrow_scale=0.5, arrow_headwidth=8, arrow_headlength=9)
suptitle('Amundsen Sea', fontsize=30)
fig.show()
fig.savefig('amundsen.png')

# Bellingshausen Sea
x_min_tmp = -20.25
x_max_tmp = -15.5
y_min_tmp = -4.5
y_max_tmp = 7.6
fig = figure(figsize=(8,14))
fig.patch.set_facecolor('white')
# Melt
gs_a = GridSpec(1,3)
gs_a.update(left=0.05, right=0.9, bottom=0.68, top=0.9, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.72, 0.025, 0.15])
cbar_ticks = arange(0, 12+6, 6)
plot_melt(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_a, cbaxes_tmp, cbar_ticks, [0.5, 2, 4], 'a', 1.15, [-100, -80], [-70])
# Bottom water temperature
gs_b = GridSpec(1,3)
gs_b.update(left=0.05, right=0.9, bottom=0.435, top=0.655, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.475, 0.025, 0.15])
cbar_ticks = arange(-2, 0+1, 1)
plot_bwtemp(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_b, cbaxes_tmp, cbar_ticks, 'b')
# Bottom water salinity
gs_c = GridSpec(1,3)
gs_c.update(left=0.05, right=0.9, bottom=0.19, top=0.41, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.23, 0.025, 0.15])
cbar_ticks = arange(33.6, 34.4+0.4, 0.4)
plot_bwsalt(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_c, cbaxes_tmp, cbar_ticks, 'c')
# Velocity for George VI
x_min_tmp = -18.75
x_max_tmp = -15.5
y_min_tmp = 4.25
y_max_tmp = 7.6
x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr, fesom_ubin_hr, fesom_vbin_hr = make_vectors(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, 18, 18)
gs_d = GridSpec(1,3)
gs_d.update(left=0.05, right=0.9, bottom=0.025, top=0.165, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.05, 0.025, 0.08])
cbar_ticks = arange(0, 0.08+0.04, 0.04)
plot_velavg(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_d, cbaxes_tmp, cbar_ticks, x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr, fesom_ubin_hr, fesom_vbin_hr, 'd', loc_string='George VI Ice Shelf', arrow_scale=0.4, arrow_headwidth=8, arrow_headlength=9)
suptitle('Bellingshausen Sea', fontsize=30)
fig.show()
fig.savefig('bellingshausen.png')

# Larsen
x_min_tmp = -22.5
x_max_tmp = -14.5
y_min_tmp = 8.3
y_max_tmp = 13
fig = figure(figsize=(10,8))
fig.patch.set_facecolor('white')
# Melt
gs_a = GridSpec(1,3)
gs_a.update(left=0.05, right=0.9, bottom=0.61, top=0.84, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.7, 0.025, 0.15])
cbar_ticks = arange(0, 8+4, 4)
plot_melt(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_a, cbaxes_tmp, cbar_ticks, [0.5, 2, 4], 'a', 1.3, [-60], [-70])
# Velocity
x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr, fesom_ubin_hr, fesom_vbin_hr = make_vectors(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, 18, 18)
gs_b = GridSpec(1,3)
gs_b.update(left=0.05, right=0.9, bottom=0.33, top=0.56, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.4, 0.025, 0.15])
cbar_ticks = arange(0, 0.12+0.06, 0.06)
plot_velavg(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_b, cbaxes_tmp, cbar_ticks, x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr, fesom_ubin_hr, fesom_vbin_hr, 'b', loc_string='George VI Ice Shelf', arrow_scale=0.7, arrow_headwidth=8, arrow_headlength=9)
# Bottom water temperature
gs_c = GridSpec(1,3)
gs_c.update(left=0.05, right=0.9, bottom=0.05, top=0.28, wspace=0.05)
cbaxes_tmp = fig.add_axes([0.91, 0.1, 0.025, 0.15])
cbar_ticks = arange(-1.8, -1+0.4, 0.4)
plot_bwtemp(x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp, gs_c, cbaxes_tmp, cbar_ticks, 'c')
suptitle('Larsen Ice Shelves', fontsize=30)
fig.show()
fig.savefig('larsen.png')


# For the given variable in MetROMS, low-res FESOM, and high-res FESOM, find
# the max and min values across all 3 models in the given region.
# Input:
# roms_data = 2D field of data on the MetROMS grid (lat x lon)
# fesom_data_lr, fesom_data_hr = data for each 2D element on the FESOM low-res
#                and high-res meshes respectively
# x_min, x_max, y_min, y_max = bounds on x and y (using polar coordinate
#               transformation from above) for the desired region
# cavity = optional boolean indicating to only consider values in ice shelf
#          cavities (default True)
# Output:
# var_min, var_max = min and max data values in this region across all 3 models
def get_min_max (roms_data, fesom_data_lr, fesom_data_hr, x_min, x_max, y_min, y_max, cavity=True):

    # Start with ROMS
    loc = (roms_x >= x_min)*(roms_x <= x_max)*(roms_y >= y_min)*(roms_y <= y_max)
    var_min = amin(roms_data[loc])
    var_max = amax(roms_data[loc])
    # Modify with FESOM
    # Low-res
    i = 0
    for elm in elements_lr:
        if (not cavity) or (cavity and elm.cavity):
            if any(elm.x >= x_min) and any(elm.x <= x_max) and any(elm.y >= y_min) and any(elm.y <= y_max):
                if fesom_data_lr[i] < var_min:
                    var_min = fesom_data_lr[i]
                if fesom_data_lr[i] > var_max:
                    var_max = fesom_data_lr[i]
            i += 1
    # High-res
    i = 0
    for elm in elements_hr:
        if (not cavity) or (cavity and elm.cavity):
            if any(elm.x >= x_min) and any(elm.x <= x_max) and any(elm.y >= y_min) and any(elm.y <= y_max):
                if fesom_data_hr[i] < var_min:
                    var_min = fesom_data_hr[i]
                if fesom_data_hr[i] > var_max:
                    var_max = fesom_data_hr[i]
            i += 1
    return var_min, var_max


# Create a 2D vector field for vertically averaged velocity in MetROMS, low-res
# FESOM, and high-res FESOM for the given region. Average velocity over
# horizontal bins (in x and y) for easy plotting.
# Input:
# x_min, x_max, y_min, y_max = bounds on x and y (using polar coordinate
#               transformation from above) for the desired region
# num_bins_x, num_bins_y = number of bins to use for x and y dimensions
# Output:
# x_centres, y_centres = 1D arrays of length num_bins containing x and y
#                        coordinates of the bin centres
# roms_ubin, roms_vbin = 2D arrays of size num_bins x num_bins containing
#                        vector components for MetROMS at each bin
# fesom_ubin_lr, fesom_vbin_lr = same for FESOM low-res
# fesom_ubin_hr, fesom_vbin_hr = same for FESOM high-res
def make_vectors (x_min, x_max, y_min, y_max, num_bins_x, num_bins_y):

    # Set up bins (edges)
    x_bins = linspace(x_min, x_max, num=num_bins_x+1)
    y_bins = linspace(y_min, y_max, num=num_bins_y+1)
    # Calculate centres of bins (for plotting)
    x_centres = 0.5*(x_bins[:-1] + x_bins[1:])
    y_centres = 0.5*(y_bins[:-1] + y_bins[1:])
    # ROMS
    # First set up arrays to integrate velocity in each bin
    # Simple averaging of all the points inside each bin
    roms_ubin = zeros([size(y_centres), size(x_centres)])
    roms_vbin = zeros([size(y_centres), size(x_centres)])
    roms_num_pts = zeros([size(y_centres), size(x_centres)])
    # First convert to polar coordinates, rotate to account for
    # longitude in circumpolar projection, and convert back to vector
    # components
    theta_roms = arctan2(roms_v, roms_u)
    theta_circ_roms = theta_roms - roms_lon*deg2rad
    u_circ_roms = roms_speed*cos(theta_circ_roms)
    v_circ_roms = roms_speed*sin(theta_circ_roms)
    # Loop over all points (can't find a better way to do this)
    for j in range(size(roms_speed,0)):
        for i in range(size(roms_speed,1)):
            # Make sure data isn't masked (i.e. land or open ocean)
            if u_circ_roms[j,i] is not ma.masked:
                # Check if we're in the region of interest
                if roms_x[j,i] > x_min and roms_x[j,i] < x_max and roms_y[j,i] > y_min and roms_y[j,i] < y_max:
                    # Figure out which bins this falls into
                    x_index = nonzero(x_bins > roms_x[j,i])[0][0]-1
                    y_index = nonzero(y_bins > roms_y[j,i])[0][0]-1
                    # Integrate
                    roms_ubin[y_index, x_index] += u_circ_roms[j,i]
                    roms_vbin[y_index, x_index] += v_circ_roms[j,i]
                    roms_num_pts[y_index, x_index] += 1
    # Convert from sums to averages
    # First mask out points with no data
    roms_ubin = ma.masked_where(roms_num_pts==0, roms_ubin)
    roms_vbin = ma.masked_where(roms_num_pts==0, roms_vbin)
    # Divide everything else by the number of points
    flag = roms_num_pts > 0
    roms_ubin[flag] = roms_ubin[flag]/roms_num_pts[flag]
    roms_vbin[flag] = roms_vbin[flag]/roms_num_pts[flag]
    # FESOM low-res
    fesom_ubin_lr = zeros([size(y_centres), size(x_centres)])
    fesom_vbin_lr = zeros([size(y_centres), size(x_centres)])
    fesom_num_pts_lr = zeros([size(y_centres), size(x_centres)])
    theta_fesom_lr = arctan2(node_v_lr, node_u_lr)
    theta_circ_fesom_lr = theta_fesom_lr - fesom_lon_lr*deg2rad
    u_circ_fesom_lr = node_speed_lr*cos(theta_circ_fesom_lr)
    v_circ_fesom_lr = node_speed_lr*sin(theta_circ_fesom_lr)
    # Loop over 2D nodes to fill in the velocity bins as before
    for n in range(fesom_n2d_lr):
        if fesom_cavity_lr[n]:
            if fesom_x_lr[n] > x_min and fesom_x_lr[n] < x_max and fesom_y_lr[n] > y_min and fesom_y_lr[n] < y_max:
                x_index = nonzero(x_bins > fesom_x_lr[n])[0][0]-1
                y_index = nonzero(y_bins > fesom_y_lr[n])[0][0]-1
                fesom_ubin_lr[y_index, x_index] += u_circ_fesom_lr[n]
                fesom_vbin_lr[y_index, x_index] += v_circ_fesom_lr[n]
                fesom_num_pts_lr[y_index, x_index] += 1
    fesom_ubin_lr = ma.masked_where(fesom_num_pts_lr==0, fesom_ubin_lr)
    fesom_vbin_lr = ma.masked_where(fesom_num_pts_lr==0, fesom_vbin_lr)
    flag = fesom_num_pts_lr > 0
    fesom_ubin_lr[flag] = fesom_ubin_lr[flag]/fesom_num_pts_lr[flag]
    fesom_vbin_lr[flag] = fesom_vbin_lr[flag]/fesom_num_pts_lr[flag]
    # FESOM high-res
    fesom_ubin_hr = zeros([size(y_centres), size(x_centres)])
    fesom_vbin_hr = zeros([size(y_centres), size(x_centres)])
    fesom_num_pts_hr = zeros([size(y_centres), size(x_centres)])
    theta_fesom_hr = arctan2(node_v_hr, node_u_hr)
    theta_circ_fesom_hr = theta_fesom_hr - fesom_lon_hr*deg2rad
    u_circ_fesom_hr = node_speed_hr*cos(theta_circ_fesom_hr)
    v_circ_fesom_hr = node_speed_hr*sin(theta_circ_fesom_hr)
    for n in range(fesom_n2d_hr):
        if fesom_cavity_hr[n]:
            if fesom_x_hr[n] > x_min and fesom_x_hr[n] < x_max and fesom_y_hr[n] > y_min and fesom_y_hr[n] < y_max:
                x_index = nonzero(x_bins > fesom_x_hr[n])[0][0]-1
                y_index = nonzero(y_bins > fesom_y_hr[n])[0][0]-1
                fesom_ubin_hr[y_index, x_index] += u_circ_fesom_hr[n]
                fesom_vbin_hr[y_index, x_index] += v_circ_fesom_hr[n]
                fesom_num_pts_hr[y_index, x_index] += 1
    fesom_ubin_hr = ma.masked_where(fesom_num_pts_hr==0, fesom_ubin_hr)
    fesom_vbin_hr = ma.masked_where(fesom_num_pts_hr==0, fesom_vbin_hr)
    flag = fesom_num_pts_hr > 0
    fesom_ubin_hr[flag] = fesom_ubin_hr[flag]/fesom_num_pts_hr[flag]
    fesom_vbin_hr[flag] = fesom_vbin_hr[flag]/fesom_num_pts_hr[flag]

    return x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr, fesom_ubin_hr, fesom_vbin_hr


# Plot ice shelf draft for MetROMS, low-res FESOM, and high-res FESOM for the
# given region, in the given axes.
# Input:
# x_min, x_max, y_min, y_max = bounds on x and y (using polar coordinate
#               transformation from above) for the desired region
# gs = GridSpec object of size 1x3 to plot in
# cbaxes = Axes object for location of colourbar
# cbar_ticks = 1D array containing values for ticks on colourbar
# letter = 'a', 'b', 'c', etc. to add before the ice shelf draft title, for
#          use in a figure showing multiple variables
def plot_draft (x_min, x_max, y_min, y_max, gs, cbaxes, cbar_ticks, letter):

    # Set up a grey square for FESOM to fill the background with land
    x_reg_fesom, y_reg_fesom = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
    land_square = zeros(shape(x_reg_fesom))
    # Find bounds on variable in this region
    var_min, var_max = get_min_max(roms_draft, fesom_draft_lr, fesom_draft_hr, x_min, x_max, y_min, y_max)

    # MetROMS
    ax = subplot(gs[0,0], aspect='equal')
    # First shade land and zice in grey
    grey_cmap = ListedColormap([(0.6, 0.6, 0.6)])
    pcolor(roms_x, roms_y, land_zice, cmap=grey_cmap)
    # Fill in the missing circle
    pcolor(x_reg_roms, y_reg_roms, land_circle, cmap=grey_cmap)
    # Now shade the data
    pcolor(roms_x, roms_y, roms_draft, vmin=var_min, vmax=var_max, cmap='jet')
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    axis('off')

    # FESOM low-res
    ax = subplot(gs[0,1], aspect='equal')
    # Start with land background
    contourf(x_reg_fesom, y_reg_fesom, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    # Add ice shelf elements
    img = PatchCollection(patches_lr, cmap='jet')
    img.set_array(array(fesom_draft_lr))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    # Mask out the open ocean in white
    overlay = PatchCollection(mask_patches_lr, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    axis('off')
    # Main title
    title(letter + ') Ice shelf draft (m)', fontsize=20)

    # FESOM high-res
    ax = subplot(gs[0,2], aspect='equal')
    contourf(x_reg_fesom, y_reg_fesom, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches_hr, cmap='jet')
    img.set_array(array(fesom_draft_hr))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    overlay = PatchCollection(mask_patches_hr, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    axis('off')
    # Colourbar on the right
    cbar = colorbar(img, cax=cbaxes, ticks=cbar_ticks)


# Plot ice shelf melt rate for MetROMS, low-res FESOM, and high-res FESOM for
# the given region, in the given axes.
# Input:
# x_min, x_max, y_min, y_max = bounds on x and y (using polar coordinate
#               transformation from above) for the desired region
# gs = GridSpec object of size 1x3 to plot in
# cbaxes = Axes object for location of colourbar
# cbar_ticks = 1D array containing values for ticks on colourbar
# change_points = list of size 3 containing values where the colourmap should
#                 transition (1) from yellow to orange, (2) from orange to red,
#                 (3) from red to magenta. Should not include the minimum value,
#                 0, or the maximum value. This way the custom colourmap can be 
#                 adjusted so that all melt rates are visible, particularly
#                 for ice shelves with strong spatial variations in melt.
# letter = 'a', 'b', 'c', etc. to add before the ice shelf melt rate title, for
#          use in a figure showing multiple variables
# y0 = y-coordinate of model titles for the entire plot, assuming melt rate is
#      always at the top (i.e. letter='a'). Play around between 1.15 and 1.35.
# lon_lines = list of longitudes to overlay as dotted lines on the first panel
#             (-180 to 180)
# lat_lines = list of latitudes to overlay (-90 to 90)
def plot_melt (x_min, x_max, y_min, y_max, gs, cbaxes, cbar_ticks, change_points, letter, y0, lon_lines, lat_lines):

    # Set up a grey square for FESOM to fill the background with land
    x_reg_fesom, y_reg_fesom = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
    land_square = zeros(shape(x_reg_fesom))
    # Find bounds on variable in this region
    var_min, var_max = get_min_max(roms_melt, fesom_melt_lr, fesom_melt_hr, x_min, x_max, y_min, y_max)
    # Special colour map
    if var_min < 0:
        # There is refreezing here; include blue for elements < 0
        cmap_vals = array([var_min, 0, change_points[0], change_points[1], change_points[2], var_max])
        cmap_colors = [(0.26, 0.45, 0.86), (1, 1, 1), (1, 0.9, 0.4), (0.99, 0.59, 0.18), (0.5, 0.0, 0.08), (0.96, 0.17, 0.89)]
        cmap_vals_norm = (cmap_vals - var_min)/(var_max - var_min)
        cmap_vals_norm[-1] = 1
        cmap_list = []
        for i in range(size(cmap_vals)):
            cmap_list.append((cmap_vals_norm[i], cmap_colors[i]))
        mf_cmap = LinearSegmentedColormap.from_list('melt_freeze', cmap_list)
    else:
        # No refreezing
        cmap_vals = array([0, change_points[0], change_points[1], change_points[2], var_max])
        cmap_colors = [(1, 1, 1), (1, 0.9, 0.4), (0.99, 0.59, 0.18), (0.5, 0.0, 0.08), (0.96, 0.17, 0.89)]
        cmap_vals_norm = cmap_vals/var_max
        cmap_vals_norm[-1] = 1
        cmap_list = []
        for i in range(size(cmap_vals)):
            cmap_list.append((cmap_vals_norm[i], cmap_colors[i]))
        mf_cmap = LinearSegmentedColormap.from_list('melt_freeze', cmap_list)

    # Make sure longitudes to overlay are between 0 and 360, to suit ROMS
    # convention
    for i in range(len(lon_lines)):
        if lon_lines[i] < 0:
            lon_lines[i] += 360

    # Plot MetROMS
    ax = subplot(gs[0,0], aspect='equal')
    # First shade land and zice in grey
    grey_cmap = ListedColormap([(0.6, 0.6, 0.6)])
    pcolor(roms_x, roms_y, land_zice, cmap=grey_cmap)
    # Fill in the missing circle
    pcolor(x_reg_roms, y_reg_roms, land_circle, cmap=grey_cmap)
    # Now shade the data
    pcolor(roms_x, roms_y, roms_melt, vmin=var_min, vmax=var_max, cmap=mf_cmap)
    # Overlay longitudes
    contour(roms_x, roms_y, roms_lon, lon_lines, colors='black', linestyles='dashed')
    # Overlay latitudes
    contour(roms_x, roms_y, roms_lat, lat_lines, colors='black', linestyles='dashed')
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    axis('off')
    # Melt rate is always at the top, so add model labels
    text(0.5, y0, 'MetROMS', fontsize=18, horizontalalignment='center', transform=ax.transAxes) 

    # FESOM low-res
    ax = subplot(gs[0,1], aspect='equal')
    # Start with land background
    contourf(x_reg_fesom, y_reg_fesom, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    # Add ice shelf elements
    img = PatchCollection(patches_lr, cmap=mf_cmap)
    img.set_array(array(fesom_melt_lr))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    # Mask out the open ocean in white
    overlay = PatchCollection(mask_patches_lr, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    axis('off')
    title(letter + ') Ice shelf melt rate (m/y)', fontsize=20)
    text(0.5, y0, 'FESOM (low-res)', fontsize=18, horizontalalignment='center', transform=ax.transAxes) 

    # FESOM high-res
    ax = subplot(gs[0,2], aspect='equal')
    contourf(x_reg_fesom, y_reg_fesom, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches_hr, cmap=mf_cmap)
    img.set_array(array(fesom_melt_hr))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    overlay = PatchCollection(mask_patches_hr, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    axis('off')
    # Colourbar on the right
    cbar = colorbar(img, cax=cbaxes, ticks=cbar_ticks)
    text(0.5, y0, 'FESOM (high-res)', fontsize=18, horizontalalignment='center', transform=ax.transAxes) 


# Plot bottom water temperature for MetROMS, low-res FESOM, and high-res FESOM
# for the given region, in the given axes.
# Input:
# x_min, x_max, y_min, y_max = bounds on x and y (using polar coordinate
#               transformation from above) for the desired region
# gs = GridSpec object of size 1x3 to plot in
# cbaxes = Axes object for location of colourbar
# cbar_ticks = 1D array containing values for ticks on colourbar
# letter = 'a', 'b', 'c', etc. to add before the bottom water temp title, for
#          use in a figure showing multiple variables
def plot_bwtemp (x_min, x_max, y_min, y_max, gs, cbaxes, cbar_ticks, letter):

    # Set up a grey square for FESOM to fill the background with land
    x_reg_fesom, y_reg_fesom = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
    land_square = zeros(shape(x_reg_fesom))
    # Find bounds on variable in this region
    var_min, var_max = get_min_max(roms_bwtemp, fesom_bwtemp_lr, fesom_bwtemp_hr, x_min, x_max, y_min, y_max, cavity=False)

    # Plot MetROMS
    ax = subplot(gs[0,0], aspect='equal')
    # First shade land and zice in grey
    grey_cmap = ListedColormap([(0.6, 0.6, 0.6)])
    pcolor(roms_x, roms_y, land_zice, cmap=grey_cmap)
    # Fill in the missing circle
    pcolor(x_reg_roms, y_reg_roms, land_circle, cmap=grey_cmap)
    # Now shade the data
    pcolor(roms_x, roms_y, roms_bwtemp, vmin=var_min, vmax=var_max, cmap='jet')
    # Add a black contour for the ice shelf front
    rcParams['contour.negative_linestyle'] = 'solid'
    contour(roms_x, roms_y, zice_contour, levels=[min_zice], colors=('black'))
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    axis('off')

    # FESOM low-res
    ax = subplot(gs[0,1], aspect='equal')
    # Start with land background
    contourf(x_reg_fesom, y_reg_fesom, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    # Add ice shelf elements
    img = PatchCollection(patches_all_lr, cmap='jet')
    img.set_array(array(fesom_bwtemp_lr))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    # Mask out the open ocean in white
    #overlay = PatchCollection(mask_patches_lr, facecolor=(1,1,1))
    #overlay.set_edgecolor('face')
    #ax.add_collection(overlay)
    # Add ice shelf front contour lines
    fesom_contours_lr = LineCollection(contour_lines_lr, edgecolor='black', linewidth=1)
    ax.add_collection(fesom_contours_lr)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    axis('off')
    # Main title
    title(letter + r') Bottom water temperature ($^{\circ}$C)', fontsize=20)

    # FESOM high-res
    ax = subplot(gs[0,2], aspect='equal')
    contourf(x_reg_fesom, y_reg_fesom, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches_all_hr, cmap='jet')
    img.set_array(array(fesom_bwtemp_hr))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    #overlay = PatchCollection(mask_patches_hr, facecolor=(1,1,1))
    #overlay.set_edgecolor('face')
    #ax.add_collection(overlay)
    fesom_contours_hr = LineCollection(contour_lines_hr, edgecolor='black', linewidth=1)
    ax.add_collection(fesom_contours_hr)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    axis('off')
    # Colourbar on the right
    cbar = colorbar(img, cax=cbaxes, ticks=cbar_ticks)
    

# Plot bottom water salinity for MetROMS, low-res FESOM, and high-res FESOM
# for the given region, in the given axes.
# Input:
# x_min, x_max, y_min, y_max = bounds on x and y (using polar coordinate
#               transformation from above) for the desired region
# gs = GridSpec object of size 1x3 to plot in
# cbaxes = Axes object for location of colourbar
# cbar_ticks = 1D array containing values for ticks on colourbar
# letter = 'a', 'b', 'c', etc. to add before the bottom water salinity title,
#          for use in a figure showing multiple variables
def plot_bwsalt (x_min, x_max, y_min, y_max, gs, cbaxes, cbar_ticks, letter):

    # Set up a grey square for FESOM to fill the background with land
    x_reg_fesom, y_reg_fesom = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
    land_square = zeros(shape(x_reg_fesom))
    # Find bounds on variable in this region
    var_min, var_max = get_min_max(roms_bwsalt, fesom_bwsalt_lr, fesom_bwsalt_hr, x_min, x_max, y_min, y_max, cavity=False)

    # Plot MetROMS
    ax = subplot(gs[0,0], aspect='equal')
    # First shade land and zice in grey
    grey_cmap = ListedColormap([(0.6, 0.6, 0.6)])
    pcolor(roms_x, roms_y, land_zice, cmap=grey_cmap)
    # Fill in the missing circle
    pcolor(x_reg_roms, y_reg_roms, land_circle, cmap=grey_cmap)
    # Now shade the data
    pcolor(roms_x, roms_y, roms_bwsalt, vmin=var_min, vmax=var_max, cmap='jet')
    # Add a black contour for the ice shelf front
    rcParams['contour.negative_linestyle'] = 'solid'
    contour(roms_x, roms_y, zice_contour, levels=[min_zice], colors=('black'))
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    axis('off')

    # FESOM low-res
    ax = subplot(gs[0,1], aspect='equal')
    # Start with land background
    contourf(x_reg_fesom, y_reg_fesom, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    # Add ice shelf elements
    img = PatchCollection(patches_all_lr, cmap='jet')
    img.set_array(array(fesom_bwsalt_lr))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    # Mask out the open ocean in white
    #overlay = PatchCollection(mask_patches_lr, facecolor=(1,1,1))
    #overlay.set_edgecolor('face')
    #ax.add_collection(overlay)
    # Add ice shelf front contour lines
    fesom_contours_lr = LineCollection(contour_lines_lr, edgecolor='black', linewidth=1)
    ax.add_collection(fesom_contours_lr)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    axis('off')
    # Main title
    title(letter + ') Bottom water salinity (psu)', fontsize=20)

    # FESOM high-res
    ax = subplot(gs[0,2], aspect='equal')
    contourf(x_reg_fesom, y_reg_fesom, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches_all_hr, cmap='jet')
    img.set_array(array(fesom_bwsalt_hr))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    #overlay = PatchCollection(mask_patches_hr, facecolor=(1,1,1))
    #overlay.set_edgecolor('face')
    #ax.add_collection(overlay)
    fesom_contours_hr = LineCollection(contour_lines_hr, edgecolor='black', linewidth=1)
    ax.add_collection(fesom_contours_hr)
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    axis('off')
    # Colourbar on the right
    cbar = colorbar(img, cax=cbaxes, ticks=cbar_ticks)


# Plot vertically averaged velocity for MetROMS, low-res FESOM, and high-res
# FESOM for the given region, in the given axes.
# Input:
# x_min, x_max, y_min, y_max = bounds on x and y (using polar coordinate
#               transformation from above) for the desired region
# gs = GridSpec object of size 1x3 to plot in
# cbaxes = Axes object for location of colourbar
# cbar_ticks = 1D array containing values for ticks on colourbar
# x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr,
#            fesom_ubin_hr, fesom_vbin_hr = output variables from the 
#                           "make_vectors" function for the given region
# letter = 'a', 'b', 'c', etc. to add before the bottom water salinity title,
#          for use in a figure showing multiple variables
# loc_string = optional string containing location title, in the case of
#              velocity being zoomed into a smaller region to the other
#              variables (eg "George VI Ice Shelf" for the Bellingshausen plot)
# arrow_scale, arrow_headwidth, arrow_headlength = optional parameters for
#              arrows on vector overlay
def plot_velavg (x_min, x_max, y_min, y_max, gs, cbaxes, cbar_ticks, x_centres, y_centres, roms_ubin, roms_vbin, fesom_ubin_lr, fesom_vbin_lr, fesom_ubin_hr, fesom_vbin_hr, letter, loc_string=None, arrow_scale=0.9, arrow_headwidth=8, arrow_headlength=9):

    # Set up a grey square for FESOM to fill the background with land
    x_reg_fesom, y_reg_fesom = meshgrid(linspace(x_min, x_max, num=100), linspace(y_min, y_max, num=100))
    land_square = zeros(shape(x_reg_fesom))
    # Find bounds on variable in this region
    var_min, var_max = get_min_max(roms_speed, fesom_speed_lr, fesom_speed_hr, x_min, x_max, y_min, y_max)

    # Plot MetROMS
    ax = subplot(gs[0,0], aspect='equal')
    # First shade land and zice in grey
    grey_cmap = ListedColormap([(0.6, 0.6, 0.6)])
    pcolor(roms_x, roms_y, land_zice, cmap=grey_cmap)
    # Fill in the missing circle
    pcolor(x_reg_roms, y_reg_roms, land_circle, cmap=grey_cmap)
    # Now shade the data
    pcolor(roms_x, roms_y, roms_speed, vmin=var_min, vmax=var_max, cmap='cool')
    # Overlay vectors
    quiver(x_centres, y_centres, roms_ubin, roms_vbin, scale=arrow_scale, headwidth=arrow_headwidth, headlength=arrow_headlength, color='black')
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    axis('off')

    # FESOM low-res
    ax = subplot(gs[0,1], aspect='equal')
    # Start with land background
    contourf(x_reg_fesom, y_reg_fesom, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    # Add ice shelf elements
    img = PatchCollection(patches_lr, cmap='cool')
    img.set_array(array(fesom_speed_lr))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    # Mask out the open ocean in white
    overlay = PatchCollection(mask_patches_lr, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    # Overlay vectors
    quiver(x_centres, y_centres, fesom_ubin_lr, fesom_vbin_lr, scale=arrow_scale, headwidth=arrow_headwidth, headlength=arrow_headlength, color='black')
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    axis('off')
    # Main title
    if loc_string is None:
        title(letter + ') Vertically averaged ocean velocity (m/s)', fontsize=20)
    else:
        title(letter + ') Vertically averaged ocean velocity (m/s): '+loc_string, fontsize=15)

    # FESOM high-res
    ax = subplot(gs[0,2], aspect='equal')
    contourf(x_reg_fesom, y_reg_fesom, land_square, 1, colors=(('0.6', '0.6', '0.6')))
    img = PatchCollection(patches_hr, cmap='cool')
    img.set_array(array(fesom_speed_hr))
    img.set_edgecolor('face')
    img.set_clim(vmin=var_min, vmax=var_max)
    ax.add_collection(img)
    overlay = PatchCollection(mask_patches_hr, facecolor=(1,1,1))
    overlay.set_edgecolor('face')
    ax.add_collection(overlay)
    quiver(x_centres, y_centres, fesom_ubin_hr, fesom_vbin_hr, scale=arrow_scale, headwidth=arrow_headwidth, headlength=arrow_headlength, color='black')
    xlim([x_min, x_max])
    ylim([y_min, y_max])
    axis('off')
    # Colourbar on the right
    cbar = colorbar(img, cax=cbaxes, ticks=cbar_ticks)


# Plot zonal slices (latitude vs depth) of temperature and salinity for
# MetROMS, low-res FESOM, and high-res FESOM through the given longitude, in the
# given axes.
# Input:
# lon0 = longitude to plot (-180 to 180)
# lat_min, lat_max = bounds on latitude to plot (-90 to 90)
# lat_ticks = array of floats containing latitude values to tick
# lat_labels = array of strings containing label for each tick
# depth_min, depth_max = bounds on depth to plot (negative, in metres)
# depth_ticks, depth_labels = tick locations and labels as for latitude
# gs1, gs2 = GridSpec objects of size 1x3, in which to plot temperature and
#            salinity respectively
# cbaxes1, cbaxes2 = Axes objects for locations of temperature and salinity
#                    colourbars
# cbar_ticks1, cbar_ticks2 = 1D arrays containing values for ticks on
#                            temperature and salinity colourbars
# letter1, letter2 = 'a', 'b', 'c', etc. to add before the temperature and 
#                    salinity titles
# loc_string = optional string containing location title, in the case of
#              plotting slices through a single ice shelf in a large region
#              (eg "Fimbul Ice Shelf" for the Eastern Weddell plot)
def plot_zonal_ts (lon0, lat_min, lat_max, lat_ticks, lat_labels, depth_min, depth_max, depth_ticks, depth_labels, gs1, gs2, cbaxes1, cbaxes2, cbar_ticks1, cbar_ticks2, letter1, letter2, loc_string=None):

    # Figure out what to write on the title about longitude
    if lon0 < 0:
        lon_string = str(-lon0)+r'$^{\circ}$W'
    else:
        lon_string = str(lon0)+r'$^{\circ}$E'

    # ROMS
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
    roms_temp, roms_z_1d, roms_lat_1d = interp_lon_roms(roms_temp_3d, roms_z_3d, roms_lat, roms_lon, lon0)
    roms_salt, roms_z_1d, roms_lat_1d = interp_lon_roms(roms_salt_3d, roms_z_3d, roms_lat, roms_lon, lon0)
    # Switch back to range -180-180
    if lon0 > 180:
        lon0 -= 360

    # FESOM low-res
    # Read temperature and salinity
    id = Dataset(fesom_file_lr_o, 'r')
    fesom_temp_nodes_lr = id.variables['temp'][0,:]
    fesom_salt_nodes_lr = id.variables['salt'][0,:]
    id.close()
    # Build arrays of SideElements making up zonal slices
    selements_temp_lr = fesom_sidegrid(elm2D_lr, fesom_temp_nodes_lr, lon0, lat_max)
    selements_salt_lr = fesom_sidegrid(elm2D_lr, fesom_salt_nodes_lr, lon0, lat_max)
    # Build array of quadrilateral patches for the plots, and data values
    # corresponding to each SideElement
    spatches_lr = []
    fesom_temp_lr = []
    for selm in selements_temp_lr:
        # Make patch
        coord = transpose(vstack((selm.y, selm.z)))
        spatches_lr.append(Polygon(coord, True, linewidth=0.))
        # Save data value
        fesom_temp_lr.append(selm.var)
    fesom_temp_lr = array(fesom_temp_lr)
    # Salinity has same patches but different values
    fesom_salt_lr = []
    for selm in selements_salt_lr:
        fesom_salt_lr.append(selm.var)
    fesom_salt_lr = array(fesom_salt_lr)

    # FESOM high-res
    id = Dataset(fesom_file_hr_o, 'r')
    fesom_temp_nodes_hr = id.variables['temp'][0,:]
    fesom_salt_nodes_hr = id.variables['salt'][0,:]
    id.close()
    selements_temp_hr = fesom_sidegrid(elm2D_hr, fesom_temp_nodes_hr, lon0, lat_max)
    selements_salt_hr = fesom_sidegrid(elm2D_hr, fesom_salt_nodes_hr, lon0, lat_max)
    spatches_hr = []
    fesom_temp_hr = []
    for selm in selements_temp_hr:
        coord = transpose(vstack((selm.y, selm.z)))
        spatches_hr.append(Polygon(coord, True, linewidth=0.))
        fesom_temp_hr.append(selm.var)
    fesom_temp_hr = array(fesom_temp_hr)
    fesom_salt_hr = []
    for selm in selements_salt_hr:
        fesom_salt_hr.append(selm.var)
    fesom_salt_hr = array(fesom_salt_hr)

    # Find bounds on each variable
    flag = (roms_lat_1d >= lat_min)*(roms_lat_1d <= lat_max)
    temp_min = amin(array([amin(roms_temp[flag]), amin(fesom_temp_lr), amin(fesom_temp_hr)]))
    temp_max = amax(array([amax(roms_temp[flag]), amax(fesom_temp_lr), amax(fesom_temp_hr)]))
    salt_min = amin(array([amin(roms_salt[flag]), amin(fesom_salt_lr), amin(fesom_salt_hr)]))
    salt_max = amax(array([amax(roms_salt[flag]), amax(fesom_salt_lr), amax(fesom_salt_hr)]))

    # Plot temperature
    # MetROMS
    ax = subplot(gs1[0,0])
    pcolor(roms_lat_1d, roms_z_1d, roms_temp, vmin=temp_min, vmax=temp_max, cmap='jet')
    ylabel('Depth (m)', fontsize=12)
    xlabel('Latitude', fontsize=12)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels)

    # FESOM low-res
    ax = subplot(gs1[0,1])
    img = PatchCollection(spatches_lr, cmap='jet')
    img.set_array(fesom_temp_lr)
    img.set_edgecolor('face')
    img.set_clim(vmin=temp_min, vmax=temp_max)
    ax.add_collection(img)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if loc_string is None:
        title(letter1 + r') Temperature ($^{\circ}$C) through ' + lon_string, fontsize=20)
    else:
        title(letter1 + r') Temperature ($^{\circ}$C) through ' + lon_string + ': ' + loc_string, fontsize=20)

    # FESOM high-res
    ax = subplot(gs1[0,2])
    img = PatchCollection(spatches_hr, cmap='jet')
    img.set_array(fesom_temp_hr)
    img.set_edgecolor('face')
    img.set_clim(vmin=temp_min, vmax=temp_max)
    ax.add_collection(img)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    # Add a colorbar
    cbar = colorbar(img, cax=cbaxes1, ticks=cbar_ticks1)

    # Plot salinity
    # MetROMS
    ax = subplot(gs2[0,0])
    pcolor(roms_lat_1d, roms_z_1d, roms_salt, vmin=salt_min, vmax=salt_max, cmap='jet')
    ylabel('Depth (m)', fontsize=12)
    xlabel('Latitude', fontsize=12)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(lat_labels)
    ax.set_yticks(depth_ticks)
    ax.set_yticklabels(depth_labels)
    # FESOM low-res
    ax = subplot(gs2[0,1])
    img = PatchCollection(spatches_lr, cmap='jet')
    img.set_array(fesom_salt_lr)
    img.set_edgecolor('face')
    img.set_clim(vmin=salt_min, vmax=salt_max)
    ax.add_collection(img)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if loc_string is None:
        title(letter2 + ') Salinity (psu) through ' + lon_string, fontsize=20)
    else:
        title(letter2 + ') Salinity (psu) through ' + lon_string + ': ' + loc_string, fontsize=20)
    # FESOM high-res
    ax = subplot(gs2[0,2])
    img = PatchCollection(spatches_hr, cmap='jet')
    img.set_array(fesom_salt_hr)
    img.set_edgecolor('face')
    img.set_clim(vmin=salt_min, vmax=salt_max)
    ax.add_collection(img)
    xlim([lat_min, lat_max])
    ylim([depth_min, depth_max])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    # Add a colorbar
    cbar = colorbar(img, cax=cbaxes2, ticks=cbar_ticks2)
    
    
    
    

    
