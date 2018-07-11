from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from cartesian_grid_3d import *
from rotate_vector_roms import *
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from unrotate_grid import *
from unrotate_vector import *

def slope_current ():

# File paths
roms_grid = '/short/m68/kaa561/metroms_iceshelf/apps/common/grid/circ30S_quarterdegree.nc'
roms_file = '/short/m68/kaa561/metroms_iceshelf/tmproms/run/intercomparison/2002_2016_avg.nc'
fesom_mesh_path_lr = '/short/y99/kaa561/FESOM/mesh/meshA/'
fesom_mesh_path_hr = '/short/y99/kaa561/FESOM/mesh/meshB/'
fesom_file_lr = '/short/y99/kaa561/FESOM/intercomparison_lowres/output/oce_2002_2016_avg.nc'
fesom_file_hr = '/short/y99/kaa561/FESOM/intercomparison_highres/output/oce_2002_2016_avg.nc'
# ROMS vertical grid parameters
theta_s = 7.0
theta_b = 2.0
hc = 250
N = 31
# FESOM mesh parameters
circumpolar = False
cross_180 = False
# Spacing of longitude bins
dlon = 1
# Parameters for continental shelf selection
lat0 = -64  # Maximum latitude to consider
h0 = 2500  # Deepest depth to consider

# Set up longitude bins
# Start with edges
lon_bins = arange(-180, 180+dlon, dlon)
# Centres for plotting
lon_centres = 0.5*(lon_bins[:-1] + lon_bins[1:])
num_bins = size(lon_centres)
# Set up arrays to store maximum barotropic speed in each bin
current_roms = zeros(num_bins)
current_fesom_lr = zeros(num_bins)
current_fesom_hr = zeros(num_bins)

print 'Processing MetROMS'

print 'Reading grid'
id = Dataset(roms_grid, 'r')
roms_lon = id.variables['lon_rho'][:,:]
roms_lat = id.variables['lat_rho'][:,:]
roms_h = id.variables['h'][:,:]
roms_zice = id.variables['zice'][:,:]
roms_angle = id.variables['angle'][:,:]
id.close()
print 'Reading data'
# Read full 3D u and v
id = Dataset(roms_file, 'r')
u_3d_tmp = id.variables['u'][0,:,:,:]
v_3d_tmp = id.variables['v'][0,:,:,:]
id.close()
print 'Vertically averaging velocity'
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
    u_k, v_k = rotate_vector_roms(u_3d_tmp[k,:,:], v_3d_tmp[k,:,:], roms_angle)
    u_3d[k,:,:] = u_k
    v_3d[k,:,:] = v_k
# Vertically average u and v
roms_u = sum(u_3d*dz, axis=0)/sum(dz, axis=0)
roms_v = sum(v_3d*dz, axis=0)/sum(dz, axis=0)
# Calculate speed
roms_speed = sqrt(roms_u**2 + roms_v**2)
print 'Selecting slope current'
# First make sure longitude is between -180 and 180
index = roms_lon > 180
roms_lon[index] = roms_lon[index] - 360
for j in range(size(roms_speed,0)):
    for i in range(size(roms_speed,1)):
        # Check if we care about this point
        if roms_lat[j,i] <= lat0 and roms_h[j,i] <= h0 and roms_zice[j,i] == 0:
            # Find longitude bin
            lon_index = nonzero(lon_bins > roms_lon[j,i])[0][0] - 1
            # Update slope current speed in this bin if needed
            if roms_speed[j,i] > current_roms[lon_index]:
                current_roms[lon_index] = roms_speed[j,i]

print 'Processing low-res FESOM'

print 'Building mesh'
# We only care about nodes, not elements, so don't need to use the
# fesom_grid function.
# Read cavity flag for each 2D surface node
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
# Read lists of which nodes are directly below which
f = open(fesom_mesh_path_lr + 'aux3d.out', 'r')
max_num_layers_lr = int(f.readline())
node_columns_lr = zeros([fesom_n2d_lr, max_num_layers_lr])
for n in range(fesom_n2d_lr):
    for k in range(max_num_layers_lr):
        node_columns_lr[n,k] = int(f.readline())
node_columns_lr = node_columns_lr.astype(int)
f.close()
# Now figure out the bottom depth of each 2D node
bottom_depth_lr = zeros(fesom_n2d_lr)
for n in range(fesom_n2d_lr):
    node_id = node_columns_lr[n,0] - 1
    for k in range(1, max_num_layers_lr):
        if node_columns_lr[n,k] == -999:
            # Reached the bottom
            break
        node_id = node_columns_lr[n,k] - 1
    # Save the last valid depth
    bottom_depth_lr[n] = node_depth_lr[n]
print 'Reading data'
# Read full 3D field for both u and v
id = Dataset(fesom_file_lr, 'r')
node_ur_3d_lr = id.variables['u'][0,:]
node_vr_3d_lr = id.variables['v'][0,:]
id.close()
print 'Vertically averaging velocity'
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
print 'Selecting slope current'
for n in range(fesom_n2d_lr):
    # Check if we care about this node
    if fesom_lat_lr[n] <= lat0 and bottom_depth_lr[n] <= h0 and not fesom_cavity_lr[n]:
        # Find longitude bin
        lon_index = nonzero(lon_bins > fesom_lon_lr[n])[0][0] - 1
        # Update slope current speed in this bin if needed
        if node_speed_lr[n] > current_fesom_lr[lon_index]:
            current_fesom_lr[lon_index] = node_speed_lr[n]

print 'Processing high-res FESOM'

print 'Building mesh'
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
f = open(fesom_mesh_path_hr + 'aux3d.out', 'r')
max_num_layers_hr = int(f.readline())
node_columns_hr = zeros([fesom_n2d_hr, max_num_layers_hr])
for n in range(fesom_n2d_hr):
    for k in range(max_num_layers_hr):
        node_columns_hr[n,k] = int(f.readline())
node_columns_hr = node_columns_hr.astype(int)
f.close()
bottom_depth_hr = zeros(fesom_n2d_hr)
for n in range(fesom_n2d_hr):
    node_id = node_columns_hr[n,0] - 1
    for k in range(1, max_num_layers_hr):
        if node_columns_hr[n,k] == -999:
            break
        node_id = node_columns_hr[n,k] - 1
    bottom_depth_hr[n] = node_depth_hr[n]
print 'Reading data'
id = Dataset(fesom_file_hr, 'r')
node_ur_3d_hr = id.variables['u'][0,:]
node_vr_3d_hr = id.variables['v'][0,:]
id.close()
print 'Vertically averaging velocity'
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
print 'Selecting slope current'
for n in range(fesom_n2d_hr):
    if fesom_lat_hr[n] <= lat0 and bottom_depth_hr[n] <= h0 and not fesom_cavity_hr[n]:
        lon_index = nonzero(lon_bins > fesom_lon_hr[n])[0][0] - 1
        if node_speed_hr[n] > current_fesom_hr[lon_index]:
            current_fesom_hr[lon_index] = node_speed_hr[n]

print 'Plotting'
fig = figure(figsize=(12,8))
plot(lon_centres, current_roms, color='blue', label='MetROMS')
plot(lon_centres, current_fesom_lr, color='green', label='FESOM low-res')
plot(lon_centres, current_fesom_hr, color='magenta', label='FESOM high-res')
grid(True)
title('Slope current speed', fontsize=20)
xlabel('Longitude', fontsize=14)
ylabel('m/s', fontsize=14)
xlim([-180, 180])
legend()
fig.savefig('slope_current.png')

print 'Mean slope current in MetROMS: ' + str(mean(current_roms)) + ' m/s'
print 'Mean slope current in low-res FESOM: ' + str(mean(current_fesom_lr)) + ' m/s'
print 'Mean slope current in high-res FESOM: ' + str(mean(current_fesom_hr)) + ' m/s'


# Command-line interface
if __name__ == "__main__":

    coastal_current()
    
        
    
    
