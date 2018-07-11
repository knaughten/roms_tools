from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from scipy.interpolate import griddata
from cartesian_grid_2d import *
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from fesom_grid import *

def mip_seaice_tamura ():

    # File paths
    # ROMS grid (just for bathymetry)
    roms_grid = '/short/m68/kaa561/metroms_iceshelf/apps/common/grid/circ30S_quarterdegree.nc'
    # FESOM mesh paths
    fesom_mesh_path_lr = '/short/y99/kaa561/FESOM/mesh/meshA/'
    fesom_mesh_path_hr = '/short/y99/kaa561/FESOM/mesh/meshB/'
    # CICE 1992-2013 mean ice production (precomputed in calc_ice_prod.py)
    cice_file = '/short/m68/kaa561/metroms_iceshelf/tmproms/run/intercomparison/ice_prod_1992_2013.nc'
    # FESOM 1992-2013 mean ice production (precomputed in calc_annual_ice_prod.py in fesomtools)
    fesom_lr_file = '/short/y99/kaa561/FESOM/intercomparison_lowres/output/ice_prod_1992_2013.nc'
    fesom_hr_file = '/short/y99/kaa561/FESOM/intercomparison_highres/output/ice_prod_1992_2013.nc'
    # Tamura's 1992-2013 mean ice production (precomputed on desktop with Matlab)
    tamura_file = '/short/m68/kaa561/tamura_1992_2013_monthly_climatology.nc'
    # Output ASCII file
    output_file = 'seaice_prod_bins.log'
    # Size of longitude bin
    dlon_bin = 1.0
    # Definition of continental shelf: everywhere south of lat0 with
    # bathymetry shallower than h0
    lat0 = -60
    h0 = 1500
    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0

    # Set up longitude bins
    bin_edges = arange(-180, 180+dlon_bin, dlon_bin)
    bin_centres = 0.5*(bin_edges[:-1] + bin_edges[1:])
    num_bins = len(bin_centres)

    print 'Processing MetROMS'
    # Read CICE grid
    id = Dataset(cice_file, 'r')
    cice_lon = id.variables['TLON'][:,:]
    cice_lat = id.variables['TLAT'][:,:]
    # Read sea ice production
    cice_data = id.variables['ice_prod'][:,:]
    id.close()
    # Get area integrands
    dx, dy = cartesian_grid_2d(cice_lon, cice_lat)
    dA = dx*dy
    # Make sure longitude is in the range [-180, 180]
    index = cice_lon > 180
    cice_lon[index] = cice_lon[index] - 360
    # Read bathymetry (ROMS grid file) and trim to CICE grid
    id = Dataset(roms_grid, 'r')
    cice_bathy = id.variables['h'][1:-1,1:-1]
    id.close()
    # Set up integral
    cice_data_bins = zeros(num_bins)
    # Loop over all cells
    num_lon = size(cice_lon,1)
    num_lat = size(cice_lat,0)
    for j in range(num_lat):
        for i in range(num_lon):
            # Check for land mask or ice shelves
            if cice_data[j,i] is ma.masked:
                continue
            # Check for continental shelf
            if cice_lat[j,i] < lat0 and cice_bathy[j,i] < h0:
                # Find the right bin
                bin_index = nonzero(bin_edges > cice_lon[j,i])[0][0] - 1
                # Integrate (m^3/y)
                cice_data_bins[bin_index] += cice_data[j,i]*dA[j,i]
    # Convert to 10^9 m^3/y
    cice_data_bins *= 1e-9

    print 'Processing low-res FESOM'
    # Build mesh
    elements_lr = fesom_grid(fesom_mesh_path_lr, circumpolar=True, cross_180=False)
    # Read sea ice production
    id = Dataset(fesom_lr_file, 'r')
    fesom_data_lr = id.variables['ice_prod'][:]
    id.close()
    # Set up integral
    fesom_data_bins_lr = zeros(num_bins)
    # Loop over elements
    for elm in elements_lr:
        # Exclude ice shelf cavities
        if not elm.cavity:
            # Check for continental shelf in 2 steps
            if all(elm.lat < lat0):
                elm_bathy = mean([elm.nodes[0].find_bottom().depth, elm.nodes[1].find_bottom().depth, elm.nodes[2].find_bottom().depth])
                if elm_bathy < h0:
                    # Get element-averaged sea ice production
                    elm_data = mean([fesom_data_lr[elm.nodes[0].id], fesom_data_lr[elm.nodes[1].id], fesom_data_lr[elm.nodes[2].id]])
                    # Find the right bin
                    elm_lon = mean(elm.lon)
                    if elm_lon < -180:
                        elm_lon += 360
                    elif elm_lon > 180:
                        elm_lon -= 360
                    bin_index = nonzero(bin_edges > elm_lon)[0][0] - 1
                    # Integrate (m^3/y)
                    fesom_data_bins_lr[bin_index] += elm_data*elm.area()
    # Convert to 10^9 m^3/y
    fesom_data_bins_lr *= 1e-9

    print 'Processing high-res FESOM'
    elements_hr = fesom_grid(fesom_mesh_path_hr, circumpolar=True, cross_180=False)
    id = Dataset(fesom_hr_file, 'r')
    fesom_data_hr = id.variables['ice_prod'][:]
    id.close()
    fesom_data_bins_hr = zeros(num_bins)
    for elm in elements_hr:
        if not elm.cavity:
            if all(elm.lat < lat0):
                elm_bathy = mean([elm.nodes[0].find_bottom().depth, elm.nodes[1].find_bottom().depth, elm.nodes[2].find_bottom().depth])
                if elm_bathy < h0:
                    elm_data = mean([fesom_data_hr[elm.nodes[0].id], fesom_data_hr[elm.nodes[1].id], fesom_data_hr[elm.nodes[2].id]])
                    elm_lon = mean(elm.lon)
                    if elm_lon < -180:
                        elm_lon += 360
                    elif elm_lon > 180:
                        elm_lon -= 360
                    bin_index = nonzero(bin_edges > elm_lon)[0][0] - 1
                    fesom_data_bins_hr[bin_index] += elm_data*elm.area()
    fesom_data_bins_hr *= 1e-9

    print 'Processing Tamura obs'
    id = Dataset(tamura_file, 'r')
    # Read grid and data
    tamura_lon = id.variables['longitude'][:,:]
    tamura_lat = id.variables['latitude'][:,:]
    # Read sea ice formation
    tamura_data = id.variables['ice_prod'][:,:]
    id.close()
    # Interpolate to a regular grid so we can easily integrate over area
    dlon_reg = 0.2
    dlat_reg = 0.1
    lon_reg_edges = arange(-180, 180+dlon_reg, dlon_reg)
    lon_reg = 0.5*(lon_reg_edges[:-1] + lon_reg_edges[1:])
    lat_reg_edges = arange(-80, -60+dlat_reg, dlat_reg)
    lat_reg = 0.5*(lat_reg_edges[:-1] + lat_reg_edges[1:])
    lon_reg_2d, lat_reg_2d = meshgrid(lon_reg, lat_reg)
    dx_reg = r*cos(lat_reg_2d*deg2rad)*dlon_reg*deg2rad
    dy_reg = r*dlat_reg*deg2rad
    dA_reg = dx_reg*dy_reg
    # Be careful with the periodic boundary here
    num_pts = size(tamura_lon)
    num_wrap1 = count_nonzero(tamura_lon < -179)
    num_wrap2 = count_nonzero(tamura_lon > 179)
    points = empty([num_pts+num_wrap1+num_wrap2,2])
    values = empty(num_pts+num_wrap1+num_wrap2)
    points[:num_pts,0] = ravel(tamura_lon)
    points[:num_pts,1] = ravel(tamura_lat)
    values[:num_pts] = ravel(tamura_data)
    # Wrap the periodic boundary on both sides
    index = tamura_lon < -179
    points[num_pts:num_pts+num_wrap1,0] = tamura_lon[index] + 360
    points[num_pts:num_pts+num_wrap1,1] = tamura_lat[index]
    values[num_pts:num_pts+num_wrap1] = tamura_data[index]
    index = tamura_lon > 179
    points[num_pts+num_wrap1:,0] = tamura_lon[index] - 360
    points[num_pts+num_wrap1:,1] = tamura_lat[index]
    values[num_pts+num_wrap1:] = tamura_data[index]
    values = ma.masked_where(isnan(values), values)
    xi = empty([size(lon_reg_2d),2])
    xi[:,0] = ravel(lon_reg_2d)
    xi[:,1] = ravel(lat_reg_2d)
    result = griddata(points, values, xi)
    tamura_data_reg = reshape(result, shape(lon_reg_2d))
    # Now, regrid the MetROMS bathymetry to this regular grid
    num_pts = size(cice_lon)
    num_wrap1 = count_nonzero(cice_lon < -179)
    num_wrap2 = count_nonzero(cice_lon > 179)
    points = empty([num_pts+num_wrap1+num_wrap2,2])
    values = empty(num_pts+num_wrap1+num_wrap2)
    points[:num_pts,0] = ravel(cice_lon)
    points[:num_pts,1] = ravel(cice_lat)
    values[:num_pts] = ravel(cice_bathy)
    index = cice_lon < -179
    points[num_pts:num_pts+num_wrap1,0] = cice_lon[index] + 360
    points[num_pts:num_pts+num_wrap1,1] = cice_lat[index]
    values[num_pts:num_pts+num_wrap1] = cice_bathy[index]
    index = cice_lon > 179
    points[num_pts+num_wrap1:,0] = cice_lon[index] - 360
    points[num_pts+num_wrap1:,1] = cice_lat[index]
    values[num_pts+num_wrap1:] = cice_bathy[index]
    values = ma.masked_where(isnan(values), values)
    xi = empty([size(lon_reg_2d),2])
    xi[:,0] = ravel(lon_reg_2d)
    xi[:,1] = ravel(lat_reg_2d)
    result = griddata(points, values, xi)
    bathy_reg = reshape(result, shape(lon_reg_2d))
    # Mask everything but the continental shelf from dA_reg
    dA_reg = ma.masked_where(lat_reg_2d > lat0, dA_reg)
    dA_reg = ma.masked_where(bathy_reg > h0, dA_reg)
    # Mask the land mask (and ice shelves) from tamura_data_reg
    tamura_data_reg = ma.masked_where(isnan(tamura_data_reg), tamura_data_reg)
    # Set up integral
    tamura_data_bins = zeros(num_bins)
    # Loop over longitude only
    for i in range(len(lon_reg)):
        # Find the right bin
        bin_index = nonzero(bin_edges > lon_reg[i])[0][0] - 1
        # Integrate (m^3/y)
        tamura_data_bins[bin_index] += sum(tamura_data_reg[:,i]*dA_reg[:,i])
    # Convert to 10^9 m^3/y
    tamura_data_bins *= 1e-9

    # Write data to ASCII file
    print 'Writing to file'
    f = open(output_file, 'w')
    f.write('Longitude:\n')
    for val in bin_centres:
        f.write(str(val) + '\n')
    f.write('MetROMS sea ice production (10^9 m^3/y):\n')
    for val in cice_data_bins:
        f.write(str(val) + '\n')
    f.write('FESOM (low-res) sea ice production (10^9 m^3/y):\n')
    for val in fesom_data_bins_lr:
        f.write(str(val) + '\n')
    f.write('FESOM (high-res) sea ice production (10^9 m^3/y):\n')
    for val in fesom_data_bins_hr:
        f.write(str(val) + '\n')
    f.write('Tamura sea ice production (10^9 m^3/y):\n')
    for val in tamura_data_bins:
        f.write(str(val) + '\n')
    f.close()
    
    
