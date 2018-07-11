from numpy import *
from netCDF4 import Dataset

from cartesian_grid_2d import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from fesom_grid import *

def total_iceshelf_area (roms_grid_file, fesom_mesh_path_lr, fesom_mesh_path_hr):

    id = Dataset(roms_grid_file, 'r')
    lon = id.variables['lon_rho'][:-15,1:-1]
    lat = id.variables['lat_rho'][:-15,1:-1]
    zice = id.variables['zice'][:-15,1:-1]
    id.close()
    dx, dy = cartesian_grid_2d(lon, lat)
    dA = ma.masked_where(zice==0, dx*dy)
    print 'MetROMS: ' + str(sum(dA)) + ' m^2'

    elements_lr = fesom_grid(fesom_mesh_path_lr, circumpolar=True, cross_180=False)
    area_elm_lr = zeros(len(elements_lr))
    for i in range(len(elements_lr)):
        elm = elements_lr[i]
        if elm.cavity:
            area_elm_lr[i] = elm.area()
    print 'FESOM (low-res): ' + str(sum(area_elm_lr)) + ' m^2'

    elements_hr = fesom_grid(fesom_mesh_path_hr, circumpolar=True, cross_180=False)
    area_elm_hr = zeros(len(elements_hr))
    for i in range(len(elements_hr)):
        elm = elements_hr[i]
        if elm.cavity:
            area_elm_hr[i] = elm.area()
    print 'FESOM (high-res): ' + str(sum(area_elm_hr)) + ' m^2'


if __name__ == "__main__":

    roms_grid_file = raw_input('Enter path to ROMS grid file: ')
    fesom_mesh_path_lr = raw_input('Enter path to FESOM low-res mesh directory: ')
    fesom_mesh_path_hr = raw_input('Enter path to FESOM high-res mesh directory: ')
    total_iceshelf_area (roms_grid_file, fesom_mesh_path_lr, fesom_mesh_path_hr)
    
    
