from netCDF4 import Dataset
from numpy import *

def find_isolated_points (cice_kmt_file):

    start_j = 50
    end_j = 250

    id = Dataset(cice_kmt_file, 'r')
    kmt = id.variables['kmt'][:,:]
    id.close()

    num_i = size(kmt,1)

    for j in range(start_j, end_j):
        for i in range(num_i-1):
            if kmt[j,i] == 1:
                if i == num_i-1:
                    neighbours = array([kmt[j,i-1], kmt[j,0], kmt[j-1,i], kmt[j+1,i]])
                else:
                    neighbours = array([kmt[j,i-1], kmt[j,i+1], kmt[j-1,i], kmt[j+1,i]])
                if sum(neighbours) < 2:
                    print "i=" + str(i+1) + ', j=' + str(j+1)


if __name__ == "__main__":

    cice_kmt_file = raw_input("Path to CICE land mask file: ")
    find_isolated_points(cice_kmt_file)
