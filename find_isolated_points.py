from netCDF4 import Dataset
from numpy import *

# Find all of the CICE grid points which are land (or ice shelf) on 3 sides.
# Sea ice can grow in these isolated points but cannot escape due to CICE's
# coastal boundary conditions, so it gets crazy thick (like 2 km thick).
# Print the indices of these points to the screen. This script assumes a
# periodic boundary in the longitude direction.
# Input: cice_kmt_file = path to CICE land mask file, created using cice_grid.py
def find_isolated_points (cice_kmt_file):

    # j-indices that might have sea ice
    start_j = 50
    end_j = 250

    # Read land mask
    id = Dataset(cice_kmt_file, 'r')
    kmt = id.variables['kmt'][:,:]
    id.close()

    num_i = size(kmt,1)

    # Double loop, can't find a cleaner way to do this
    for j in range(start_j, end_j):
        for i in range(num_i):
            # Check for unmasked points
            if kmt[j,i] == 1:
                # Count the number of neighbours which are unmasked
                if i == num_i-1:
                    # Loop back to the beginning for neighbour on the left
                    neighbours = array([kmt[j,i-1], kmt[j,0], kmt[j-1,i], kmt[j+1,i]])
                else:
                    neighbours = array([kmt[j,i-1], kmt[j,i+1], kmt[j-1,i], kmt[j+1,i]])
                if sum(neighbours) < 2:
                    # Blocked on at least 3 sides
                    print "i=" + str(i+1) + ', j=' + str(j+1)


# Command-line interface
if __name__ == "__main__":

    cice_kmt_file = raw_input("Path to CICE land mask file: ")
    find_isolated_points(cice_kmt_file)
