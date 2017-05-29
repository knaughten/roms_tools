from numpy import *

# Given ROMS grid variables, calculate s-coordinates, stretching curves, and
# z-coordinates. Assumes Vtransform = 2.

# Input (all straight out of grid file and *.in file):
# h, zice = 2D arrays containing values for bathymetry and ice shelf draft.
#           Both have dimension latitude x longitude.
# theta_s, theta_b, hc, N = scalar parameters
# zeta = optional 2D array containing values for sea surface height
# Vstretching = optional integer showing stretching transfomration, 2 or 4
# Output:
# z = 3D array containing negative z-coordinate values for depth on the rho 
#     grid; dimension depth x latitude x longitude
# s = 1D array of s-coordinate values
# C = 1D array of stretching curve values
def calc_z (h, zice, theta_s, theta_b, hc, N, zeta=None, Vstretching=4):

    # Follows the method of scoord_zice.m and stretching.m on katabatic
    # (in /ds/projects/iomp/matlab_scripts/ROMS_NetCDF/iomp_IAF/)
    # which is also explained on the ROMS wiki:
    # https://www.myroms.org/wiki/Vertical_S-coordinate.

    alpha = 1.0
    beta = 1.0
    ds = 1.0/N
    lev = arange(1,N+1)-0.5
    s = (lev-N)*ds

    if Vstretching == 2:
        Csur = (-cosh(theta_s*s) + 1)/(cosh(theta_s) - 1)
        Cbot = sinh(theta_b*(s+1))/sinh(theta_b) - 1
        weight = (s+1)**alpha*(1 + (alpha/beta)*(1 - (s+1)**beta))
        C = weight*Csur + (1-weight)*Cbot
    elif Vstretching == 4:
        C = (1-cosh(theta_s*s))/(cosh(theta_s)-1)
        C = (exp(theta_b*C)-1)/(1-exp(-theta_b))
        
    h = h - abs(zice)

    num_lon = size(h, 1)
    num_lat = size(h, 0)
    z = zeros((N, num_lat, num_lon))
    for k in range(N):
        z0 = (h*C[k] + hc*s[k])/(h + hc)
        if zeta is None:
            z[k,:,:] = h*z0 - abs(zice)
        else:
            z[k,:,:] = (zeta+h)*z0 + zeta - abs(zice)

    return z, s, C
