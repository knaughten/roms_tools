from numpy import *

# Compute the UNESCO equation of state for seawater density
# given temp = temperature (C), salt = salinity (psu), and press = pressure 
# (bar - note that p = depth (in m) / 10 is acceptable). Fully vectorized.
# See http://unesdoc.unesco.org/images/0004/000461/046148eb.pdf; all
# intermediate variables are named the same as in this document.

# Input:
# temp = array of any dimension, containing temperature values (C)
# salt = array of same dimension as temp, containing salinity values (psu)
# press = array of same dimension as temp, containing pressure values (bar)
#         for which 0.1*depth (in m) is close enough
# Output: rho = array of same dimension as temp, containing density values 
#               (kg/m^3)
def unesco (temp, salt, press):

    # Set constants
    b0 = 8.24493e-1
    b1 = -4.0899e-3
    b2 = 7.6438e-5
    b3 = -8.2467e-7
    b4 = 5.3875e-9
    c0 = -5.72466e-3
    c1 = 1.0227e-4
    c2 = -1.6546e-6
    d0 = 4.8314e-4
    a0 = 999.842594
    a1 = 6.793952e-2
    a2 = -9.095290e-3
    a3 = 1.001685e-4
    a4 = -1.120083e-6
    a5 = 6.536332e-9
    f0 = 54.6746
    f1 = -0.603459
    f2 = 1.09987e-2
    f3 = -6.1670e-5
    g0 = 7.944e-2
    g1 = 1.6483e-2
    g2 = -5.3009e-4
    i0 = 2.2838e-3
    i1 = -1.0981e-5
    i2 = -1.6078e-6
    j0 = 1.91075e-4
    m0 = -9.9348e-7
    m1 = 2.0816e-8
    m2 = 9.1697e-10
    e0 = 19652.21
    e1 = 148.4206
    e2 = -2.327105
    e3 = 1.360477e-2
    e4 = -5.155288e-5
    h0 = 3.239908
    h1 = 1.43713e-3
    h2 =  1.16092e-4
    h3 = -5.77905e-7
    k0 = 8.50935e-5
    k1 = -6.12293e-6
    k2 = 5.2787e-8

    rho_0 = a0 + a1*temp + a2*temp**2 + a3*temp**3 + a4*temp**4 + a5*temp**5 + \
            (b0 + b1*temp + b2*temp**2 + b3*temp**3 + b4*temp**4)*salt + \
            (c0 + c1*temp + c2*temp**2)*salt**1.5 + d0*salt**2
    A = h0 + h1*temp + h2*temp**2 + h3*temp**3 + \
        (i0 + i1*temp + i2*temp**2)*salt + j0*salt**1.5
    B = k0 + k1*temp + k2*temp**2 + (m0 + m1*temp + m2*temp**2)*salt
    K = e0 + e1*temp + e2*temp**2 + e3*temp**3 + e4*temp**4 + \
        (f0 + f1*temp + f2*temp**2 + f3*temp**3)*salt + \
        (g0 + g1*temp + g2*temp**2)*salt**1.5 + A*press + B*press**2
    rho = rho_0/(1 - press/K)

    return rho

    
