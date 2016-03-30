from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from calc_z import *

# Analyse a ROMS spinup by calculating and plotting 8 timeseries:
# Total heat content
# Total salt content
# Area-averaged ice shelf melt rate
# Ice shelf basal mass loss
# Total kinetic energy
# Maximum velocity
# Drake Passage transport
# Total sea ice extent


# Given the path to a ROMS grid file, calculate differentials for later
# integration.
# Input: grid_path = string containing path to ROMS grid file
# Output:
# dA = differential of area on the 2D rho-grid, masked with zice
# dV = differential of volume on the 3D rho-grid (depth x lat x lon), masked
#      with land mask
# dydz = differential of area in the y-z direction for each cell on the 3D
#        rho-grid, masked with land mask
def calc_grid (grid_path):

    # Grid parameters
    theta_s = 0.9
    theta_b = 4.0
    hc = 40
    N = 31
    # Radius of the Earth in m
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Northern boundary of ROMS grid
    nbdry_val = -38

    # Read grid variables
    id = Dataset(grid_path, 'r')
    h = id.variables['h'][:,:]
    zice = id.variables['zice'][:,:]
    lon = id.variables['lon_rho'][:,:]
    lat = id.variables['lat_rho'][:,:]
    id.close()
    # Save dimensions
    num_lat = size(lon, 0)
    num_lon = size(lon, 1)

    # Add or subtract 360 from longitude values which wrap around
    # so that longitude increases monotonically from west to east
    i = tile(arange(1, num_lon+1), (num_lat, 1))
    index1 = nonzero((i > 1200)*(lon < 100))
    lon[index1] = lon[index1] + 360
    index2 = nonzero((i < 200)*(lon > 300))
    lon[index2] = lon[index2] - 360

    # Interpolate to get longitude at the edges of each cell
    w_bdry = 0.5*(lon[:,0] + lon[:,-1] - 360)
    middle_lon = 0.5*(lon[:,0:-1] + lon[:,1:])
    e_bdry = 0.5*(lon[:,0] + 360 + lon[:,-1])
    lon_edges = ma.concatenate((w_bdry[:,None], middle_lon, e_bdry[:,None]), axis=1)
    # Subtract to get the change in longitude over each cell
    dlon = abs(lon_edges[:,1:] - lon_edges[:,0:-1])

    # Similarly for latitude
    s_bdry = lat[0,:]
    middle_lat = 0.5*(lat[0:-1,:] + lat[1:,:])
    n_bdry = lat[-1,:]*0 + nbdry_val
    lat_edges = ma.concatenate((s_bdry[None,:], middle_lat, n_bdry[None,:]))
    dlat = lat_edges[1:,:] - lat_edges[0:-1,:]

    # Convert from spherical to Cartesian coordinates
    # dy = r*dlat where dlat is converted to radians
    dy_2d = r*dlat*pi/180.0    
    # dx = r*cos(lat)*dlon where lat and dlon are converted to radians
    dx_2d = r*cos(pi*lat/180.0)*dlon*pi/180.0

    # Calculate dA and mask with zice
    dA = dx_2d*dy_2d
    id = Dataset(grid_path, 'r')
    mask = id.variables['mask_zice'][:,:]
    id.close()
    dA = ma.masked_where(mask==0, dA)

    # Copy dx and dy into 3D arrays, same at each depth level
    dy = tile(dy_2d, (N,1,1))
    dx = tile(dx_2d, (N,1,1))

    # Get a 3D array of z-coordinates; sc_r and Cs_r are unused in this script
    z, sc_r, Cs_r = calc_z(h, zice, lon, lat, theta_s, theta_b, hc, N)
    # We have z at the midpoint of each cell, now find it on the top and
    # bottom edges of each cell
    z_edges = zeros((size(z,0)+1, size(z,1), size(z,2)))
    z_edges[1:-1,:,:] = 0.5*(z[0:-1,:,:] + z[1:,:,:])
    # At surface, z = zice; at bottom, extrapolate
    z_edges[-1,:,:] = zice[:,:]
    z_edges[0,:,:] = 2*z[0,:,:] - z_edges[1,:,:]
    # Now find dz
    dz = z_edges[1:,:,:] - z_edges[0:-1,:,:]

    # Read land mask
    id = Dataset(grid_path, 'r')
    mask = id.variables['mask_rho'][:,:]
    id.close()
    mask = tile(mask, (N,1,1))
    # Calculate dV and mask with land mask
    dV = ma.masked_where(mask==0, dx*dy*dz)
    # Similarly for dydz
    dydz = ma.masked_where(mask==0, dy*dz)

    return dA, dV, dydz


# Read and return density.
# Input:
# file_path = path to ocean history/averages file
# t = timestep index in file_path
# Output: rho = density field at timestep t
def get_rho (file_path, t):

    # Reference density
    # New users: edit this based on the value in your .in ROMS configuration file
    rho0 = 1025.0

    id = Dataset(file_path, 'r')
    # Read density anomalies, add rho0 to get absolute density
    # Convert to float128 to prevent overflow later
    rho = ma.asarray(id.variables['rho'][t,:,:,:], dtype=float128) + rho0
    id.close()
    return rho


# Calculate ocean heat content at the given timestep t.
# Input:
# file_path = path to ocean history/averages file
# dV = elements of volume on the rho grid, masked with land mask
# rho = density on the rho grid at timestep t
# t = timestep index in file_path
# Output: ohc = ocean heat content (J)
def calc_ohc (file_path, dV, rho, t):

    # Specific heat of polar seawater (J/K/kg)
    cp = 3974.0
    # Celsius to Kelvin conversion constant
    celsius2kelvin = 273.15

    # Read temperature, converting to float128 to prevent overflow during
    # integration
    id = Dataset(file_path, 'r')
    temp = ma.asarray(id.variables['temp'][t,:,:,:], dtype=float128)
    # Convert from Celsius to Kelvin
    temp = temp + celsius2kelvin
    id.close()

    # Integrate temp*rho*cp over volume to get OHC
    ohc = sum(temp*rho*cp*dV)
    return ohc


# Calculate total salt content at the given timestep t.
# Input:
# file_path = path to ocean history/averages file
# dV = elements of volume on the rho grid, masked with land mask
# rho = density on the rho grid at timestep t
# t = timestep index in file_path
# Output: totalsalt = total salt content (kg)
def calc_totalsalt (file_path, dV, rho, t):

    # Read salinity, converting to float128 to prevent overflow during
    # integration
    id = Dataset(file_path, 'r')
    salt = ma.asarray(id.variables['salt'][t,:,:,:], dtype=float128)
    id.close()

    # Integrate 1e-3*salt*rho over volume to get total mass of salt
    totalsalt = sum(1e-3*salt*rho*dV)
    return totalsalt    


# Calculate area-averaged ice shelf melt rate at the given timestep t.
# Input:
# file_path = path to ocean history/averages file
# dA = elements of area on the rho grid, masked with zice
# t = timestep index in file_path
# Output: avgismr = area-averaged ice shelf melt rate (m/y)
#         ismr = 2D ice shelf melt rate field (m/y) at this timestep
def calc_avgismr (file_path, dA, t):

    # Read ice shelf melt rate, converting to float128 to prevent overflow
    # during integration
    id = Dataset(file_path, 'r')
    ismr = ma.asarray(id.variables['m'][t,:,:], dtype=float128)
    # Convert from m/s to m/y
    ismr = ismr*365.25*24*60*60
    id.close()    

    # Integrate ismr over area and divide by total area to get average
    avgismr = sum(ismr*dA)/sum(dA)
    return avgismr, ismr


# Calculate net basal mass loss based on the given ice shelf melt rate field.
# Input:
# ismr = 2D ice shelf melt rate field (m/y)
# dA = elements of area on the rho grid, masked with zice
# Output: massloss = net basal mass loss (Gt/y)
def calc_massloss (ismr, dA):

    # Density of ice in kg/m^3
    rho_ice = 916

    # Integrate over area to get volume loss
    volumeloss = sum(ismr*dA)
    # Convert to mass loss in Gt/y
    massloss = 1e-12*rho_ice*volumeloss
    return massloss


# Calculate total kinetic energy at the given timestep t.
# Input:
# file_path = path to ocean history/averages file
# dV = elements of volume on the rho grid, masked with land mask
# rho = density on the rho grid at timestep t
# t = timestep index in file_path
# Output: tke = total kinetic energy (J)
def calc_tke (file_path, dV, rho, t):

    # Read u and v, converting to float 128 to prevent overflow during
    # integration
    id = Dataset(file_path, 'r')
    u = ma.asarray(id.variables['u'][t,:,:,:], dtype=float128)
    v = ma.asarray(id.variables['v'][t,:,:,:], dtype=float128)
    id.close()

    # Interpolate u onto the rho-grid
    w_bdry_u = 0.5*(u[:,:,0] + u[:,:,-1])
    middle_u = 0.5*(u[:,:,0:-1] + u[:,:,1:])
    e_bdry_u = w_bdry_u[:,:]
    u_rho = ma.concatenate((w_bdry_u[:,:,None], middle_u, e_bdry_u[:,:,None]), axis=2)

    # Interpolate v onto the rho-grid
    s_bdry_v = v[:,0,:]
    middle_v = 0.5*(v[:,0:-1,:] + v[:,1:,:])
    n_bdry_v = v[:,-1,:]
    v_rho = ma.concatenate((s_bdry_v[:,None,:], middle_v, n_bdry_v[:,None,:]), axis=1)

    # Integrate 0.5*rho*(u^2 + v^2) over volume to get TKE
    tke = sum(0.5*rho*(u_rho**2 + v_rho**2)*dV)
    return tke, u_rho, v_rho


# Calculate the maximum velocity.
# Input: u_rho, v_rho = u and v at timestep t, interpolated to the rho-grid
# Output: maxvel = maximum velocity (m/s)
def calc_maxvel (u_rho, v_rho):

    return amax(sqrt(u_rho**2 + v_rho**2))


# Calculate zonal transport through the Drake Passage.
# Input:
# file_path = path to ocean history/averages file
# dydz = elements of area in the y-z direction for each cell in the 3D
#        rho-grid, masked with land mask
# u_rho = u at timestep t, interpolated to the rho-grid
# Output: drakepsg_trans = zonal transport through the Drake Passage (60W),
#                          integrated over depth and latitude
def calc_drakepsgtrans (file_path, dydz, u_rho):

    # Bounds on Drake Passage; edit for new grids
    # i-index of single north-south line to plot (representing a zonal slice);
    # it doesn't really matter which slice of the Drake Passage this is, due
    # to volume conservation
    i_DP = 1175    
    # j-indices of the southern tip of South America (j_min) and the northern
    # tip of the Antarctic Peninsula (j_max); make sure these are far enough
    # north/south to be land points, but not so far that they pass through the
    # land and become ocean again (eg Weddell Sea)
    j_min = 220
    j_max = 300

    # Trim arrays to these bounds
    u_rho_DP = u_rho[:,j_min:j_max,i_DP]
    dydz_DP = dydz[:,j_min:j_max,i_DP]

    # Calculate transport
    transport = sum(u_rho_DP*dydz_DP)

    # Divide by 1e6 to convert to Sv
    return transport*1e-6


# Calculate total sea ice extent at the given timestep t.
# Input:
# cice_path = path to CICE history file
# dA = elements of area on the 2D rho grid (any mask will be removed)
# t = timestep index in file_path
# Output: totalice = total sea ice extent (m^2)
def calc_totalice (cice_path, dA, t):

    id = Dataset(cice_path, 'r')
    # Read sea ice area fraction at each grid cell
    aice = ma.asarray(id.variables['aice'][t,:,:], dtype=float128)
    id.close()

    # Remove masks on aice and dA, and fill aice with zeros on land mask
    # (numpy was throwing weird masking errors originally, and it doesn't
    # matter if dA is unmasked because we are integrating not averaging)
    aice_nomask = aice.data
    aice_nomask[aice.mask] = 0.0
    dA_nomask = dA.data

    # Find the cells with at least 15% sea ice
    extent_flag = aice_nomask >= 0.15

    # Integrate area of these cells
    totalice = sum(dA_nomask*extent_flag)
    # Convert to million km^2 and return
    return totalice*1e-12   


# Main routine
# Input:
# grid_path = path to ROMS grid file
# file_path = path to ocean history/averages file
# cice_path = path to CICE history file
# log_path = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
def spinup_plots (grid_path, file_path, cice_path, log_path):

    time = []
    ohc = []
    totalsalt = []
    avgismr = []
    massloss = []
    tke = []
    maxvel = []
    drakepsgtrans = []
    totalice = []
    # Check if the log file exists
    if exists(log_path):
        print 'Reading previously calculated values'
        f = open(log_path, 'r')
        # Skip the first line (header for time array)
        f.readline()
        for line in f:
            try:
                time.append(float(line))
            except (ValueError):
                # Reached the header for the next variable
                break
        for line in f:
            try:
                ohc.append(float(line))
            except (ValueError):
                break
        for line in f:
            try:
                totalsalt.append(float(line))
            except (ValueError):
                break
        for line in f:
            try:
                avgismr.append(float(line))
            except (ValueError):
                break
        for line in f:
            try:
                massloss.append(float(line))
            except (ValueError):
                break
        for line in f:
            try:
                tke.append(float(line))
            except (ValueError):
                break
        for line in f:
            try:
                maxvel.append(float(line))
            except (ValueError):
                break
        for line in f:
            try:
                drakepsgtrans.append(float(line))
            except (ValueError):
                break
        for line in f:
            totalice.append(float(line))
        f.close()

    # Calculate differentials
    print 'Analysing grid'
    dA, dV, dydz = calc_grid(grid_path)
    # Read time data and convert from seconds to years
    id = Dataset(file_path, 'r')
    new_time = id.variables['ocean_time'][:]/(365*24*60*60)
    id.close()
    # Concatenate with time values from log file
    for t in range(size(new_time)):
        time.append(new_time[t])

    # Process each timestep separately to prevent memory overflow
    for t in range(size(new_time)):
        print 'Processing timestep '+str(t+1)+' of '+str(size(new_time))
        rho = get_rho(file_path, t)
        print 'Calculating ocean heat content'
        ohc.append(calc_ohc(file_path, dV, rho, t))
        print 'Calculating total salt content'
        totalsalt.append(calc_totalsalt(file_path, dV, rho, t))
        print 'Calculating average ice shelf melt rate'
        avgismr_tmp, ismr = calc_avgismr(file_path, dA, t)
        avgismr.append(avgismr_tmp)
        print 'Calculating basal mass loss'
        massloss.append(calc_massloss(ismr, dA))
        print 'Calculating total kinetic energy'
        tke_tmp, u_rho, v_rho = calc_tke(file_path, dV, rho, t)
        tke.append(tke_tmp)
        print 'Calculating maximum velocity'
        maxvel.append(calc_maxvel(u_rho, v_rho))
        print 'Calculating Drake Passage transport'
        drakepsgtrans.append(calc_drakepsgtrans(file_path, dydz, u_rho))
        print 'Calculating total sea ice extent'
        totalice.append(calc_totalice(cice_path, dA, t))

    # Plot each timeseries in sequence
    print 'Plotting ocean heat content'
    clf()
    plot(time, ohc)
    xlabel('Years')
    ylabel('Southern Ocean Heat Content (J)')
    savefig('ohc.png')
    print 'Plotting total salt content'
    clf()
    plot(time, totalsalt)
    xlabel('Years')
    ylabel('Southern Ocean Salt Content (kg)')
    savefig('totalsalt.png')
    print 'Plotting average ice shelf melt rate'
    clf()
    plot(time, avgismr)
    xlabel('Years')
    ylabel('Area-averaged Ice Shelf Melt Rate (m/y)')
    savefig('avgismr.png')
    print 'Plotting basal mass loss'
    clf()
    plot(time, massloss)
    xlabel('Years')
    ylabel('Ice Shelf Basal Mass Loss (Gt/y)')
    savefig('massloss.png')
    print 'Plotting total kinetic energy'
    clf()
    plot(time, tke)
    xlabel('Years')
    ylabel('Southern Ocean Total Kinetic Energy (J)')
    savefig('tke.png')
    print 'Plotting maximum velocity'
    clf()
    plot(time, maxvel)
    xlabel('Years')
    ylabel('Maximum Southern Ocean Velocity (m/s)')
    savefig('maxvel.png')
    print 'Plotting Drake Passage transport'
    clf()
    plot(time, drakepsgtrans)
    xlabel('Years')
    ylabel('Drake Passage Transport (Sv)')
    savefig('drakepsgtrans.png')
    print 'Plotting total sea ice extent'
    clf()
    plot(time, totalice)
    xlabel('Years')
    ylabel(r'Total Sea Ice Extent (million km$^2$)')
    savefig('totalice.png')

    print 'Saving results to log file'
    f = open(log_path, 'w')
    f.write('Time (years):\n')
    for elm in time:
        f.write(str(elm) + '\n')
    f.write('Southern Ocean Heat Content (J):\n')
    for elm in ohc:
        f.write(str(elm) + '\n')
    f.write('Southern Ocean Salt Content (kg):\n')
    for elm in totalsalt:
        f.write(str(elm) + '\n')
    f.write('Area-averaged Ice Shelf Melt Rate (m/y):\n')
    for elm in avgismr:
        f.write(str(elm) + '\n')
    f.write('Ice Shelf Basal Mass Loss (Gt/y):\n')
    for elm in massloss:
        f.write(str(elm) + '\n')
    f.write('Southern Ocean Total Kinetic Energy (J):\n')
    for elm in tke:
        f.write(str(elm) + '\n')
    f.write('Maximum Southern Ocean Velocity (m/s):\n')
    for elm in maxvel:
        f.write(str(elm) + '\n')
    f.write('Drake Passage Transport (Sv):\n')
    for elm in drakepsgtrans:
        f.write(str(elm) + '\n')
    f.write('Total Sea Ice Extent (million km^2):\n')
    for elm in totalice:
        f.write(str(elm) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    grid_path = raw_input('Enter path to grid file: ')
    file_path = raw_input('Enter path to ocean history/averages file: ')
    cice_path = raw_input('Enter path to CICE history file: ')
    log_path = raw_input('Enter path to log file to save values and/or read previously calculated values: ')

    spinup_plots(grid_path, file_path, cice_path, log_path)


