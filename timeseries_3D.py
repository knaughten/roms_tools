from netCDF4 import Dataset
from numpy import *
from matplotlib.pyplot import *
from os.path import *
from cartesian_grid_3d import *
from rotate_vector_roms import *

# Calculate and plot timeseries of total ocean heat content, average salinity, 
# and total kinetic energy during a ROMS simulation.
# Takes 32 GB of memory for Kaitlin's circumpolar quarter-degree grid to 30S.
# Input:
# grid_path = path to ROMS grid file
# file_path = path to ROMS history/averages file
# log_path = path to log file (if it exists, previously calculated values will
#            be read from it; regardless, it will be overwritten with all
#            calculated values following computation)
def timeseries_3D (grid_path, file_path, log_path):

    # Grid parameters
    theta_s = 0.9
    theta_b = 4.0
    hc = 40
    N = 31
    rho0 = 1000.0    # Reference density (kg/m^3)
    Cp = 3974        # Specific heat of polar seawater (J/K/kg)
    C2K = 273.15     # Celsius to Kelvin conversion

    time = []
    ohc = []
    avgsalt = []
    tke = []
    # Check if the log file exists
    if exists(log_path):
        print 'Reading previously calculated values'
        f = open(log_path, 'r')
        # Skip first line (header for time array)
        f.readline()
        for line in f:
            try:
                time.append(float(line))
            except(ValueError):
                # Reached the header for the next variable
                break
        for line in f:
            try:
                ohc.append(float(line))
            except(ValueError):
                break
        for line in f:
            try:
                avgsalt.append(float(line))
            except(ValueError):
                break
        for line in f:
            tke.append(float(line))
        f.close()

    print 'Analysing grid'
    id = Dataset(grid_path, 'r')
    h = id.variables['h'][:-15,:-3]
    zice = id.variables['zice'][:-15,:-3]    
    lon = id.variables['lon_rho'][:-15,:-3]
    lat = id.variables['lat_rho'][:-15,:-3]
    mask = id.variables['mask_rho'][:-15,:-3]
    # Keep the overlapping periodic boundary on "angle" for now
    angle = id.variables['angle'][:-15,:]
    id.close()

    id = Dataset(file_path, 'r')
    # Read time values and convert from seconds to years
    new_time = id.variables['ocean_time'][:]/(60*60*24*365.25)
    num_time = size(new_time)
    # Concatenate with time values from log file
    for t in range(num_time):        
        time.append(new_time[t])

    # Process 10 time indices at a time so we don't use too much memory
    start_t = 0
    while True:
        end_t = min(start_t+10, num_time)
        print 'Processing time indices ' + str(start_t+1) + ' to ' + str(end_t)
        num_time_curr = end_t-start_t

        print 'Calculating time-dependent dV'
        # Read time-dependent sea surface height
        zeta = id.variables['zeta'][start_t:end_t,:-15,:-3]
        # Calculate time-dependent dz
        dz = ma.empty([num_time_curr, N, size(lon,0), size(lon,1)])
        for t in range(num_time_curr):
            # dx and dy will be recomputed unnecessarily each timestep
            # but that's ok
            dx, dy, dz_tmp, z = cartesian_grid_3d(lon, lat, h, zice, theta_s, theta_b, hc, N, zeta[t,:,:])
            dz[t,:,:,:] = dz_tmp
        # Calculate time-dependent dV and mask with land mask
        # Here mask, dx, dy are all copied into arrays of dimension
        # time x depth x lat x lon
        dV = ma.masked_where(tile(mask, (num_time_curr,N,1,1))==0, tile(dx, (num_time_curr,1,1,1))*tile(dy, (num_time_curr,1,1,1))*dz)

        print 'Reading data'
        temp = id.variables['temp'][start_t:end_t,:,:-15,:-3]
        salt = id.variables['salt'][start_t:end_t,:,:-15,:-3]
        rho = id.variables['rho'][start_t:end_t,:,:-15,:-3] + rho0
        # Keep overlapping periodic boundary for u and v
        u_xy = id.variables['u'][start_t:end_t,:,:-15,:]
        v_xy = id.variables['v'][start_t:end_t,:,:-15,:]

        print 'Interpolating velocities onto rho-grid'
        # We are actually rotating them at the same time as interpolating
        # which is a bit of unnecessary work (sum of squares won't change with
        # rotation) but not much extra work, and it's conveneint
        u = ma.empty(shape(temp))
        v = ma.empty(shape(temp))
        for t in range(num_time_curr):
            for k in range(N):
                u_tmp, v_tmp = rotate_vector_roms(u_xy[t,k,:,:], v_xy[t,k,:,:], angle)
                u[t,k,:,:] = u_tmp[:,:-3]
                v[t,k,:,:] = v_tmp[:,:-3]

        print 'Building timeseries'
        for t in range(num_time_curr):
            # Integrate temp*rho*Cp*dV to get OHC
            ohc.append(sum((temp[t,:,:,:]+C2K)*rho[t,:,:,:]*Cp*dV[t,:,:,:]))
            # Average salinity (weighted with rho*dV)
            avgsalt.append(sum(salt[t,:,:,:]*rho[t,:,:,:]*dV[t,:,:,:])/sum(rho[t,:,:,:]*dV[t,:,:,:]))
            # Integrate 0.5*rho*speed^2*dV to get TKE
            tke.append(sum(0.5*rho[t,:,:,:]*(u[t,:,:,:]**2 + v[t,:,:,:]**2)*dV[t,:,:,:]))

        # Get ready for next 10 time indices
        if end_t == num_time:
            break
        start_t = end_t

    id.close()

    print 'Plotting ocean heat content'
    clf()
    plot(time, ohc)
    xlabel('Years')
    ylabel('Southern Ocean Heat Content (J)')
    grid(True)
    savefig('ohc.png')

    print 'Plotting average salinity'
    clf()
    plot(time, avgsalt)
    xlabel('Years')
    ylabel('Southern Ocean Average Salinity (psu)')
    grid(True)
    savefig('avgsalt.png')

    print 'Plotting total kinetic energy'
    clf()
    plot(time, tke)
    xlabel('Years')
    ylabel('Southern Ocean Total Kinetic Energy (J)')
    grid(True)
    savefig('tke.png')

    print 'Saving results to log file'
    f = open(log_path, 'w')
    f.write('Time (years):\n')
    for elm in time:
        f.write(str(elm) + '\n')
    f.write('Southern Ocean Heat Content (J):\n')
    for elm in ohc:
        f.write(str(elm) + '\n')
    f.write('Southern Ocean Average Salinity (psu):\n')
    for elm in avgsalt:
        f.write(str(elm) + '\n')
    f.write('Southern Ocean Total Kinetic Energy (J):\n')
    for elm in tke:
        f.write(str(elm) + '\n')
    f.close()


# Command-line interface
if __name__ == "__main__":

    grid_path = raw_input("Path to ROMS grid file: ")
    file_path = raw_input("Path to ROMS history/averages file: ")
    log_path = raw_input("Path to logfile to save values and/or read previously calculated values: ")
    timeseries_3D(grid_path, file_path, log_path)

    

    
