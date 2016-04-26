from cmip5_paths import *
from numpy import *
from netCDF4 import Dataset
from matplotlib.pyplot import *

# Calculate root-mean square errors (adapted from Gleckler et al., 2008 to be
# in depth-longitude space instead of latitude-longitude space) for each of 39
# CMIP5 models and the multi-model mean, with respect to 4 ECCO2 ocean
# variables at the northern boundary of ROMS (currently 30S). Use the monthly
# climatology averaged over 1992-2005 inclusive. Also calculate the relative
# errors as in Gleckler et al. and make a "portrait plot" of coloured tiles
# in a model vs variable matrix. Save both rms errors and relative errors into
# text files. Note that to run this script, you must have previously run
# ecco2_climatology_netcdf.py, cmip5_ocean_climatology_netcdf.py, and
# mmm_ocean_netcdf.py.
def cmip5_ecco2_rms_errors ():

    # Path to directory containing climatolgoy NetCDF files (created by
    # ecco2_climatology_netcdf.py and cmip5_ocean_climatology_netcdf.py)
    directory = '/short/y99/kaa561/CMIP5_forcing/ocean/'
    # Variable names in NetCDF files
    var_names = ['temp', 'salt', 'u', 'v']
    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    # Latitude of the northern boundary
    nbdry = -30.0

    # Read ECCO2 grid
    id = Dataset(directory + 'ECCO2.nc', 'r')
    lon = id.variables['longitude'][:]
    depth = id.variables['depth'][:]
    id.close()

    # Calculate dx, dz, dt
    # ECCO2 has constant spacing for longitude; find this value (should be 0.25
    # degrees)
    dlon = mean(lon[1:] - lon[0:-1])
    # dx = r*cos(lat)*dlon where lat and dlon are converted to radians
    dx = abs(r*cos(nbdry*deg2rad)*dlon*deg2rad)
    # We have depth values in the middle of each cell; interpolate to the edges
    z_edges = zeros(size(depth)+1)
    z_edges[1:-1] = 0.5*(depth[0:-1] + depth[1:])
    # Assume surface is 0
    z_edges[0] = 0.0
    # Extrapolate to bottom
    z_edges[-1] = 2*depth[-1] - z_edges[-2]
    # Now calculate dz
    dz = z_edges[1:] - z_edges[0:-1]
    # dt = number of seconds in each month
    dt = array([31.0, 28.25, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0, 31.0])*24*60*60
    # Copy dz and dt to 3D arrays with dimension time x depth x longitude
    # We don't have to do this for dx because it's a scalar
    dz = transpose(tile(dz, (size(lon), 1)))
    dz = tile(dz, (12, 1, 1))
    dt = tile(dt, (size(lon), size(depth), 1))
    dt = transpose(dt)    

    # Build a list of CMIP5 Model objects
    models = build_model_list()
    # Build a corresponding list of model names
    model_names = []
    for model in models:
        model_names.append(model.name)
    # Add the multi-model mean to the list
    model_names.append('MMM')

    # Set up arrays for RMS errors and relative errors, of dimension
    # variable x model
    rms_errors = ma.empty([len(var_names), len(model_names)])
    relative_errors = ma.empty(shape(rms_errors))

    # Loop over variables
    for i in range(len(var_names)):
        var = var_names[i]
        print 'Variable ' + var

        # Read ECCO2 data
        print 'Processing ECCO2'
        id = Dataset(directory + 'ECCO2.nc', 'r')
        ecco_data = id.variables[var][:,:,:]
        id.close()

        if i == 0:
            # Mask dz with the land mask on the first iteration
            # If at least one of the integrands are masked, land points
            # will be automatically excluded from the integrals
            dz_mask = ma.masked_where(ecco_data.mask, dz)

        # Loop over models
        for j in range(len(model_names)):
            model_name = model_names[j]
            print 'Processing ' + model_name

            # Read model data
            # This has already been interpolated to the ECCO2 grid at the
            # northern boundary of ROMS
            id = Dataset(directory + model_name + '.nc', 'r')
            model_data = id.variables[var][:,:,:]
            id.close()

            if all(model_data.mask):
                # This variable is missing; mask out the RMS error
                rms_errors[i,j] = ma.masked
            else:
                # Calculate RMS error: square root of the sum of the squares of
                # the difference between this model and ECCO2 at every point,
                # weighted by dx*dz*dt, with land points excluded since dz
                # has been masked
                rms_errors[i,j] = sqrt(sum((model_data - ecco_data)**2*dx*dz_mask*dt)/sum(dx*dz_mask*dt))

        # Calculate relative error as (E-E')/E where E is the RMS error for each
        # model, and E' is the median error across all models for this variable
        median_error = median(rms_errors[i,:-1])
        for j in range(len(models)+1):
            relative_errors[i,j] = (rms_errors[i,j] - median_error)/rms_errors[i,j]

    # Make the portrait plot
    fig = figure(figsize=(8,20))
    ax = fig.add_subplot(111)
    img = ax.pcolormesh(transpose(relative_errors), vmin=-1.5, vmax=1.5, cmap='RdBu_r')

    # Configure plot
    # Add model names and variable names to the axes
    xticks(arange(0.5, len(var_names)+0.5), var_names)
    yticks(arange(0.5, len(model_names)+0.5), model_names)
    xlim([0, len(var_names)])
    ylim([len(model_names), 0])
    title('Relative error', fontsize=16)
    # Move the plot over a bit to make room for all the model names
    box = ax.get_position()
    ax.set_position([box.x0+box.width*0.1, box.y0, box.width*0.9, box.height])
    cbaxes = fig.add_axes([0.3,0.03,0.5,0.025])
    colorbar(img, orientation='horizontal', ticks=arange(-1.5, 1.75, 0.5), cax=cbaxes)
    savefig('relative_errors_ecco2.png')

    # Write RMS errors to a text file, tabulated nicely
    f = open('rms_errors_ecco2.txt', 'w')
    f.write('{0:16}'.format('Model'))
    for var in var_names:
        f.write('{0:16}'.format(var))
    f.write('\n')
    for j in range(len(model_names)):
        f.write('{0:16}'.format(model_names[j]))
        for i in range(len(var_names)):
            f.write('{0:16}'.format(str(round(rms_errors[i,j],4))))
        f.write('\n')
    f.close()

    # Write relative errors to a text file, tabulated nicely
    f = open('relative_errors_ecco2.txt', 'w')
    f.write('{0:16}'.format('Model'))
    for var in var_names:
        f.write('{0:16}'.format(var))
    f.write('\n')
    for j in range(len(model_names)):
        f.write('{0:16}'.format(model_names[j]))
        for i in range(len(var_names)):
            f.write('{0:16}'.format(str(round(relative_errors[i,j],4))))
        f.write('\n')
    f.close()


# Command-line interface
if __name__ == '__main__':

    cmip5_ecco2_rms_errors()
        
        

        
