from cmip5_paths import *
from numpy import *
from netCDF4 import Dataset
from matplotlib.pyplot import *

# Calculate root-mean-square errors (as in Gleckler et al., 2008) for each of
# 39 CMIP5 models and the multi-model mean, with respect to 11 ERA-Interim
# atmospheric variables. The domain is the Southern Ocean (all longitudes, and
# latitudes from the northern boundary of ROMS to the southernmost ocean point
# not in an ice shelf cavity) and the monthly climatology averaged over
# 1992-2005 inclusive. Also calculate the relative errors as in Gleckler et al.
# and make a "portrait plot" of coloured tiles in a model vs variable matrix.
# Save both rms errors and relative errors into text files.
# Note that to run this script, you must have previously run 
# eraint_climatology_netcdf.py, cmip5_atmos_climatology_netcdf.py, and
# mmm_atmos_netcdf.py.
def cmip5_eraint_rms_errors():

    # Path to directory containing climatology NetCDF files (created by
    # eraint_climatology_netcdf.py and cmip5_atmos_climatology_netcdf.py)
    directory = '/short/y99/kaa561/CMIP5_forcing/atmos/'
    # Variable names in NetCDF files
    var_names = ['Pair', 'Tair', 'Hair', 'cloud', 'Uwind', 'Vwind', 'precip', 'snow', 'evap', 'swrad', 'lwrad']
    # Path to ROMS grid
    roms_grid = '/short/m68/kaa561/ROMS-CICE-MCT/apps/common/grid/circ30S_quarterdegree_10m.nc'
    # Radius of the Earth in metres
    r = 6.371e6
    # Degrees to radians conversion factor
    deg2rad = pi/180.0

    # Read ROMS grid
    id = Dataset(roms_grid, 'r')
    lat_rho = id.variables['lat_rho'][:,:]
    mask_rho = id.variables['mask_rho'][:,:]
    mask_zice = id.variables['mask_zice'][:,:]
    id.close()
    # Determine latitude bounds of ocean points not in ice shelf cavities
    lat_roms = ma.masked_where(mask_rho-mask_zice==0, lat_rho)
    latS = amin(lat_roms)
    latN = amax(lat_roms)

    # Read ERA-Interim grid
    id = Dataset(directory + 'ERA-Interim.nc', 'r')
    lon = id.variables['longitude'][:]
    lat = id.variables['latitude'][:]
    id.close()

    # Find bounds on j to extend just past latS and latN
    j_min = nonzero(lat < latN)[0][0] - 1
    j_max = nonzero(lat < latS)[0][0] + 1
    lat = lat[j_min:j_max]

    # Calculate dx, dy, dt
    # ERA-Interim has constant spacing for both lat and lon; find these values
    # (should be 0.75 degrees)
    dlon = mean(lon[1:] - lon[0:-1])
    dlat = mean(lat[1:] - lat[0:-1])
    # dx = r*cos(lat)*dlon where lat and dlon are converted to radians
    dx = abs(r*cos(lat*deg2rad)*dlon*deg2rad)
    # dy = r*dlat where dlat is converted to radians
    dy = abs(r*dlat*deg2rad)
    # dt = number of seconds in each month
    dt = array([31.0, 28.25, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0, 31.0])*24*60*60
    # Copy dx and dt to 3D arrays with dimension time x latitude x longitude
    # We don't have to do this for dy because it's a scalar
    dx = transpose(tile(dx, (size(lon), 1)))
    dx = tile(dx, (12, 1, 1))
    dt = tile(dt, (size(lon), size(lat), 1))
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

        # Read ERA-Interim data, trimmed to latitude bounds
        print 'Processing ERA-Interim'
        id = Dataset(directory + 'ERA-Interim.nc', 'r')
        era_data = id.variables[var][:,j_min:j_max,:]
        id.close()

        # Loop over models
        for j in range(len(model_names)):
            model_name = model_names[j]
            print 'Processing ' + model_name

            # Read model data, trimmed to latitude bounds
            # This has already been interpolated to the ERA-Interim grid
            id = Dataset(directory + model_name + '.nc', 'r')
            model_data = id.variables[var][:,j_min:j_max,:]
            id.close()

            # Figure out if this variable is missing
            try:
                mask = model_data.mask
            except(AttributeError):
                # There is no mask at all
                mask = False
            if all(mask):
                # Everything is masked, i.e. the variable is missing
                # Mask out the RMS error too
                rms_errors[i,j] = ma.masked
            else:
                # Calculate RMS error: square root of the sum of the squares
                # of the difference between this model and ERA-Interim at every
                # point, weighted by dx*dy*dt
                rms_errors[i,j] = sqrt(sum((model_data - era_data)**2*dx*dy*dt)/sum(dx*dy*dt))

        # Calculate relative error as (E-E')/E where E is the RMS error for each
        # model, and E' is the median error across all models for this variable
        median_error = median(rms_errors[i,:-1])
        for j in range(len(models)+1):
            relative_errors[i,j] = (rms_errors[i,j] - median_error)/rms_errors[i,j]

    # Make the portrait plot
    fig = figure(figsize=(8,20))
    ax = fig.add_subplot(111)
    img = ax.pcolormesh(transpose(relative_errors), vmin=-0.8, vmax=0.8, cmap='RdBu_r')

    # Configure plot
    # Add model names and variable names to the axes
    xticks(arange(0.5, len(var_names)+0.5), var_names, rotation=-45)
    yticks(arange(0.5, len(model_names)+0.5), model_names)
    xlim([0, len(var_names)])
    ylim([len(model_names), 0])
    title('Relative error', fontsize=16)
    # Move the plot over a bit to make room for all the model names
    box = ax.get_position()
    ax.set_position([box.x0+box.width*0.1, box.y0, box.width*0.9, box.height])
    cbaxes = fig.add_axes([0.3,0.03,0.5,0.025])
    colorbar(img, orientation='horizontal', ticks=arange(-0.75, 1, 0.25), cax=cbaxes)
    savefig('relative_errors_eraint.png')

    # Write RMS errors to a text file, tabulated nicely
    f = open('rms_errors_eraint.txt', 'w')
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
    f = open('relative_errors_eraint.txt', 'w')
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

    cmip5_eraint_rms_errors()
