from numpy import *
from netCDF4 import Dataset
from matplotlib.pyplot import *

# Plot precipitation (1992-2005 average) in the ERA-Interim versus GPCP datasets
# over the Southern Ocean.
def era_gpcp ():

    # Beginning and end of paths to ERA-Interim precip files, monthly averaged,
    # on original grid
    era_head = '/short/y99/kaa561/CMIP5_forcing/atmos/climatology/ERA_Interim_monthly/FC_'
    era_tail = '_monthly_orig.nc'
    # Beginning of paths to GPCP files
    gpcp_head = '/short/m68/kaa561/gpcp/gpcp_cdr_v23rB1_y'
    # Days per month (leap years will be handled later)
    days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    start_year = 1992
    end_year = 2005
    deg2rad = pi/180.0

    # Read each grid
    id = Dataset(gpcp_head + str(start_year) + '_m01.nc', 'r')
    gpcp_lat_1d = id.variables['latitude'][:]
    gpcp_lon_1d = id.variables['longitude'][:]
    id.close()
    id = Dataset(era_head + str(start_year) + era_tail, 'r')
    era_lat_1d = id.variables['latitude'][:]
    era_lon_1d = id.variables['longitude'][:]
    id.close()

    # Start integrating precipitation data
    gpcp_precip = zeros([size(gpcp_lat_1d), size(gpcp_lon_1d)])
    era_precip = zeros([size(era_lat_1d), size(era_lon_1d)])
    # Number of days so far
    ndays = 0
    # Loop over years
    for year in range(start_year, end_year+1):
        # Check if it's a leap year
        leap_year = False
        if mod(year, 4) == 0:
            leap_year = True
            if mod(year, 100) == 0:
                leap_year = False
                if mod(year, 400) == 0:
                    leap_year = True
        if leap_year:
            days_per_month[1] = 29
        else:
            days_per_month[1] = 28
        # Loop over months
        for month in range(12):
            if month + 1 < 10:
                month_string = '0' + str(month+1)
            else:
                month_string = str(month+1)
            # Read GPCP data
            id = Dataset(gpcp_head + str(year) + '_m' + month_string + '.nc', 'r')
            gpcp_precip += id.variables['precip'][:,:]*days_per_month[month]
            id.close()
            # Read ERA-Interim data
            id = Dataset(era_head + str(year) + era_tail, 'r')
            era_precip += id.variables['tp'][month,:,:]*days_per_month[month]
            id.close()
            ndays += days_per_month[month]
    # Convert from integral to average
    gpcp_precip /= ndays
    era_precip /= ndays
    # Convert to m/y
    gpcp_precip *= 1e-3*365.25  # originally mm/day
    era_precip *= 2*365.25  # originally m/12h

    # Make 2D versions of lat and lon arrays
    gpcp_lon, gpcp_lat = meshgrid(gpcp_lon_1d, gpcp_lat_1d)
    era_lon, era_lat = meshgrid(era_lon_1d, era_lat_1d)

    # Polar coordinates for plotting
    gpcp_x = -(gpcp_lat+90)*cos(gpcp_lon*deg2rad+pi/2)
    gpcp_y = (gpcp_lat+90)*sin(gpcp_lon*deg2rad+pi/2)
    era_x = -(era_lat+90)*cos(era_lon*deg2rad+pi/2)
    era_y = (era_lat+90)*sin(era_lon*deg2rad+pi/2)

    # Colour levels
    lev = linspace(0, 2, num=50)
    # Northern boundary 40S
    bdry = -40+90

    # Make figure
    fig = figure(figsize=(20,9))
    # GPCP
    fig.add_subplot(1,2,1, aspect='equal')
    contourf(gpcp_x, gpcp_y, gpcp_precip, lev, extend='max')
    title('GPCP', fontsize=24)
    xlim([-bdry, bdry])
    ylim([-bdry, bdry])
    axis('off')
    # ERA-Interim
    fig.add_subplot(1,2,2, aspect='equal')
    img = contourf(era_x, era_y, era_precip, lev, extend='max')
    title('ERA-Interim', fontsize=24)
    xlim([-bdry, bdry])
    ylim([-bdry, bdry])
    axis('off')
    # Add colourbar on the bottom
    cbaxes = fig.add_axes([0.3, 0.04, 0.4, 0.04])
    cbar = colorbar(img, orientation='horizontal', cax=cbaxes)
    # Main title
    suptitle('Precipitation (m/y), 1992-2005 average', fontsize=30)

    fig.savefig('era_gpcp.png')


# Command-line interface
if __name__ == "__main__":

    era_gpcp()
    
