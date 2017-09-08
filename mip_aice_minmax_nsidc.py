from netCDF4 import Dataset, num2date
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from monthly_avg_cice import *
# Import FESOM scripts (have to modify path first)
import sys
sys.path.insert(0, '/short/y99/kaa561/fesomtools')
from patches import *
from monthly_avg import *

# Make a 4x2 plot showing February (top) and September (bottom) sea ice
# concentration (1992-2015 average) for NSIDC, MetROMS, and FESOM. The fourth
# column shows timeseries of February and September sea ice extent for NSIDC,
# MetROMS, and FESOM.
# Input:
# cice_file = path to CICE output file containing 5-day averages for the entire
#             simulation
# cice_log = path to CICE logfile from timeseries_seaice_extent.py
# fesom_mesh_path = path to FESOM mesh directory
# fesom_output_dir = path to FESOM output directory containing one ice.mean.nc
#                    file for each year (5-day averages)
# fesom_log = path to FESOM logfile from timeseries_seaice_extent.py in
#             fesomtools
def mip_aice_minmax_nsidc (cice_file, cice_log, fesom_mesh_path, fesom_output_dir, fesom_log):

    # Range of years to process (exclude 2016 because no NSIDC data)
    start_year = 1992
    end_year = 2015
    # Starting and ending days for each month
    # Ignore leap years, they will be dealt with later
    start_day = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    end_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    # Beginning of FESOM output filenames
    fesom_file_head = 'MK44005'
    # Paths to reconstruct NSIDC files for each month
    nsidc_head1 = '/short/m68/kaa561/nsidc_aice/seaice_conc_monthly_sh_f11_'
    nsidc_head2 = '/short/m68/kaa561/nsidc_aice/seaice_conc_monthly_sh_f13_'
    nsidc_head3 = '/short/m68/kaa561/nsidc_aice/seaice_conc_monthly_sh_f17_'
    nsidc_tail = '_v02r00.nc'
    # Path to NSIDC monthly extent timeseries for February and September
    nsidc_feb_ts_file = '/short/m68/kaa561/nsidc_aice/S_02_extent_v2.1.csv'
    nsidc_sep_ts_file = '/short/m68/kaa561/nsidc_aice/S_09_extent_v2.1.csv'
    # FESOM plotting parameters
    circumpolar = True
    mask_cavities = True
    # Degrees to radians conversion factor
    deg2rad = pi/180.0
    num_years = end_year - start_year + 1

    print 'Processing NSIDC'

    # First read the grid
    id = Dataset(nsidc_head1 + '199201' + nsidc_tail, 'r')
    nsidc_lon = id.variables['longitude'][:,:]
    nsidc_lat = id.variables['latitude'][:,:]
    id.close()
    # Read February and September concentration for each year
    nsidc_feb = ma.empty([num_years, size(nsidc_lon,0), size(nsidc_lat,1)])
    nsidc_sep = ma.empty([num_years, size(nsidc_lon,0), size(nsidc_lat,1)])
    # Loop over years
    for year in range(start_year, end_year+1):
        # Reconstruct file paths
        if year < 1996:
            feb_file = nsidc_head1 + str(year) + '02' + nsidc_tail
            sep_file = nsidc_head1 + str(year) + '09' + nsidc_tail
        elif year < 2008:
            feb_file = nsidc_head2 + str(year) + '02' + nsidc_tail
            sep_file = nsidc_head2 + str(year) + '09' + nsidc_tail
        else:
            feb_file = nsidc_head3 + str(year) + '02' + nsidc_tail
            sep_file = nsidc_head3 + str(year) + '09' + nsidc_tail
        # Read February data and mask
        id = Dataset(feb_file, 'r')
        feb_data_tmp = id.variables['seaice_conc_monthly_cdr'][0,:,:]
        # (This variable is masked but sea ice concentration isn't)
        nsidc_mask = id.variables['stdev_of_seaice_conc_monthly_cdr'][0,:,:]
        id.close()
        # Apply mask
        feb_data = ma.empty(shape(feb_data_tmp))
        feb_data[:,:] = 0.0
        feb_data[~nsidc_mask.mask] = feb_data_tmp[~nsidc_mask.mask]
        feb_data[nsidc_mask.mask] = ma.masked
        # Save result
        nsidc_feb[year-start_year,:,:] = feb_data[:,:]
        # Repeat for September
        id = Dataset(sep_file, 'r')
        sep_data_tmp = id.variables['seaice_conc_monthly_cdr'][0,:,:]
        nsidc_mask = id.variables['stdev_of_seaice_conc_monthly_cdr'][0,:,:]
        id.close()
        sep_data = ma.empty(shape(sep_data_tmp))
        sep_data[:,:] = 0.0
        sep_data[~nsidc_mask.mask] = sep_data_tmp[~nsidc_mask.mask]
        sep_data[nsidc_mask.mask] = ma.masked
        nsidc_sep[year-start_year,:,:] = sep_data[:,:]
    # Average over years
    nsidc_feb = mean(nsidc_feb, axis=0)
    nsidc_sep = mean(nsidc_sep, axis=0)
    # Make sure mask is still there
    nsidc_feb[nsidc_mask.mask] = ma.masked
    nsidc_sep[nsidc_mask.mask] = ma.masked
    # Polar coordinates for plotting
    nsidc_x = -(nsidc_lat+90)*cos(nsidc_lon*deg2rad+pi/2)
    nsidc_y = (nsidc_lat+90)*sin(nsidc_lon*deg2rad+pi/2)
    # Choose boundaries based on extent of NSIDC grid
    bdry1 = amax(nsidc_x[:,0])
    bdry2 = amin(nsidc_x[:,-1])
    bdry3 = amin(nsidc_y[:,0])
    bdry4 = amax(nsidc_y[:,-1])

    # Now read extent timeseries
    nsidc_feb_extent = []
    f = open(nsidc_feb_ts_file, 'r')
    f.readline()
    # Skip the years we don't care about
    for year in range(1979, start_year):
        f.readline()
    # Read the years we care about
    for year in range(start_year, end_year+1):
        tmp = f.readline().split(',')
        # Extract the extent (second last column)
        nsidc_feb_extent.append(float(tmp[-2]))
    f.close()
    # Repeat for September
    nsidc_sep_extent = []
    f = open(nsidc_sep_ts_file, 'r')
    f.readline()
    for year in range(1979, start_year):
        f.readline()
    for year in range(start_year, end_year+1):
        tmp = f.readline().split(',')
        nsidc_sep_extent.append(float(tmp[-2]))
    f.close()

    print 'Processing MetROMS'

    # First read the grid
    id = Dataset(cice_file, 'r')
    cice_lon_tmp = id.variables['TLON'][:,:]
    cice_lat_tmp = id.variables['TLAT'][:,:]
    id.close()
    # Wrap the periodic boundary by 1 cell
    cice_lon = ma.empty([size(cice_lon_tmp,0), size(cice_lon_tmp,1)+1])
    cice_lat = ma.empty([size(cice_lat_tmp,0), size(cice_lat_tmp,1)+1])
    cice_lon[:,:-1] = cice_lon_tmp
    cice_lon[:,-1] = cice_lon_tmp[:,0]
    cice_lat[:,:-1] = cice_lat_tmp
    cice_lat[:,-1] = cice_lat_tmp[:,0]
    # Get averages for February and September
    # Start with first year just to initialise the arrays with the right size
    print '...monthly average for ' + str(start_year)
    cice_feb_tmp = monthly_avg_cice(cice_file, 'aice', shape(cice_lon_tmp), 1, instance=1)
    cice_sep_tmp = monthly_avg_cice(cice_file, 'aice', shape(cice_lon_tmp), 8, instance=1)
    # Loop over the rest of the years
    for year in range(start_year+1, end_year+1):
        print '...monthly average for ' + str(year)
        cice_feb_tmp = cice_feb_tmp + monthly_avg_cice(cice_file, 'aice', shape(cice_lon_tmp), 1, instance=year-start_year+1)
        cice_sep_tmp = cice_sep_tmp + monthly_avg_cice(cice_file, 'aice', shape(cice_lon_tmp), 8, instance=year-start_year+1)
    # Convert from integrals to averages
    cice_feb_tmp = cice_feb_tmp/num_years
    cice_sep_tmp = cice_sep_tmp/num_years
    # Wrap the periodic boundary
    cice_feb = ma.empty(shape(cice_lon))
    cice_sep = ma.empty(shape(cice_lon))
    cice_feb[:,:-1] = cice_feb_tmp
    cice_sep[:,:-1] = cice_sep_tmp
    cice_feb[:,-1] = cice_feb_tmp[:,0]
    cice_sep[:,-1] = cice_sep_tmp[:,0]
    # Polar coordinates for plotting
    cice_x = -(cice_lat+90)*cos(cice_lon*deg2rad+pi/2)
    cice_y = (cice_lat+90)*sin(cice_lon*deg2rad+pi/2)

    # Now get extent timeseries
    # Read 5-day logfile
    cice_time_vals = []
    cice_extent_5day = []
    f = open(cice_log, 'r')
    f.readline()
    for line in f:
        try:
            cice_time_vals.append(float(line))
        except(ValueError):
            break
    for line in f:
        cice_extent_5day.append(float(line))
    f.close()
    # Convert time to Date objects
    cice_time = num2date(array(cice_time_vals)*365.25, units='days since 1992-01-01 00:00:00', calendar='gregorian')
    # Initialise integral arrays for monthly averages
    # Add an extra year because the simulation goes to the end of 2016
    cice_extent = zeros((num_years+1)*12)
    cice_ndays = zeros((num_years+1)*12)
    for t in range(size(cice_time)):
        # 5-day averages marked with the next day's date
        year = cice_time[t].year
        month = cice_time[t].month-1  # Convert to 0-indexed
        day = cice_time[t].day
        # Check for leap years
        leap_year = False
        if mod(year, 4) == 0:
            leap_year = True
            if mod(year, 100) == 0:
                leap_year = False
                if mod(year, 400) == 0:
                    leap_year = True
        if leap_year:
            end_day[1] = 29
        else:
            end_day[1] = 28
        if day-5 < start_day[month]:
            # Spills over into the previous month
            prev_month = mod(month-1, 12)
            # How many days does it spill over by?
            spill_days = start_day[month]-day+5
            # Should be between 1 and 5
            if spill_days < 1 or spill_days > 5:
                print 'Problem: spill_days is ' + str(spill_days)
                print 'Timestep ' + str(t+1)
                print 'Year ' + str(year)
                print 'Month ' + str(month+1)
                print 'Day ' + str(day)
                #return
            # Split between previous month and this month
            # First find indices to update
            if prev_month+1 == 12:
                # Spilled into previous year
                index_prev = (year-1-start_year)*12 + prev_month
            else:
                index_prev = (year-start_year)*12 + prev_month
            index = (year-start_year)*12 + month
            # Integrate
            cice_extent[index_prev] += cice_extent_5day[t]*spill_days
            cice_ndays[index_prev] += spill_days
            cice_extent[index] += cice_extent_5day[t]*(5-spill_days)
            cice_ndays[index] += 5-spill_days
        else:
            # Entirely within the month
            index = (year-start_year)*12 + month
            cice_extent[index] += cice_extent_5day[t]*5
            cice_ndays[index] += 5
    # Convert from integrals to averages
    cice_extent /= cice_ndays
    # Extract February and September
    cice_feb_extent = []
    cice_sep_extent = []
    for year in range(start_year, end_year+1):
        cice_feb_extent.append(cice_extent[(year-start_year)*12+1])
        cice_sep_extent.append(cice_extent[(year-start_year)*12+8])

    print 'Processing FESOM'

    # First build the grid
    elements, patches = make_patches(fesom_mesh_path, circumpolar, mask_cavities)
    # Get averages for February and September
    # Start with first year just to initialise the arrays with the right size
    print '...monthly average for ' + str(start_year)
    fesom_feb_nodes = monthly_avg(fesom_output_dir + fesom_file_head + '.' + str(start_year) + '.ice.mean.nc', 'area', 1)
    fesom_sep_nodes = monthly_avg(fesom_output_dir + fesom_file_head + '.' + str(start_year) + '.ice.mean.nc', 'area', 8)
    # Loop over the rest of the years
    for year in range(start_year+1, end_year+1):
        print '...monthly average for ' + str(year)
        fesom_feb_nodes = fesom_feb_nodes + monthly_avg(fesom_output_dir + fesom_file_head + '.' + str(year) + '.ice.mean.nc', 'area', 1)
        fesom_sep_nodes = fesom_sep_nodes + monthly_avg(fesom_output_dir + fesom_file_head + '.' + str(year) + '.ice.mean.nc', 'area', 8)
    # Convert from integrals to averages
    fesom_feb_nodes = fesom_feb_nodes/num_years
    fesom_sep_nodes = fesom_sep_nodes/num_years
    # Find element-averages
    fesom_feb = []
    fesom_sep = []
    for elm in elements:
        if not elm.cavity:
            # Average over 3 component nodes
            fesom_feb.append(mean(array([fesom_feb_nodes[elm.nodes[0].id], fesom_feb_nodes[elm.nodes[1].id], fesom_feb_nodes[elm.nodes[2].id]])))
            fesom_sep.append(mean(array([fesom_sep_nodes[elm.nodes[0].id], fesom_sep_nodes[elm.nodes[1].id], fesom_sep_nodes[elm.nodes[2].id]])))
    fesom_feb = array(fesom_feb)
    fesom_sep = array(fesom_sep)

    # Get extent timeseries
    # Read 5-day logfile
    fesom_extent_5day = []
    f = open(fesom_log, 'r')    
    f.readline()
    for line in f:
        fesom_extent_5day.append(float(line))
    f.close()
    # Initialise monthly arrays
    fesom_feb_extent = zeros(num_years)
    fesom_sep_extent = zeros(num_years)
    for year in range(start_year, end_year+1):
        # First timestep of year in 5-day logfile
        t0 = (year-start_year)*73
        # Monthly averages are hard-coded and ugly
        # Feburary: 4/5 of index 7, indices 8-11, and 4/5 of index 12
        fesom_feb_extent[year-start_year] = (fesom_extent_5day[t0+6]*4 + sum(fesom_extent_5day[t0+7:t0+11]*5) + fesom_extent_5day[t0+11]*4)/28.0
        # September: 2/5 of index 49, indices 50-54, 3/5 of index 55
        fesom_sep_extent[year-start_year] = (fesom_extent_5day[t0+48]*2 + sum(fesom_extent_5day[t0+49:t0+54]*5) + fesom_extent_5day[t0+54]*3)/30.0

    time_axis = arange(start_year, end_year+1)

    print 'Plotting'
    fig = figure(figsize=(20,10))
    gs1 = GridSpec(2, 3)
    gs1.update(left=0.12, right=0.7, wspace=0.05, hspace=0.05)
    # NSIDC, February
    ax = subplot(gs1[0, 0], aspect='equal')
    img = pcolor(nsidc_x, nsidc_y, nsidc_feb, vmin=0, vmax=1, cmap='jet')
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    ax.set_xticks([])
    ax.set_yticks([])
    title('NSIDC', fontsize=24)
    text(-39, 0, 'February', fontsize=24, ha='right')
    # MetROMS, February
    ax = subplot(gs1[0, 1], aspect='equal')
    img = pcolor(cice_x, cice_y, cice_feb, vmin=0, vmax=1, cmap='jet')
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    ax.set_xticks([])
    ax.set_yticks([])
    title('MetROMS', fontsize=24)
    # FESOM, February
    ax = subplot(gs1[0, 2], aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(fesom_feb)
    img.set_clim(vmin=0, vmax=1)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    ax.set_xticks([])
    ax.set_yticks([])
    title('FESOM (high-res)', fontsize=24)
    # Main title
    text(-170, 47, 'a) Sea ice concentration ('+str(start_year)+'-'+str(end_year)+' average)', fontsize=30)
    # NSIDC, September
    ax = subplot(gs1[1, 0], aspect='equal')
    img = pcolor(nsidc_x, nsidc_y, nsidc_sep, vmin=0, vmax=1, cmap='jet')
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    ax.set_xticks([])
    ax.set_yticks([])
    text(-39, 0, 'September', fontsize=24, ha='right')
    # MetROMS, September
    ax = subplot(gs1[1, 1], aspect='equal')
    img = pcolor(cice_x, cice_y, cice_sep, vmin=0, vmax=1, cmap='jet')
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    ax.set_xticks([])
    ax.set_yticks([])    
    # FESOM, September
    ax = subplot(gs1[1, 2], aspect='equal')
    img = PatchCollection(patches, cmap='jet')
    img.set_array(fesom_sep)
    img.set_clim(vmin=0, vmax=1)
    img.set_edgecolor('face')
    ax.add_collection(img)
    xlim([bdry1, bdry2])
    ylim([bdry3, bdry4])
    ax.set_xticks([])
    ax.set_yticks([])
    # Add a colourbar at the bottom
    cbaxes = fig.add_axes([0.12, 0.04, 0.3, 0.04])
    cbar = colorbar(img, orientation='horizontal', ticks=arange(0,1+0.25,0.25), cax=cbaxes, extend='min')
    cbar.ax.tick_params(labelsize=20)
    # Add extent timeseries on rightmost column, with more space for labels
    gs2 = GridSpec(2, 1)
    gs2.update(left=0.73, right=0.95, wspace=0.1, hspace=0.15)
    # February
    ax = subplot(gs2[0, 0])
    ax.plot(time_axis, nsidc_feb_extent, color='black', linewidth=2, linestyle='dashed')
    ax.plot(time_axis, cice_feb_extent, color='blue', linewidth=1.5)
    ax.plot(time_axis, fesom_feb_extent, color='green', linewidth=1.5)
    xlim([start_year, end_year])
    ax.set_yticks(arange(0,4+0.5,0.5))
    ax.set_yticklabels(['0', '', '1', '', '2', '', '3', '', '4'])
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    grid(True)
    title('b) Sea ice extent\n'+r'(million km$^2$)', fontsize=26)
    # Extent timeseries, September
    ax = subplot(gs2[1, 0])
    ax.plot(time_axis, nsidc_sep_extent, color='black', label='NSIDC', linewidth=2, linestyle='dashed')
    ax.plot(time_axis, cice_sep_extent, color='blue', label='MetROMS', linewidth=1.5)
    ax.plot(time_axis, fesom_sep_extent, color='green', label='FESOM (high-res)', linewidth=1.5)
    xlim([start_year, end_year])
    ax.set_yticks(arange(16,23+1,1))
    ax.set_yticklabels(['16', '', '18', '', '20', '', '22', ''])
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    grid(True)
    ax.legend(bbox_to_anchor=(1.04,-0.07), ncol=3, fontsize=20)
    fig.show()
    fig.savefig('aice_minmax_nsidc.png')


# Command-line interface
if __name__ == "__main__":

    cice_file = raw_input("Path to CICE output file containing data for the entire simulation: ")
    cice_log = raw_input("Path to CICE sea ice extent logfile: ")
    fesom_mesh_path = raw_input("Path to FESOM mesh directory: ")
    fesom_output_dir = raw_input("Path to FESOM output directory containing one ice.mean.nc file for each year: ")
    fesom_log = raw_input("Path to FESOM sea ice extent logfile: ")
    mip_aice_minmax_nsidc(cice_file, cice_log, fesom_mesh_path, fesom_output_dir, fesom_log)
    
         
    
        
        
