from numpy import *
from matplotlib.pyplot import *

# Plot MetROMS and FESOM timeseries together, for Drake Passage transport, total
# Antarctic sea ice area and volume, Antarctic sea ice extent, and basal mass
# loss for major ice shelves.
# Include the range of observations for Drake Passage transport and ice shelf
# mass loss.
# Input:
# roms_dir = path to ROMS directory containing logfiles from timeseries_dpt.py,
#            timeseries_seaice.py, timeserires_seaice_extent.py, and 
#            timeseries_massloss.py. It is assumed they are saved with the
#            filenames dpt.log, seaice.log, seaice_extent.log, and massloss.log.
# fesom_dir = path to FESOM directory containing logfiles from the equivalent
#             "fesomtools" scripts with the same names. It is assumed they are
#             also saved with the filenames dpt.log, seaice.log,
#             seaice_extent.log, and massloss.log.
# annual = optional boolean indicating to calculate annual averages of Drake
#          Passage transport and ice shelf mass loss, and to not bother with
#          sea ice area, volume, and extent.
def mip_timeseries (roms_dir, fesom_dir, annual=False):

    # Averaging period (days)
    days_per_output = 5
    year_start = 1992
    year_end = 2016
    # Titles for plotting
    model_titles = ['MetROMS', 'FESOM']
    # Colours for plotting
    model_colours = ['blue', 'green']

    # Titles for each ice shelf
    shelf_names = ['All Ice Shelves', 'Larsen D Ice Shelf', 'Larsen C Ice Shelf', 'Wilkins & George VI & Stange Ice Shelves', 'Ronne-Filchner Ice Shelf', 'Abbot Ice Shelf', 'Pine Island Glacier Ice Shelf', 'Thwaites Ice Shelf', 'Dotson Ice Shelf', 'Getz Ice Shelf', 'Nickerson Ice Shelf', 'Sulzberger Ice Shelf', 'Mertz Ice Shelf', 'Totten & Moscow University Ice Shelves', 'Shackleton Ice Shelf', 'West Ice Shelf', 'Amery Ice Shelf', 'Prince Harald Ice Shelf', 'Baudouin & Borchgrevink Ice Shelves', 'Lazarev Ice Shelf', 'Nivl Ice Shelf', 'Fimbul & Jelbart & Ekstrom Ice Shelves', 'Brunt & Riiser-Larsen Ice Shelves', 'Ross Ice Shelf']
    # Beginnings of figure names for each ice shelf
    fig_names = ['total_massloss', 'larsen_d', 'larsen_c', 'wilkins_georgevi_stange', 'ronne_filchner', 'abbot', 'pig', 'thwaites', 'dotson', 'getz', 'nickerson', 'sulzberger', 'mertz', 'totten_moscowuni', 'shackleton', 'west', 'amery', 'princeharald', 'baudouin_borchgrevink', 'lazarev', 'nivl', 'fimbul_jelbart_ekstrom', 'brunt_riiserlarsen', 'ross']
    num_shelves = len(shelf_names)

    # Bounds of observations for Drake Passage transport
    dpt_low = 107  # Lower bound of range of estimates from Cunningham et al 2003, doi:10.1029/2001JC001147: 134 +/- 15-27 Sv
    dpt_high = 173.3  # Donohue et al 2016, doi:10.1002/2016GL070319
    # Observed mass loss (Rignot 2013) and uncertainty for each ice shelf, in
    # Gt/y
    obs_massloss = [1325, 1.4, 20.7, 135.4, 155.4, 51.8, 101.2, 97.5, 45.2, 144.9, 4.2, 18.2, 7.9, 90.6, 72.6, 27.2, 35.5, -2, 21.6, 6.3, 3.9, 26.8, 9.7, 47.7]
    obs_massloss_error = [235, 14, 67, 40, 45, 19, 8, 7, 4, 14, 2, 3, 3, 8, 15, 10, 23, 3, 18, 2, 2, 14, 16, 34]

    # Bounds for plotting Drake Passage transport
    dpt_bounds = [80, 200]
    # Bounds for plotting each ice shelf mass loss
    massloss_bounds_low = [0, -15, -100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, 0, 0, 0, 0, -20, 0]
    massloss_bounds_high = [2000, 30, 300, 250, 600, 100, 120, 110, 55, 250, 35, 50, 30, 100, 120, 75, 850, 25, 150, 40, 25, 120, 175, 500]

    # Drake Passage transport
    # Read ROMS timeseries
    roms_time = []
    roms_dpt = []
    f = open(roms_dir + 'dpt.log', 'r')
    # Skip first line (header for time array)
    f.readline()
    for line in f:
        try:
            roms_time.append(float(line))
        except(ValueError):
            # Reached the header for the next variable
            break
    for line in f:
        roms_dpt.append(float(line))
    f.close()
    # Add start year to ROMS time array
    roms_time = array(roms_time) + year_start
    # Read FESOM timeseries
    fesom_dpt = []
    f = open(fesom_dir + 'dpt.log', 'r')
    # Skip header
    f.readline()
    for line in f:
        fesom_dpt.append(float(line))
    f.close()
    # Make FESOM time array (note that FESOM neglects leap years in its output)
    fesom_time = arange(len(fesom_dpt))*days_per_output/365. + year_start
    if annual:
        # Annually average ROMS data
        roms_dpt_avg = roms_annual_avg(roms_time, roms_dpt, year_start)
        # Annually average FESOM data
        fesom_dpt_avg = fesom_annual_avg(fesom_dpt, days_per_output)
        # Make new time arrays
        roms_time_annual = arange(len(roms_dpt_avg)) + year_start
        fesom_time_annual = arange(len(fesom_dpt_avg)) + year_start
        # Find last timestep
        max_time = max(roms_time_annual[-1], fesom_time_annual[-1])
    else:
        max_time = max(roms_time[-1], fesom_time[-1])
    # Plot
    fig, ax = subplots(figsize=(10,6))
    if annual:
        ax.plot(roms_time_annual, roms_dpt_avg, label=model_titles[0], color=model_colours[0], linewidth=2)
        ax.plot(fesom_time_annual, fesom_dpt_avg, label=model_titles[1], color=model_colours[1], linewidth=2)
    else:
        ax.plot(roms_time, roms_dpt, label=model_titles[0], color=model_colours[0], linewidth=1)
        ax.plot(fesom_time, fesom_dpt, label=model_titles[1], color=model_colours[1], linewidth=1)
    # Add lines for range of observations
    ax.axhline(dpt_low, color='red', linestyle='dashed', linewidth=2, label='observations')
    ax.axhline(dpt_high, color='red', linewidth=2, linestyle='dashed')
    if annual:
        title('Drake Passage Transport (annually averaged)', fontsize=18)
    else:
        title('Drake Passage Transport', fontsize=18)
    xlabel('Year', fontsize=14)
    ylabel('Sv', fontsize=14)
    xlim([year_start, max_time])
    ylim([dpt_bounds[0], dpt_bounds[1]])
    grid(True)
    # Move plot over to make room for legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
    # Make legend
    ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    if annual:
        fig.savefig('drakepsgtrans_avg.png')
    else:
        fig.savefig('drakepsgtrans.png')

    # Sea ice area and volume
    if not annual:
        # Read ROMS timeseries
        roms_icearea = []
        roms_icevolume = []
        f = open(roms_dir + 'seaice.log', 'r')
        # Skip the first line (header for time array)
        f.readline()
        for line in f:
            # Skip the time values (already have them from DPT)
            try:
                tmp = float(line)
            except(ValueError):
                # Reached the header for the next variable
                break
        # Sea ice area
        for line in f:
            try:
                roms_icearea.append(float(line))
            except(ValueError):
                break
        # Sea ice volume
            roms_icevolume.append(float(line))
        f.close()
        # Read FESOM timeseries
        fesom_icearea = []
        fesom_icevolume = []
        f = open(fesom_dir + 'seaice.log', 'r')
        f.readline()
        for line in f:
            try:
                fesom_icearea.append(float(line))
            except(ValueError):
                break
        for line in f:
            fesom_icevolume.append(float(line))
        f.close()
        # Plot sea ice area
        fig, ax = subplots(figsize=(10,6))
        ax.plot(roms_time, roms_icearea, label=model_titles[0], color=model_colours[0], linewidth=1)
        ax.plot(fesom_time, fesom_icearea, label=model_titles[1], color=model_colours[1], linewidth=1)
        title('Antarctic Sea Ice Area', fontsize=18)
        xlabel('Year', fontsize=14)
        ylabel(r'million km$^2$', fontsize=14)
        xlim([year_start, max_time])
        grid(True)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
        fig.savefig('seaice_area.png')
        # Plot sea ice volume
        fig, ax = subplots(figsize=(10,6))
        ax.plot(roms_time, roms_icevolume, label=model_titles[0], color=model_colours[0], linewidth=1)
        ax.plot(fesom_time, fesom_icevolume, label=model_titles[1], color=model_colours[1], linewidth=1)
        title('Antarctic Sea Ice Volume', fontsize=18)
        xlabel('Year', fontsize=14)
        ylabel(r'million km$^3$', fontsize=14)
        xlim([year_start, max_time])
        grid(True)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
        fig.savefig('seaice_volume.png')

    # Sea ice extent
    if not annual:
        # Read ROMS timeseries
        roms_iceextent = []
        f = open(roms_dir + 'seaice_extent.log', 'r')
        # Skip the first line (header for time array)
        f.readline()
        for line in f:
            # Skip the time values (already have them from DPT)
            try:
                tmp = float(line)
            except(ValueError):
                # Reached the header for the next variable
                break
        # Sea ice extent
        for line in f:
            roms_iceextent.append(float(line))
        f.close()
        # Read FESOM timeseries
        fesom_iceextent = []
        f = open(fesom_dir + 'seaice_extent.log', 'r')
        f.readline()
        for line in f:
            fesom_iceextent.append(float(line))
        f.close()
        # Plot
        fig, ax = subplots(figsize=(10,6))
        ax.plot(roms_time, roms_iceextent, label=model_titles[0], color=model_colours[0], linewidth=1)
        ax.plot(fesom_time, fesom_iceextent, label=model_titles[1], color=model_colours[1], linewidth=1)
        title('Antarctic Sea Ice Extent', fontsize=18)
        xlabel('Year', fontsize=14)
        ylabel(r'million km$^2$', fontsize=14)
        xlim([year_start, max_time])
        grid(True)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
        fig.savefig('seaice_extent.png')

    # Ice shelf mass loss
    # Read ROMS timeseries
    roms_massloss = empty([num_shelves, size(roms_time)])
    f = open(roms_dir + 'massloss.log', 'r')
    # Skip the first line (header for time array)
    f.readline()
    for line in f:
        # Skip the time values (already have them from DPT)
        try:
            tmp = float(line)
        except(ValueError):
            # Reached the header for the next variable
            break
    # Loop over ice shelves
    index = 0
    while index < num_shelves:
        t = 0
        for line in f:
            try:
                roms_massloss[index, t] = float(line)
                t += 1
            except(ValueError):
                # Reached the header for the next ice shelf
                break
        index += 1
    f.close()
    # Read FESOM timeseries
    fesom_massloss = empty([num_shelves, size(fesom_time)])
    f = open(fesom_dir + 'massloss.log', 'r')
    f.readline()
    # Loop over ice shelves
    index = 0
    while index < num_shelves:
        t = 0
        for line in f:
            try:
                fesom_massloss[index, t] = float(line)
                t += 1
            except(ValueError):
                # Reached the header for the next ice shelf
                break
        index += 1
    f.close()
    if annual:
        # Annually average ROMS data for each ice shelf
        roms_massloss_avg = empty([num_shelves, size(roms_time_annual)])
        for index in range(num_shelves):
            roms_massloss_avg[index,:] = roms_annual_avg(roms_time, roms_massloss[index,:], year_start)
        # Annually average FESOM data for each ice shelf
        fesom_massloss_avg = empty([num_shelves, size(fesom_time_annual)])
        for index in range(num_shelves):
            fesom_massloss_avg[index,:] = fesom_annual_avg(fesom_massloss[index,:], days_per_output)
    # One plot for each ice shelf
    for index in range(num_shelves):
        # Calculate range of observations
        massloss_low = obs_massloss[index] - obs_massloss_error[index]
        massloss_high = obs_massloss[index] + obs_massloss_error[index]
        fig, ax = subplots(figsize=(10,6))
        if annual:
            ax.plot(roms_time_annual, roms_massloss_avg[index,:], label=model_titles[0], color=model_colours[0], linewidth=2)
            ax.plot(fesom_time_annual, fesom_massloss_avg[index,:], label=model_titles[1], color=model_colours[1], linewidth=2)
        else:
            ax.plot(roms_time, roms_massloss[index,:], label=model_titles[0], color=model_colours[0], linewidth=1)
            ax.plot(fesom_time, fesom_massloss[index,:], label=model_titles[1], color=model_colours[1], linewidth=1)
        # Add lines for range of observations
        ax.axhline(massloss_low, color='red', linestyle='dashed', linewidth=2, label='observations')
        ax.axhline(massloss_high, color='red', linewidth=2, linestyle='dashed')
        if annual:
            title(shelf_names[index] + '\nBasal Mass Loss (Annually Averaged)', fontsize=18)
        else:
            title(shelf_names[index] + '\nBasal Mass Loss', fontsize=18)
        xlabel('Year', fontsize=14)
        ylabel('Gt/y', fontsize=14)
        xlim([year_start, max_time])
        if not annual:
            ylim([massloss_bounds_low[index], massloss_bounds_high[index]])
        grid(True)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
        if annual:
            fig.savefig(fig_names[index] + '_avg.png')
        else:
            fig.savefig(fig_names[index] + '.png')
    
        
# Annually average the given ROMS data. Note that ROMS takes leap years into
# account for its output, i.e. the 5-day averaging period holds true throughout
# the entire simulation with no days skipped.
# Input:
# time = 1D array of ROMS time values, in years, with year_start added (eg
#        1992.001, 1993.50) assuming year length of 365.25 days
# data = 1D array of ROMS data corresponding to time array
# year_start = integer containing the first year to process
# Output: data_avg = 1D array of annually averaged data
def roms_annual_avg (time, data, year_start):

    data_avg = []
    year = year_start
    # Loop over years until we run out of data
    while True:
        # Find the first index after the beginning of this year
        tmp1 = nonzero(time >= year)[0]
        if len(tmp1)==0:
            # No such index, we have run out of data
            break
        t_start = tmp1[0]
        # Find the first index after the beginning of next year
        tmp2 = nonzero(time >= year+1)[0]
        if len(tmp2)==0:
            # No such index, but we might not have run out of data
            # eg simulation that ends on 31 December 2016
            t_end = len(time)
        else:
            t_end = tmp2[0]
        if t_end - t_start < 73:
            # This is a partial year, don't count it
            break
        # Calculate mean between these bounds
        data_avg.append(mean(array(data[t_start:t_end])))
        # Increment year
        year += 1
    return data_avg


# Annually average the given FESOM data. Note that FESOM neglects leap years
# in its output, i.e. the 5-day averaging period will skip 1 day every 4 years,
# so every year has exactly 365/5 = 73 batches of output.
# Input:
# data = 1D array of FESOM data 
# days_per_output = averaging period for data
# Output: data_avg = 1D array of annually averaged data
def fesom_annual_avg (data, days_per_output):

    peryear = 365/days_per_output
    data_avg = []
    year = 0
    # Loop over years
    while True:
        if len(data) < peryear*(year+1):
            # Run out of data
            # FESOM output comes in annual files so we don't have to worry about
            # partial years like in ROMS
            break
        # Calculate mean for this year
        data_avg.append(mean(array(data[peryear*year:peryear*(year+1)])))
        # Increment year
        year += 1
    return data_avg


# Command-line interface
if __name__ == "__main__":

    roms_dir = raw_input("Path to directory containing MetROMS logfiles: ")
    fesom_dir = raw_input("Path to directory containing FESOM logfiles: ")
    time_flag = int(raw_input("Full timeseries (1) or annual averages (2)? "))
    if time_flag == 1:
        annual = False
    elif time_flag == 2:
        annual = True
    mip_timeseries(roms_dir, fesom_dir, annual)
        

    
            
        
    
    
