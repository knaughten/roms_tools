from numpy import *
from matplotlib.pyplot import *

def mip_scatterplot (roms_logfile, fesom_logfile_lr, fesom_logfile_hr):

    # Year simulations start
    year_start = 1992
    # Years to average over
    calc_start = 2002
    calc_end = 2016
    # Number of output steps per year in FESOM
    peryear = 365/5
    # Name of each ice shelf
    names = ['Larsen D', 'Larsen C', 'Wilkins & George VI & Stange', 'Filchner-Ronne', 'Abbot', 'Pine Island', 'Thwaites', 'Dotson', 'Getz', 'Nickerson', 'Sulzberger', 'Mertz', 'Totten & Moscow University', 'Shackleton', 'West', 'Amery', 'Prince Harald', 'Baudouin & Borchgrevink', 'Lazarev', 'Nivl', 'Fimbul & Jelbart & Ekstrom', 'Brunt & Riiser-Larsen', 'Ross']
    # Observed mass loss (Rignot 2013) and uncertainty for each ice shelf, in Gt/y
    obs_massloss = [1.4, 20.7, 135.4, 155.4, 51.8, 101.2, 97.5, 45.2, 144.9, 4.2, 18.2, 7.9, 90.6, 72.6, 27.2, 35.5, -2, 21.6, 6.3, 3.9, 26.8, 9.7, 47.7]
    obs_massloss_error = [14, 67, 40, 45, 19, 8, 7, 4, 14, 2, 3, 3, 8, 15, 10, 23, 3, 18, 2, 2, 14, 16, 34]
    num_shelves = len(obs_massloss)
    # Order of indices for the ice shelves to be plotted on the x-axis (0-based)
    order = [3, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 22, 10, 9, 8, 7, 6, 5, 4, 2, 1, 0]

    # Read ROMS logfile
    roms_time = []
    f = open(roms_logfile, 'r')
    # Skip the first line (header for time array)
    f.readline()
    for line in f:
        try:
            roms_time.append(float(line))
        except(ValueError):
            # Reached the header for the next variable
            break       
    # Set up array for mass loss values at each ice shelf
    roms_massloss_ts = empty([num_shelves, len(roms_time)])
    # Skip total mass loss
    for line in f:
        try:
            tmp = float(line)
        except(ValueError):
            # Reaced the header for the first ice shelf
            break
    index = 0
    # Loop over ice shelves
    while index < num_shelves:
        t = 0
        for line in f:
            try:
                roms_massloss_ts[index, t] = float(line)
                t += 1
            except(ValueError):
                # Reached the header for the next ice shelf
                break
        index +=1
    f.close()
    # Add start year to ROMS time array
    roms_time = array(roms_time) + year_start
    # Average between given years
    t_start = nonzero(roms_time >= calc_start)[0][0]
    if calc_end == 2016:
        t_end = size(roms_time)
    else:
        t_end = nonzero(roms_time >= calc_end+1)[0][0]
    roms_massloss = mean(roms_massloss_ts[:,t_start:t_end], axis=1)

    # Read FESOM timeseries
    # Low-res
    f = open(fesom_logfile_lr, 'r')
    # Skip the first line (header)
    f.readline()
    # Skip total mass loss
    num_time = 0
    for line in f:
        try:
            tmp = float(line)
            num_time += 1
        except(ValueError):
            # Reached the header for the next variable
            break
    # Set up array for mass loss values at each ice shelf
    fesom_massloss_ts_lr = empty([num_shelves, num_time])
    # Loop over ice shelves
    index = 0
    while index < num_shelves:
        t = 0
        for line in f:
            try:
                fesom_massloss_ts_lr[index,t] = float(line)
                t += 1
            except(ValueError):
                # Reached the header for the next ice shelf
                break
        index += 1
    f.close()
    # Average between given years
    fesom_massloss_lr = mean(fesom_massloss_ts_lr[:,peryear*(calc_start-year_start):peryear*(calc_end+1-year_start)], axis=1)
    # Repeat for high-res
    f = open(fesom_logfile_hr, 'r')
    f.readline()
    num_time = 0
    for line in f:
        try:
            tmp = float(line)
            num_time += 1
        except(ValueError):
            break
    fesom_massloss_ts_hr = empty([num_shelves, num_time])
    index = 0
    while index < num_shelves:
        t = 0
        for line in f:
            try:
                fesom_massloss_ts_hr[index,t] = float(line)
                t += 1
            except(ValueError):
                break
        index += 1
    f.close()
    fesom_massloss_hr = mean(fesom_massloss_ts_hr[:,peryear*(calc_start-year_start):peryear*(calc_end+1-year_start)], axis=1)

    # Figure out error values, in correct order for plotting
    roms_error = []
    fesom_error_lr = []
    fesom_error_hr = []
    labels = []
    for index in range(num_shelves):
        obs_min = obs_massloss[order[index]] - obs_massloss_error[order[index]]
        obs_max = obs_massloss[order[index]] + obs_massloss_error[order[index]]
        if roms_massloss[order[index]] < obs_min:
            roms_error.append(roms_massloss[order[index]] - obs_min)
        elif roms_massloss[order[index]] > obs_max:
            roms_error.append(roms_massloss[order[index]] - obs_max)
        else:
            roms_error.append(0)
        if fesom_massloss_lr[order[index]] < obs_min:
            fesom_error_lr.append(fesom_massloss_lr[order[index]] - obs_min)
        elif fesom_massloss_lr[order[index]] > obs_max:
            fesom_error_lr.append(fesom_massloss_lr[order[index]] - obs_max)
        else:
            fesom_error_lr.append(0)
        if fesom_massloss_hr[order[index]] < obs_min:
            fesom_error_hr.append(fesom_massloss_hr[order[index]] - obs_min)
        elif fesom_massloss_hr[order[index]] > obs_max:
            fesom_error_hr.append(fesom_massloss_hr[order[index]] - obs_max)
        else:
            fesom_error_hr.append(0)
        labels.append(names[order[index]])

    # Plot
    fig = figure(figsize=(10,7))
    gs = GridSpec(1,1)
    gs.update(left=0.1, right=0.9, bottom=0.4, top=0.9)
    ax = subplot(gs[0,0])
    # Alternate background between white and light blue to split up regions
    axvspan(0.5, 6.5, facecolor='b', alpha=0.1)
    axvspan(7.5, 11.5, facecolor='b', alpha=0.1)
    axvspan(14.5, 18.5, facecolor='b', alpha=0.1)
    axvspan(20.5, 23, facecolor='b', alpha=0.1)
    # Region labels
    text(0, 30, 'FR', fontsize=14, ha='center')
    text(3.5, 30, 'EWed', fontsize=14, ha='center')
    text(7, -20, 'Am', fontsize=14, ha='center')
    text(9.5, 30, 'Aus', fontsize=14, ha='center')
    text(13, 30, 'RS', fontsize=14, ha='center')
    text(16.5, 30, 'AS', fontsize=14, ha='center')
    text(19.5, 30, 'BS', fontsize=14, ha='center')
    text(21.5, 30, 'Lr', fontsize=14, ha='center')
    plot(range(-1,num_shelves+1), zeros(num_shelves+2), color='black', linewidth=3)
    plot(range(num_shelves), roms_error, 'o', color=(0.08, 0.4, 0.79), ms=10, label='MetROMS')
    plot(range(num_shelves), fesom_error_lr, 'o', color=(0.73, 0.06, 0.69), ms=10, label='FESOM (low-res)')
    plot(range(num_shelves), fesom_error_hr, 'o', color=(0.06, 0.73, 0.1), ms=10, label='FESOM (high-res)')
    grid(True)    
    xlim([-0.5, num_shelves-0.5])
    xticks(range(num_shelves), labels, rotation=90)
    ylabel('Gt/y', fontsize=14)
    title('Bias in Ice Shelf Basal Mass Loss', fontsize=18)
    # Make legend
    legend(numpoints=1,loc='lower left')
    setp(gca().get_legend().get_texts(), fontsize='13')
    fig.show()
    fig.savefig('scatterplot.png')


# Command-line interface
if __name__ == "__main__":

    roms_logfile = raw_input("Path to ROMS logfile from timeseries_massloss.py: ")
    fesom_logfile_lr = raw_input("Path to FESOM low-res logfile from timeseries_massloss.py: ")
    fesom_logfile_hr = raw_input("Path to FESOM high-res logfile from timeseries_massloss.py: ")
    mip_scatterplot(roms_logfile, fesom_logfile_lr, fesom_logfile_hr)

