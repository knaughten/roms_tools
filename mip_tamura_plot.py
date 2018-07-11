from numpy import *
from matplotlib.pyplot import *

def mip_tamura_plot (log_file):

    # Longitude to split the plot
    lon_split = -65

    # Read log file
    f = open(log_file, 'r')
    # Skip first line (header for longitude array)
    f.readline()
    lon = []
    for line in f:
        try:
            lon.append(float(line))
        except(ValueError):
            # Reached the header for the next variable
            break
    metroms_data = []
    for line in f:
        try:
            metroms_data.append(float(line))
        except(ValueError):
            break
    fesom_data_lr = []
    for line in f:
        try:
            fesom_data_lr.append(float(line))
        except(ValueError):
            break
    fesom_data_hr = []
    for line in f:
        try:
            fesom_data_hr.append(float(line))
        except(ValueError):
            break
    tamura_data = []
    for line in f:
        tamura_data.append(float(line))
    lon = array(lon)
    metroms_data = array(metroms_data)
    fesom_data_lr = array(fesom_data_lr)
    fesom_data_hr = array(fesom_data_hr)
    tamura_data = array(tamura_data)

    # Split at 65W
    split_index = nonzero(lon > lon_split)[0][0]
    lon = concatenate((lon[split_index:], lon[:split_index]+360))
    metroms_data = concatenate((metroms_data[split_index:], metroms_data[:split_index]))
    fesom_data_lr = concatenate((fesom_data_lr[split_index:], fesom_data_lr[:split_index]))
    fesom_data_hr = concatenate((fesom_data_hr[split_index:], fesom_data_hr[:split_index]))
    tamura_data = concatenate((tamura_data[split_index:], tamura_data[:split_index]))

    # Plot
    fig = figure(figsize=(13,6))
    gs = GridSpec(1,1)
    gs.update(left=0.07, right=0.97, bottom=0.1, top=0.9)
    ax = subplot(gs[0,0])
    plot(lon, metroms_data, color='blue', label='MetROMS')
    plot(lon, fesom_data_lr, color='green', label='FESOM (low-res)')
    plot(lon, fesom_data_hr, color='magenta', label='FESOM (high-res)')
    plot(lon, tamura_data, color=(0.5,0.5,0.5), label='Observations (Tamura)')
    grid(True)
    legend()
    xlim([lon[0], lon[-1]])
    ax.set_xticks(arange(-60, -90+360+30, 30))
    ax.set_xticklabels([r'60$^{\circ}$W', r'30$^{\circ}$W', r'0$^{\circ}$', r'30$^{\circ}$E', r'60$^{\circ}$E', r'90$^{\circ}$E', r'120$^{\circ}$E', r'150$^{\circ}$E', r'180$^{\circ}$', r'150$^{\circ}$W', r'120$^{\circ}$W', r'90$^{\circ}$W'])
    xlabel('longitude', fontsize=16)
    ylim([0,100])
    ylabel(r'10$^9$ m$^3$/y', fontsize=16)
    title('Sea ice production on continental shelf, 1992-2013 mean', fontsize=20)
    text(-40, 60, 'Weddell\nSea', fontsize=16, ha='center', va='center')
    text(60, 70, 'Prydz\nBay', fontsize=16, ha='center', va='center')
    text(115, 58, 'Australian Sector', fontsize=16, ha='center', va='center')
    text(180, 73, 'Ross\nSea', fontsize=16, ha='center', va='center')
    text(-105+360, 50, 'Amundsen\nSea', fontsize=16, ha='center', va='center')
    fig.show()
    fig.savefig('tamura_plot.png')

    print 'MetROMS total sea ice is ' + str((mean(metroms_data)-mean(tamura_data))/mean(tamura_data)*100) + '% higher than Tamura'
    print 'FESOM low-res total sea ice is ' + str((mean(fesom_data_lr)-mean(tamura_data))/mean(tamura_data)*100) + '% higher than Tamura'
    print 'FESOM high-res total sea ice is ' + str((mean(fesom_data_hr)-mean(tamura_data))/mean(tamura_data)*100) + '% higher than Tamura'

# Command-line interface
if __name__ == "__main__":

    log_file = raw_input("Logfile from mip_tamura_binning.py: ")
    mip_tamura_plot(log_file)
    
