from os import system

def average_erainterim ():

    head_in = '../ROMS-CICE-MCT/data/ERA_Interim/originals/FC_'
    tail_in = '_unlim_orig.nc'
    head_out = '../ROMS-CICE-MCT/data/ERA_Interim/monthly/FC_'
    tail_out = '_monthly_orig.nc'
    steps_per_day = 2

    days_per_month = [31,28,31,30,31,30,31,31,30,31,30,31]
    first_year = 1992
    last_year = 2005
    leap_years = range(1976,last_year+1,4)

    for year in range(first_year,last_year+1):
        posn = 0
        for month in range(12):
            print('Processing month '+str(month+1)+' of '+str(year))
            if month == 1 and year in leap_years:
                days = 29
            else:
                days = days_per_month[month]
            num_steps = days*an_steps_per_day
            month_str = str(month+1)
            if month < 9:
                month_str = str(0) + month_str
            in_file = an_head_in + str(year) + an_tail_in
            tmp_file = 'tmp' + month_str + '.nc'
            system('ncra -d time,' + str(posn) + ',' + str(posn+num_steps-1) + ' ' + in_file + ' ' + tmp_file)
            posn = posn + num_steps
        print('Concatenating '+str(year))
        out_file = an_head_out + str(year) + an_tail_out
        system('ncrcat tmp*.nc ' + out_file)
        system('rm tmp*.nc')

    for year in range(first_year,last_year+1):
        posn = 0
        for month in range(12):
            print('Processing month '+str(month+1)+' of '+str(year))
            if month == 1 and year in leap_years:
                days = 29
            else:
                days = days_per_month[month]
            num_steps = days*fc_steps_per_day
            month_str = str(month+1)
            if month < 9:
                month_str = str(0) + month_str
            in_file = fc_head_in + str(year) + fc_tail_in
            tmp_file = 'tmp' + month_str + '.nc'
            system('ncra -d time,' + str(posn) + ',' + str(posn+num_steps-1) + ' ' + in_file + ' ' + tmp_file)
            posn = posn + num_steps
        print('Concatenating '+str(year))
        out_file = fc_head_out + str(year) + fc_tail_out
        system('ncrcat tmp*.nc ' + out_file)
        system('rm tmp*.nc')


if __name__ == "__main__":

    average_erainterim()
            
