#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib import dates
import netCDF4 as cdf
import numpy as np
import sys
import time
import datetime as dt  # Python standard library datetime  module
import pdb
import calendar # for days in months
import argparse # argument passing to script


#  Reads from input file (e.g. model file sites_2007.nc say), and provides
#  various outputs for a chosen pollutant. Also does plots for several sites
#  if wanted
        # SITES.dmax     - daily max values  
        # SITES.dmean    - daily mean values  
        # SITES.FREQhrly - all FREQ hourly values
        # SITES.vals     - all values, printed (24/FREQ) hourly values per
        #                  line
        # SITES.mmean    - monthly mean values  

# ===========================================================================
# function to print sites
# ===========================================================================
def print_ask_sites(site_names, print_sites, additional):
    "if print_sites is true, print out the list of site_names, then let \
    one be selected and return that if it exists. Flag if asking for \
    additional sites or not"
    n_sites = len(site_names)
    to_print = ['{:4} {}'.format(i, x) for i, x in enumerate(site_names)]
    if print_sites:
        print '\n'.join(to_print)
    while True:
        if additional:
            chosen_site = int(raw_input('Choose additional site\n'))
        else:
            chosen_site = int(raw_input('Choose site\n'))
        if (chosen_site < 0) or (chosen_site >= n_sites):
            print '\n'.join(to_print)
            print 'NO SITE OF THAT NAME, CHOOSE AGAIN!'
        else:
            break
    return(chosen_site)
# ===========================================================================




parser = argparse.ArgumentParser()
parser.add_argument('-i', '--ifile', help = 'Input file sites_yyyy.nc',
    required = True)
parser.add_argument('-s', '--sitenames', help = 'Name of site(s). Provide ' +
    'several sites in a comma-separated list.')
parser.add_argument('-v', '--varname', help = 'Name of variable in nc file')
parser.add_argument('-n', '--nsites', help = 'How many additional sites to ' +
    'use. You will be asked later which sites you want to add.' , type = int)
parser.add_argument('-pt', '--plottype', choices = ['all', 'daymeans', 
    'daymaxs', 'monmeans'], help = ('Type of plot wanted. Choose ' +
    'between "all": all data (e.g. hourly), "daymeans": daily ' +
    'means, "daymaxs": daily maxima, "monmeans": monthly means'))
parser.add_argument('-pr', '--plotrange', help = 'Which time range do you ' +
    'want to be plotted? Either enter "entire" or the start and end time ' +
    'separated by a "-". The format should suit the temporal resolution, ' +
    'so yyyymmddhh for "all" data, yyyymmdd for daily data and yyyymm for ' +
    'monthly data. E.g. 2012030405-2012080910')
args = parser.parse_args()
print 'ARGS:'
print args

infile = args.ifile
print 'Input argument: {}'.format(infile)

if args.sitenames is not None:
    if ',' in args.sitenames:
        wsite = args.sitenames.split(',')
        print 'Sites wanted: {}'.format(', '.join(wsite))
    else:
        wsite = [args.sitenames]
        print 'Site wanted: {}'.format(wsite[0])
else:
    wsite = []

if args.varname is not None:
    wpoll = args.varname
    print 'Poll wanted: {}'.format(wpoll)
else:
    wpoll = ''

if args.nsites is not None:
    sites_add = args.nsites
    if sites_add < 0:
        sys.exit('EXIT: Argument -n (--nsites) has to be an integer >= 0!')
    if len(wsite) > 1:
        sys.exit('EXIT: When calling the program, please either provide ' +
            'all sites to the argument -s (--sitenames) OR provide up to ' +
            'one site to -s and the number of additional sites to -n ' +
            '(--nsites). Do not mix them up.')
elif len(wsite) > 1:
    sites_add = -20
else:
    sites_add = -10

if args.plottype is None:
    which_plot = None
else:
    which_plot = args.plottype
    print 'Plot type wanted: {}'.format(which_plot)

if args.plotrange is None:
    which_time = None
    date_beg = ''
else:
    if args.plotrange == 'entire':
        which_time = 1
    else:
        which_time = 2
        date_beg = args.plotrange

inCDF = cdf.Dataset(infile,'r',format='NETCDF4')

nsites = getattr(inCDF, 'Number_of_stations_defined')
# site names are saved in dimension (10,36) with single letters as
#  elements of a list
sitename = [str(''.join(inCDF.variables['station_name'][a,:]).strip()) 
    for a in range(nsites)]

freq = int(getattr(inCDF, 'Number_of_hours_bewtween_outputs'))

# wsite: name of site
# wanted_site, all_sites: index of that site
all_sites = []
if (len(wsite) == 1) and (wsite[0] in sitename):
    wanted_site = sitename.index(wsite[0])
    all_sites.append(wanted_site)
    print 'selected site: {}'.format(wsite[0])
elif (len(wsite) > 1):
    for si in wsite:
        if si in sitename:
            wanted_site = sitename.index(si)
            all_sites.append(wanted_site)
            print 'selected site: {}'.format(si)
else:
    wanted_site = print_ask_sites(sitename, print_sites = True,
        additional = False)
    all_sites.append(wanted_site)
    wsite = [sitename[wanted_site]]
print 'NSITES {}, FREQ {}, WANTED_SITEXX {}'.format(
    nsites, freq, ('{} '*len(all_sites)).format( *['{}: {}'.format(
        all_sites[i], wsite[i]) for i in range(len(all_sites))] ))
time.sleep(1)

# in case additional sites weren't selected in the beginning, ask for those
if sites_add == -10:
    print ('Do you want to use more sites than only {}? Enter either "0" ' +
        'for only this site,\nor the number of additional sites that you ' +
        'want to be used!').format(wsite[0])
    while True:
        try:
            sites_add = int(raw_input('Your choice: '))
            break
        except ValueError:
            print 'Only integers allowed, try again!'

if sites_add > 0:
    all_sites = [wanted_site]
    print ('You selected {} additional sites. Now choose them.').format(
        sites_add)
    time.sleep(1)
    for i in range(1,sites_add+1):
        is_first = (i == 1)
        new_site = print_ask_sites(sitename, print_sites = is_first,
            additional = True)
        all_sites.append(new_site)
        print 'will also use site {}: {}'.format(new_site, 
            sitename[new_site])

# remove possible duplicate sites
all_sites = list(set(all_sites))
n_selsites = len(all_sites)

# ===========================================================================
# now looking at variables
wanted_spec = -1

nspec = len(inCDF.variables)
specname = [var for var in inCDF.variables]

# wpoll: name of variable (species, pollutant)
# wanted_spec: index
if (len(wpoll) > 0) and (wpoll in specname):
    wanted_spec = specname.index(wpoll)
else:
    to_print = ['{:4} {}'.format(i, x) for i, x in enumerate(specname)]
    print '\n'.join(to_print)
    wanted_spec = int(raw_input('Choose species\n'))
    if (wanted_spec < 0) or (wanted_spec >= nspec):
        print '\n'.join(to_print)
        raise ValueError('NO POLL OF THAT NAME, START AGAIN!')
    wpoll = specname[wanted_spec]

# this is the selected variable for the selected sites
xin = inCDF.variables[wpoll][:,all_sites]

unit = inCDF.variables[wpoll].units

print_all_sites = ['{} ({})'.format(i, sitename[i]) for i in all_sites]
print ('Sites' + (' {}'*n_selsites).format(*print_all_sites) + 
    ', Poll {} ({})'.format(wanted_spec, specname[wanted_spec]))

# Date
nctime = inCDF.variables['time'][:]
unit_time = inCDF.variables['time'].units # get unit 'days since ...'
try :
    cal_time = inCDF.variables['time'].calendar
except AttributeError : # Attribute doesn't exist
    cal_time = u"gregorian" # or standard
datevar = cdf.num2date(nctime, units = unit_time,
        calendar = cal_time)

# if first midnight is missing: copy first time step to 00:00
if datevar[0].hour == freq:
    print ("file doesn't begin at midnight: {} -> copy first " +
        "timestep").format(datevar[0].strftime('%Y%m%d%H'))
    xin = np.insert(xin, obj=0, values=xin[0,:], axis=0)
    datevar = np.insert(datevar, obj=0, 
        values = datevar[0] - dt.timedelta(hours = freq), axis=0)
elif datevar[0].hour == 0:
    pass
# otherwise remove the first time step(s) until the first midnight
else:
    print ("file begins neither at midnight nor one timestep later: {} -> " +
        'disregarding first timesteps until first midnight').format(
        datevar[0].strftime('%Y%m%d%H'))
    first_midnight = next(i for i, d in enumerate(datevar) if d.hour==0)
    datevar = datevar[first_midnight:]
    xin = xin[first_midnight:,:]
# end of data: disregard uncomplete last day
if datevar[-1].hour != (24-freq):
    print ("file doesn't end with correct hour: {}, disregarding " +
        "everything after (and including) the last midnight").format(
        datevar[-1].strftime('%Y%m%d%H'))
    for i, x in reversed(list(enumerate(datevar))):
        if x.hour == 0:
            last_midnight = i
            break
    datevar = datevar[:last_midnight]
    xin = xin[:last_midnight,:]

wanted_year = datevar[0].year

# headers for output files:
head_vals = ('idd yyyy mm dd ' + '{:13}'*(24/freq)).format(
    *['hh_{}'.format(i) for i in range(0,24,freq)])
head_hrly = 'idd yyyy mm dd hh {}hourly_val\n'.format(freq)
head_dmean = 'idd yyyy mm dd day_mean\n'
head_dmax = 'idd yyyy mm dd day_max\n'
head_mmean = 'yyyy mm mon_mean\n'


# arrays for all sites
all_day_means = []
all_day_maxs = []
all_month_means = []


# additionally will need number of days of last month
cal = calendar.Calendar()
last_month = cal.itermonthdays(wanted_year,datevar[-1].month)
days_last_month = []
[days_last_month.append(x) for x in last_month]
days_last_month = max(days_last_month)


# loop over sites
for j, site in enumerate(all_sites):

    wsite = sitename[site]
    # Output files, for each site
    # =======================================================================
    prefix = 'SITES_{}_{}'.format(wsite, wpoll)
    vals_file = open('{}.vals'.format(prefix), 'w')
    hrly_file = open('{}.{}hrly'.format(prefix, freq), 'w')
    dmean_file = open('{}.dmean'.format(prefix), 'w')
    dmax_file = open('{}.dmax'.format(prefix), 'w')
    mmean_file = open('{}.mmean'.format(prefix), 'w')
    # =======================================================================

    vals_file.write(head_vals + '\n')
    hrly_file.write(head_hrly)
    dmean_file.write(head_dmean)
    dmax_file.write(head_dmax)
    mmean_file.write(head_mmean)

    # go through variable's and site's data and write output files
    # =======================================================================
    old_day = datevar[0].day
    old_end = 0
    sum_days = 0
    daily_means = []
    daily_maxs = []
    # keep track of the days
    days_dates = []
    # also for monthly data
    old_month = datevar[0].month
    old_end_month = 0
    which_months = []
    monthly_means = []
    nmonths = 0
    
    for i, x in enumerate(xin[:,j]):
        iday = datevar[i].day
        if (iday != old_day) or (i == len(datevar) -1):
            # for daily files: need date of day before
            old_date = ('{:3d} {:4d}' + '{:3d}'*2).format(sum_days, 
                wanted_year, old_month, old_day)
            old_day = iday
            # values of day that just ended
            # special for last time step
            if i == len(datevar) -1:
                day_vals = xin[old_end:i+1,j]
            else:
                day_vals = xin[old_end:i,j]
            old_end = i
            dayline = (' {: 11.5e}'*(24/freq)).format(*day_vals)
            vals_file.write('{}{}\n'.format(old_date, dayline))
            day_mean_i = np.mean(day_vals)
            daily_means.append(day_mean_i)
            dmean_file.write('{:} {: 11.5e}\n'.format(old_date,
                day_mean_i))
            day_max_i = max(day_vals)
            daily_maxs.append(day_max_i)
            dmax_file.write('{:} {: 11.5e}\n'.format(old_date,
                day_max_i))
            # last day
            if i == len(datevar) -1:
                days_dates.append(datevar[i] - dt.timedelta(hours=23))
            else:
                days_dates.append(datevar[i] - dt.timedelta(hours=24))
                sum_days += 1

        # check if month is finished
        imonth = datevar[i].month
        if (imonth != old_month) or ((i == len(datevar) -1) and 
            (iday == days_last_month)):
            # only do monthly means if there are more than 20 days
            if (i - old_end_month)*freq / 24. > 20.:
                nmonths += 1
                which_months.append(old_month)
                # again special for very last day
                if (i == len(datevar) -1) and (iday == days_last_month):
                    month_vals = xin[old_end_month:i+1,j]
                else:
                    month_vals = xin[old_end_month:i,j]
                month_mean_i = np.mean(month_vals)
                monthly_means.append(month_mean_i)
                mmean_file.write('{:4d} {:2d} {: 11.5e}\n'.format(
                    wanted_year, old_month, month_mean_i))
            else:
                print ('no monthly statistics for {:02d} because there are '+
                'not more than 20 days of this month').format(old_month)
            old_end_month = i
            old_month = imonth
    

        today = datevar[i]
        hrly_file.write(('{:3d} {:4d}' + '{:3d}'*3 + ' {: 11.5e}\n').format(
            sum_days, today.year, today.month, today.day, today.hour, x))


    # at the end of the first site
    if j == 0:
        all_days_dates = days_dates
    else:
        if days_dates != all_days_dates:
            sys.exit(('something is wrong with the daily dates of the ' +
                'sites {} and {}').format(all_sites[0], site))

    all_day_means.append(daily_means)
    all_day_maxs.append(daily_maxs)
    all_month_means.append(monthly_means)

    vals_file.close()
    hrly_file.close()
    dmean_file.close()
    dmax_file.close()
    mmean_file.close()

# make lists of lists np.arrays with dimensions (time, nsites)
all_day_means = np.asarray(all_day_means).transpose()
all_day_maxs = np.asarray(all_day_maxs).transpose()
all_month_means = np.asarray(all_month_means).transpose()

# plotting
# ===========================================================================
# in case the type of the plot hasn't been given in the beginning, ask for it
#  now
if which_plot is None:
    print """
Now all the files have been written.
In case you do want the data to be plotted: choose between
[1] plot all data (every {} hour(s))
[2] plot daily means
[3] plot daily maxima
[4] plot monthly means
Press any other key for no plot.
    """.format(freq)
    which_plot_i = raw_input('Your choice: ')
    if (not which_plot_i.isdigit()) or (int(which_plot_i) not in [1, 2, 3, 4]):
        sys.exit('no plot')

    if int(which_plot_i) == 1:
        which_plot = 'all'
    elif int(which_plot_i) == 2:
        which_plot = 'daymeans'
    elif int(which_plot_i) == 3:
        which_plot = 'daymaxs'
    elif int(which_plot_i) == 4:
        which_plot = 'monmeans'

if which_time is None:
    print """
Do you want the data to be plotted for the entire time [1] or 
only for a specific time period [2]?
"""
    which_time_i = raw_input('Your choice: ')
    if (which_time_i.isdigit()) and (int(which_time_i) == 1):
        which_time = 1
    elif (which_time_i.isdigit()) and (int(which_time_i) == 2):
        which_time = 2

first_time = True
# case that a specific time range is to be plotted
if which_time == 2:
    # hourly plot
    if which_plot == 'all':
        # date_beg is measure if some time range has been given in the
        #  beginning already. Otherwise ask for it
        if len(date_beg) == 0:
            print 'Ok, enter the dates separated by a "-".'
            print 'Format yyyymmddhh, e.g. 2016093012-2016100721'
        while True:
            if first_time and len(date_beg) > 0:
            # date_beg will now also be checked if it makes sense.
                date_str = date_beg
                first_time = False
            else:
                date_str = raw_input('Your date: ')
            # check if date is in the correct format, later also if it fits
            #  the time range of the data
            if ((len(date_str) == 21) and (date_str[10] == '-') and
                (0 <= int(date_str[8:10]) <= 23) and 
                (0 <= int(date_str[19:]) <= 23)):
                start = dt.datetime.strptime(date_str[:10], '%Y%m%d%H')
                if (start - datevar[0]).total_seconds() < 0.:
                    print ('Starting date given is before first date of ' +
                        'the data: {}. Will simply start there.').format(
                        datevar[0].strftime('%Y%m%d%H'))
                    time.sleep(1)
                    i_start = 0
                else:
                    i_start = next((y for y, d in enumerate(datevar) 
                        if d==start), -10)
                if i_start == -10:
                    print ('Wrong start date? First date is {}, last ' +
                        'date is {}').format(datevar[0].strftime('%Y%m%d%H'),
                        datevar[-1].strftime('%Y%m%d%H'))
                else:
                    end = dt.datetime.strptime(date_str[11:], '%Y%m%d%H')
                    if (end - datevar[-1]).total_seconds() > 0.:
                        print ('Ending date given is after last date of ' +
                            'the data: {}. Will simply end there.').format(
                            datevar[-1].strftime('%Y%m%d%H'))
                        time.sleep(1)
                        i_end = len(datevar) -1
                    else:
                        i_end = next((y for y, d in enumerate(datevar)
                            if d == end), -10)
                    if (i_end == -10) or (i_end < i_start):
                        print ('Wrong end date? First date is {}, last ' +
                            'date is {}').format(datevar[0].strftime(
                            '%Y%m%d%H'),
                            datevar[-1].strftime('%Y%m%d%H'))
                    else:
                        i_end = i_end + 1
                        xplot = xin[i_start:i_end,:]
                        date_plot1 = datevar[i_start:i_end]
                        date_plot = dates.date2num(date_plot1)
                        # for title and filename
                        plot_type = ('all values (every {} hour(s)) {} - ' +
                            '{}').format(freq,
                            datevar[i_start].strftime('%Y-%m-%d %H:%M'),
                            datevar[i_end-1].strftime('%Y-%m-%d %H:%M'))
                        plot_fn = 'all_{}_{}'.format(
                            datevar[i_start].strftime('%Y%m%d%H'),
                            datevar[i_end-1].strftime('%Y%m%d%H'))
                        break
            else:
                print 'Format of date was not correct, try again!'
                if len(date_beg) > 0:
                    print 'Format yyyymmddhh, e.g. 2016093012-2016100721'

    # daily mean or max
    elif which_plot in ['daymeans', 'daymaxs']:
        # date_beg is measure if some time range has been given in the
        #  beginning already. Otherwise ask for it
        if len(date_beg) == 0:
            print 'Ok, enter the dates separated by a "-".'
            print 'Format yyyymmdd, e.g. 20160930-20161007'
        while True:
            if first_time and len(date_beg) > 0:
             # date_beg will now also be checked if it makes sense.
                date_str = date_beg
                first_time = False
            else:
                date_str = raw_input('Your date: ')
            # check if date is in the correct format, later also if it fits
            #  the time range of the data
            if (len(date_str) == 17) and (date_str[8] == '-'):
                start = dt.datetime.strptime(date_str[:8]+'00', '%Y%m%d%H')
                if (start - days_dates[0]).total_seconds() < 0.:
                    print ('Starting day given is before first day of ' +
                        'the data: {}. Will simply start there.').format(
                        days_dates[0].strftime('%Y%m%d'))
                    time.sleep(1)
                    i_start = 0
                else:
                    i_start = next((y for y, d in enumerate(days_dates) 
                        if d==start), -10)
                if i_start == -10:
                    print ('Wrong start date? First date is {}, last date ' +
                    'is {}').format(days_dates[0].strftime('%Y%m%d'),
                        days_dates[-1].strftime('%Y%m%d'))
                else:
                    end = dt.datetime.strptime(date_str[9:]+'00', '%Y%m%d%H')
                    if (end - days_dates[-1]).total_seconds() > 0.:
                        print ('Ending day given is after last day of ' +
                            'the data: {}. Will simply end there.').format(
                            days_dates[-1].strftime('%Y%m%d'))
                        time.sleep(1)
                        i_end = len(days_dates) -1
                    else:
                        i_end = next((y for y, d in enumerate(days_dates)
                            if d == end), -10)
                    if (i_end == -10) or (i_end < i_start):
                        print ('Wrong end date? First date is {}, last date '+
                            'is {}').format(days_dates[0].strftime('%Y%m%d'),
                            days_dates[-1].strftime('%Y%m%d'))
                    else:
                        i_end = i_end + 1
                        # for title and filename
                        if which_plot == 'daymeans':
                            xplot = all_day_means[i_start:i_end,:]
                            plot_ttype = 'means'
                            plot_ftype = 'means'
                        else:
                            xplot = all_day_maxs[i_start:i_end,:]
                            plot_ttype = 'maxima'
                            plot_ftype = 'maxs'
                        plot_type = 'daily {} {} - {}'.format(plot_ttype,
                            days_dates[i_start].strftime('%Y-%m-%d'),
                            days_dates[i_end-1].strftime('%Y-%m-%d'))
                        plot_fn = 'daily_{}_{}_{}'.format(plot_ftype,
                            days_dates[i_start].strftime('%Y%m%d'),
                            days_dates[i_end-1].strftime('%Y%m%d'))
                        date_plot = dates.date2num(days_dates[i_start:i_end])
                        break
            else:
                print 'Format of date was not correct, try again!'
                if len(date_beg) > 0:
                    print 'Format yyyymmdd, e.g. 20160930-20161007'

    # monthly means
    elif which_plot == 'monmeans':
        # date_beg is measure if some time range has been given in the
        #  beginning already. Otherwise ask for it
        if len(date_beg) == 0:
            print 'Ok, enter the dates separated by a "-".'
            print 'Format yyyymm, e.g. 201609-201610'
        while True:
            if first_time and len(date_beg) > 0:
            # date_beg will now also be checked if it makes sense.
                date_str = date_beg
                first_time = False
            else:
                date_str = raw_input('Your date: ')
            # check if date is in the correct format, later also if it fits
            #  the time range of the data
            if (len(date_str) == 13) and (date_str[6] == '-'):
                start = int(date_str[4:6])
                if ((int(date_str[:4]) != wanted_year) or
                    (start > which_months[-1])):
                    print ('Wrong start date? First month is {:4d}{:02d}, ' +
                        'last month is {:4d}{:02d}').format(wanted_year,
                        which_months[0],wanted_year, which_months[-1])
                else:
                    i_start = start - which_months[0]
                    if i_start < 0:
                        print ('Starting month given is before first month' +
                            ' of the data: {:02d}. Will simply start there.'
                            ).format(which_months[0])
                        time.sleep(1)
                        i_start = 0
                        start = which_months[0]
                    end = int(date_str[11:])# + 1
                    if ((end < which_months[0]) or 
                        ((end <= start) and (int(date_str[7:11]) <= wanted_year)) or 
                        (int(date_str[7:11]) < wanted_year)):
                        print ('Wrong end date? First month is ' +
                            '{:4d}{:02d}, last month is {:4d}{:02d}').format(
                            wanted_year,which_months[0],wanted_year, 
                            which_months[-1])
                    else:
                        i_end = end - which_months[0] + 1
                        if end > which_months[-1]:
                            print ('Ending month given is after last ' +
                                'month of the data: {:02d}. Will simply ' +
                                'end there.').format(which_months[-1])
                            time.sleep(1)
                            i_end = nmonths + 1
                            end = which_months[-1] #-1
                        xplot = all_month_means[i_start:i_end,:]
                        date_plot = dates.date2num([dt.datetime.strptime(
                            '{:4d}{:02d}01'.format(wanted_year,x),
                            '%Y%m%d') for x in range(start,end+1)])
                        # for title and filename
                        plot_type = \
                            'monthly means {:4d}{:02d} - {:4d}{:02d}'.format(
                            wanted_year, start, wanted_year, end)
                        plot_fn = ('monthly_means_{:4d}{:02d}_' +
                            '{:4d}{:02d}').format(
                            wanted_year, start, wanted_year, end)
                        break
            else:
                print 'Format of date was not correct, try again!'
                if len(date_beg) > 0:
                    print 'Format yyyymm, e.g. 201609-201610'

else:
    print 'Ok, will use all time steps.'
    # all (hourly) values
    if which_plot == 'all':
        # in case that there are data for more than one year
        xplot = xin
        date_plot = dates.date2num(datevar)
        # for title and filename
        plot_type = 'all values (every {} hour(s)) {} - {}'.format(
            freq, datevar[0].strftime('%Y-%m-%d %H:%M'),
            datevar[-1].strftime('%Y-%m-%d %H:%M'))
        plot_fn = 'all_{}_{}'.format(
            datevar[0].strftime('%Y%m%d%H'),
            datevar[-1].strftime('%Y%m%d%H'))
    # daily means or maxima
    elif which_plot in ['daymeans', 'daymaxs']:
        if which_plot == 'daymeans':
            xplot = all_day_means
            # for title and filename
            plot_kind = 'means'
            plot_kind_file = 'means'
        else:
            xplot = all_day_maxs
            # for title and filename
            plot_kind = 'maxima'
            plot_kind_file = 'maxs'
        date_plot = dates.date2num(days_dates)
        # for title and filename
        plot_type = 'daily {} {} - {}'.format(plot_kind,
            days_dates[0].strftime('%Y-%m-%d'),
            days_dates[-1].strftime('%Y-%m-%d'))
        plot_fn = 'daily_{}_{}_{}'.format(plot_kind_file,
            days_dates[0].strftime('%Y%m%d'),
            days_dates[-1].strftime('%Y%m%d'))
    # monthly means
    elif which_plot == 'monmeans':
        xplot = all_month_means
        date_plot = dates.date2num([dt.datetime.strptime(
            '{:4d}{:02d}01'.format(wanted_year,x),
            '%Y%m%d') for x in which_months])
        # for title and filename
        plot_type = 'monthly means {:4d}{:02d} - {:4d}{:02d}'.format(
            wanted_year, which_months[0], wanted_year, which_months[-1])
        plot_fn = 'monthly_means_{:4d}{:02d}_{:4d}{:02d}'.format(
            wanted_year, which_months[0], wanted_year, which_months[-1])


# now comes the real plotting
plt.figure(figsize = (10,5))
ax1 = plt.subplot(111)
colors = ["aqua", "blue", "fuchsia", "gray", "green", "lime", "maroon",
    "navy", "olive", "purple", "red", "silver", "teal", "yellow", "black"]

if which_plot == 'monmeans':
    # tick at every 15th
    locat = dates.MonthLocator(bymonthday = 1)
    # matplotlib date format object
    hfmt = dates.DateFormatter('%b')
    rotat = 'horizontal'
    marks = 'o'
    align_lab = 'center'
else:
    if date_plot[-1] - date_plot[0] <= 153:
        # tick at every 1st and 15th
        locat = dates.MonthLocator(bymonthday = (1,15))
    else:
        # tick at every 1st
        locat = dates.MonthLocator(bymonthday = 1)
    # matplotlib date format object
    hfmt = dates.DateFormatter('%Y-%m-%d')
    rotat='45'
    marks = None
    align_lab = 'right'

for i, site in enumerate(all_sites):
    ax1.plot(date_plot, xplot[:,i], linewidth = 1, color = colors[i],
        label = sitename[site], marker = marks)
ax1.legend(loc='best', frameon = False)


# beautify x-labels
ax1.xaxis.set_major_formatter(hfmt)
ax1.xaxis.set_major_locator(locat)
#DS ha idea from http://stackoverflow.com/questions/14852821/aligning-rotated-xticklabels-with-their-respective-xticks
plt.xticks(rotation=rotat, ha=align_lab)
ax1.set_xlim(date_plot[0],date_plot[-1])

plt.ylabel('{} in {}'.format(wpoll, unit))
# in case the variable is not temperature and the max is above 1000 or below
#  0.1 (-> scilimits), use scientific notation for y-axis
if (unit != 'K'):
    ax1.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0.1, 1000.))

# shrink plot 
plt.subplots_adjust(top = 0.94)
plt.subplots_adjust(bottom = 0.15)

fname_all_sites = '_'.join([sitename[i] for i in all_sites])
tname_all_sites = fname_all_sites.replace('_', ', ')
if len(tname_all_sites) > 100:
    tname_all_sites = tname_all_sites[:55] + '\n' + tname_all_sites[55:]

plt.title('{} for site(s) {}\n{}'.format(wpoll, tname_all_sites, plot_type))
plt.tight_layout()

out_plot = '{}_{}_{}.png'.format(fname_all_sites, wpoll, plot_fn)
plt.savefig(out_plot, bbox_inches = 'tight')






