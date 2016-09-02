#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib import dates
import matplotlib.cm as cm
import netCDF4 as cdf
import numpy as np
import sys
import time
import datetime as dt  # Python standard library datetime  module
import pdb
import calendar # for days in months
import argparse # argument passing to script


# author: Hannah Imhof, Chalmers, 2016, 
# loosely based upon David Simpson's Rd_csvsondes.f90
# but for netcdf files and adding plotting functionality


#  Reads from input file (e.g. model file sondes_2007.nc say), and provides
#  various outputs for a chosen pollutant.
#
#   SONDES.vals - provides output for each time step (by row), with here
#                 NK=20 columns representing 20 height levels.
#   SONDES.mmean - monthly means for each k level
#   SONDES.dmean  - daily mean for each k level
#   SONDES.dmax  - daily max for each k level
#   SONDES.sname  - site name


# to use this script from ipython, call it like this:
# os.system('EMEP/utils/Rd_ncsondes.py -i /home/ihannah/Documents/EMEP/Daves_output/sondes_2012.nc')

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--ifile', help = 'Input file sondes_yyyy.nc',
    required = True)
parser.add_argument('-s', '--sitename', help = 'Name of site')
parser.add_argument('-v', '--varname', help = 'Name of variable in nc file')
parser.add_argument('-k', '--klevel', help = ('Level of additional output. ' +
    'k=1 is surface'))
parser.add_argument('-mi', '--minconc', help = ('Minimum of concentration ' +
    'in plot'), type = float)
parser.add_argument('-ma', '--maxconc', help = ('Maximum of concentration ' +
    'in plot'), type = float)
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
print 'Input file: {}'.format(infile.strip())

# check for input for all the optional arguments
if args.sitename is None:
    wsite = ''
else:
    wsite = args.sitename
    print 'Site wanted: {}'.format(wsite)

if args.varname is None:
    wpoll = ''
else:
    wpoll = args.varname
    print 'Poll wanted: {}'.format(wpoll)

if args.klevel is None:
    klevel = ''
    wanted_k = 0
else:
    klevel = args.klevel
    print 'klevel wanted: {} (k = 1 is surface)'.format(klevel)
    wanted_k = int(klevel)

if args.minconc is not None: # None works fine in plot
    print 'Minimum concentration wanted: {}'.format(args.minconc)

if args.maxconc is not None:
    print 'Minimum concentration wanted: {}'.format(args.maxconc)

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


# has to be put on 0 in beginning
last_day = 0

# read input file
inCDF = cdf.Dataset(infile,'r',format='NETCDF4')

nsites = getattr(inCDF, 'Number_of_stations_defined')
# site names are saved in dimension (nsites,36) with single letters as
#  elements of a list
sitename = [str(''.join(inCDF.variables['station_name'][a,:]).strip()) 
    for a in range(nsites)]

freq = int(getattr(inCDF, 'Number_of_hours_bewtween_outputs'))

# wsite: name of site
# wanted_site: index of that site
if (len(wsite) > 0) and (wsite in sitename):
    wanted_site = sitename.index(wsite)
# ask which site to use if it hasn't been given at the beginning
else:
    to_print = ['{:4} {}'.format(i, x) for i, x in 
        enumerate(sitename)]
    print '\n'.join(to_print)
    wanted_site = int(raw_input('Choose site\n'))
    if (wanted_site < 0) or (wanted_site >= nsites):
        print '\n'.join(to_print)
        raise ValueError('NO SITE OF THAT NAME, START AGAIN!')
    wsite = sitename[wanted_site]
print 'NSITES {}, FREQ {}, WANTED_SITEXX {} {}, Chose site {}'.format(
    nsites, freq, wanted_site, wsite, 
    sitename[wanted_site])

# ===========================================================================
# now looking at variables
wanted_spec = -1

nspec = len(inCDF.variables)
specname = [var for var in inCDF.variables]

# wpoll: name of variable (species, pollutant)
# wanted_spec: index
if (len(wpoll) > 0) and (wpoll in specname):
    wanted_spec = specname.index(wpoll)
# ask which species to use if it hasn't been given at the beginning
else:
    to_print = ['{:4} {}'.format(i, x) for i, x in enumerate(specname)]
    print '\n'.join(to_print)
    wanted_spec = int(raw_input('Choose species\n'))
    if (wanted_spec < 0) or (wanted_spec >= nspec):
        print '\n'.join(to_print)
        raise ValueError('NO POLL OF THAT NAME, START AGAIN!')
    wpoll = specname[wanted_spec]

# this is the selected variable for the selected site
# Will be transposed when plotting
xin = inCDF.variables[wpoll][:,:,wanted_site]

unit = inCDF.variables[wpoll].units

print 'Site {} ({}), Poll {} ({})'.format(wanted_site, sitename[wanted_site],
    wanted_spec, specname[wanted_spec])

# Output files
# ===========================================================================
prefix = 'SONDES_{}_{}'.format(wsite,wpoll)
vals_file = open('{}.vals'.format(prefix), 'w')
snam_file = open('{}.sname'.format(prefix), 'w')
dmax_file = open('{}.dmax'.format(prefix), 'w')
dmean_file = open('{}.dmean'.format(prefix), 'w')
mmean_file = open('{}.mmean'.format(prefix), 'w')
if wanted_k > 0: # take care: will be pythonic index +1!!
    kvals_file = open('{}.k{}vals'.format(prefix,klevel), 'w')
# ===========================================================================

snam_file.write(wsite)
snam_file.close()


# looking at data now
# ===========================================================================
nk = np.shape(xin)[1]
print '{} levels'.format(nk)

# Date
nctime = inCDF.variables['time'][:]
unit_time = inCDF.variables['time'].units # get unit 'days since ...'
try :
    cal_time = inCDF.variables['time'].calendar
except AttributeError : # Attribute doesn't exist
    cal_time = u"gregorian" # or standard
datevar = cdf.num2date(nctime, units = unit_time,
        calendar = cal_time)

wanted_year = datevar[0].year

# Height
z_mid = inCDF.variables['z_mid'][:,:,wanted_site]
z_midmean = z_mid.mean(axis=0)

# as in Rd_csvsondes.f90: copy first timestep in data to 00:00
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

# go through variable's data and write output files
# ===========================================================================
old_day = datevar[0].day
sum_days = 0
old_end = 0
daily_maxs = []
daily_means = []
# keep track of the days
days_dates = []
# also for monthly data
old_month = datevar[0].month
old_end_month = 0
which_months = []
monthly_means = []
nmonths = 0
# additionally need number of days of last month
cal = calendar.Calendar()
if datevar[-1].year > wanted_year:
    last_month = cal.itermonthdays(wanted_year,12)
else:
    last_month = cal.itermonthdays(wanted_year,datevar[-1].month)
days_last_month = []
[days_last_month.append(x) for x in last_month]
days_last_month = max(days_last_month)

# write headers into output files
head_vals = '#day  mm  dd  hh ' + (' {:11}'*nk).format(
    *['k_{}_{}m'.format(k+1, int(round(z_midmean[k]))) for k in range(nk)])
vals_file.write(head_vals + '\n')
head_daily = '#day  mm  dd' + (' {:11}'*nk).format(
    *['k_{}_{}m'.format(k+1, int(round(z_midmean[k]))) for k in range(nk)])
dmax_file.write(head_daily + '\n')
dmean_file.write(head_daily + '\n')

for i, d in enumerate(datevar):
    iday = d.day
    # string for current date
    day_date = ('{:4d}'*3).format(sum_days, d.month, iday)
    # a new day has started or end of data
    if (iday != old_day) or (i == len(datevar) -1):
        # for daily files, need date of day before
        old_date = ('{:4d}'*3).format(sum_days, old_month, old_day)
        old_day = iday
        sum_days += 1
        # values of day that just ended
        day_vals = xin[old_end:i,:]
        old_end = i
        k_daymeans = np.mean(day_vals,0)
        daily_means.append(k_daymeans)
        out_daymeans = ('{:11.5e} '*nk).format(*k_daymeans)
        dmean_file.write('{} {}\n'.format(old_date, out_daymeans))

        k_daymaxs = np.max(day_vals,0)
        daily_maxs.append(k_daymaxs)
        out_daymaxs = ('{:11.5e} '*nk).format(*k_daymaxs)
        dmax_file.write('{} {}\n'.format(old_date, out_daymaxs))
        # last day
        if i == len(datevar) -1:
            days_dates.append(d - dt.timedelta(hours=23))
        else:
            days_dates.append(d - dt.timedelta(hours=24))
    # check if month is finished
    imonth = d.month
    if (imonth != old_month) or ((i == len(datevar) -1) and 
        (iday == days_last_month)):
        # only do monthly means if there are more than 20 days
        if (i - old_end_month)*freq / 24. > 20.:
            nmonths += 1
            which_months.append(old_month)
            month_vals = xin[old_end_month:i,]
            means_levels = []
            for k in range(nk):
                means_levels.append(np.mean(month_vals[:,k]))
            monthly_means.append(means_levels)
        else:
            print ('no monthly statistics for {:02d} because there are not more '+
            'than 20 days of this month').format(old_month)
        old_end_month = i
        old_month = imonth

    if d.year != wanted_year:
        print 'year {:4d} is done'.format(wanted_year)
        last_day = i-1
        break
    out_kvals = ('{:11.5e} '*nk).format(*xin[i,:])
    vals_file.write('{} {:3d} {}\n'.format(day_date, d.hour, out_kvals))
    if wanted_k > 0:
        if i == 0:
            head_k = '#day  mm  dd  hh k_{}_{}m'.format(wanted_k, int(round(
                z_midmean[wanted_k-1])))
            kvals_file.write(head_k + '\n')
        kvals_file.write('{} {:3d} {:11.5e}\n'.format(day_date, d.hour, 
            xin[i,wanted_k-1]))
# in case that there is not more than one year in the file
if last_day == 0:
    last_day = i

dmean_file.close()
dmax_file.close()
vals_file.close()
if wanted_k > 0: 
    kvals_file.close()

# write monthly mean file
head_mmean = '  k' + (' {:11}'*nmonths).format(
    *['mm_{}'.format(mm) for mm in which_months])
mmean_file.write(head_mmean + '\n')
# thus, matrix can be transposed and it's a real 2D matrix
monthly_means = np.asarray(monthly_means)
for k in range(nk):
    outmonth = ('{:3d}' + ' {:11.5e}'*nmonths).format(k+1, 
        *monthly_means[:,k])
    mmean_file.write(outmonth + '\n')
mmean_file.close()

# same: will be real 2D matrices
daily_means = np.asarray(daily_means)
daily_maxs = np.asarray(daily_maxs)



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
                    i_start = [0]
                else:
                    i_start = [y for y, d in enumerate(datevar)
                        if d == start]
                if len(i_start) != 1:
                    print ('Wrong start date? First date is {}, last ' +
                        'date is {}').format(datevar[0].strftime('%Y%m%d%H'),
                        datevar[last_day].strftime('%Y%m%d%H'))
                else:
                    i_start = i_start[0]
                    end = dt.datetime.strptime(date_str[11:], '%Y%m%d%H')
                    if (end - datevar[last_day]).total_seconds() > 0.:
                        print ('Ending date given is after last date of ' +
                            'the data: {}. Will simply end there.').format(
                            datevar[last_day].strftime('%Y%m%d%H'))
                        time.sleep(1)
                        i_end = [last_day]
                    else:
                        i_end = [y for y, d in enumerate(datevar)
                            if d == end]
                    if (len(i_end) != 1) or (i_end[0] < i_start):
                        print ('Wrong end date? First date is {}, last ' +
                            'date is {}').format(datevar[0].strftime(
                            '%Y%m%d%H'),
                            datevar[last_day].strftime('%Y%m%d%H'))
                    else:
                        i_end = i_end[0] + 1
                        xplot = xin[i_start:i_end,:]
                        date_plot = dates.date2num(datevar[i_start:i_end])
                        # for the dates to be correct on the x-axis, another
                        #  hour has to be added at the end
                        # End of the year:
                        if ((i_end >= len(datevar)) and (datevar[i_end-1].
                            strftime('%m%d%H') == '123123')):
                            date_plot = np.append(date_plot, dates.date2num(
                            dt.datetime.strptime('{:4d}010100'.format(
                            wanted_year+1), '%Y%m%d%H')))
                        else:
                            date_plot = np.append(date_plot, dates.date2num(
                            datevar[i_end-1] + dt.timedelta(hours = 1)))
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
                    i_start = [0]
                else:
                    i_start = [y for y, d in enumerate(days_dates) 
                        if d == start]
                if len(i_start) != 1:
                    print 'Wrong start date?'
                else:
                    i_start = i_start[0]
                    end = dt.datetime.strptime(date_str[9:]+'00', '%Y%m%d%H')
                    if (end - days_dates[-1]).total_seconds() > 0.:
                        print ('Ending day given is after last day of ' +
                            'the data: {}. Will simply end there.').format(
                            days_dates[-1].strftime('%Y%m%d'))
                        time.sleep(1)
                        i_end = [len(days_dates)-1]
                    else:
                        i_end = [y for y, d in enumerate(days_dates)
                            if d == end]
                    if len(i_end) != 1:
                        print 'Wrong end date?'
                    else:
                        i_end = i_end[0] + 1
                        # for title and filename
                        if which_plot == 'daymeans':
                            xplot = daily_means[i_start:i_end,:]
                            plot_ttype = 'means'
                            plot_ftype = 'means'
                        else:
                            xplot = daily_maxs[i_start:i_end,:]
                            plot_ttype = 'maxima'
                            plot_ftype = 'maxs'
                        plot_type = 'daily {} {} - {}'.format(plot_ttype,
                            days_dates[i_start].strftime('%Y-%m-%d'),
                            days_dates[i_end-1].strftime('%Y-%m-%d'))
                        plot_fn = 'daily_{}_{}_{}'.format(plot_ftype,
                            days_dates[i_start].strftime('%Y%m%d'),
                            days_dates[i_end-1].strftime('%Y%m%d'))
                        date_plot = dates.date2num(days_dates[i_start:i_end])
                        # for the dates to be correct on the x-axis, another
                        #  day has to be added at the end
                        # End of the year:
                        if ((i_end >= len(days_dates)) and
                            (days_dates[i_end-1].strftime('%m%d') == '1231')):
                            date_plot = np.append(date_plot, dates.date2num(
                            dt.datetime.strptime('{:4d}010100'.format(
                            wanted_year+1), '%Y%m%d%H')))
                        else:
                            date_plot = np.append(date_plot, dates.date2num(
                            days_dates[i_end-1] + dt.timedelta(days = 1)))
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
                        (int(date_str[7:11]) != wanted_year)):
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
                        xplot = monthly_means[i_start:i_end,:]
                        date_plot = dates.date2num([dt.datetime.strptime(
                            '{:4d}{:02d}01'.format(wanted_year,x),
                            '%Y%m%d') for x in range(start,end+1)])
                        # for the dates to be correct on the x-axis, another
                        #  month has to be added at the end
                        # End of the year:
                        if end == 12:
                            date_plot = np.append(date_plot, dates.date2num(
                            dt.datetime.strptime('{:4d}0101'.format(
                            wanted_year+1), '%Y%m%d')))
                        else:
                            date_plot = np.append(date_plot, dates.date2num(
                            dt.datetime.strptime('{:4d}{:02d}01'.format(
                            wanted_year, end+1), '%Y%m%d')))
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
        # pdb.set_trace()
        xplot = xin[:last_day+1,:]
        date_plot = dates.date2num(datevar[:last_day+1])
        # for the dates to be correct on the x-axis, another hour has to be
        #  added at the end
        # End of the year:
        if datevar[-1].strftime('%m%d%H') == '123123':
            date_plot = np.append(date_plot, dates.date2num(
            dt.datetime.strptime('{:4d}010100'.format(
            wanted_year+1), '%Y%m%d%H')))
        else:
            date_plot = np.append(date_plot, dates.date2num(
            datevar[-1] + dt.timedelta(hours = 1)))
        # for title and filename
        plot_type = 'all values (every {} hour(s)) {} - {}'.format(
            freq, datevar[0].strftime('%Y-%m-%d %H:%M'),
            datevar[last_day].strftime('%Y-%m-%d %H:%M'))
        plot_fn = 'all_{}_{}'.format(
            datevar[0].strftime('%Y%m%d%H'),
            datevar[last_day].strftime('%Y%m%d%H'))
    # daily means or maxima
    elif which_plot in ['daymeans', 'daymaxs']:
        if which_plot == 'daymeans':
            xplot = daily_means
            # for title and filename
            plot_kind = 'means'
            plot_kind_file = 'means'
        else:
            xplot = daily_maxs
            # for title and filename
            plot_kind = 'maxima'
            plot_kind_file = 'maxs'
        date_plot = dates.date2num(days_dates)
        # for the dates to be correct on the x-axis, another day has to be
        #  added at the end
        # End of the year:
        if days_dates[-1].strftime('%m%d') == '1231':
            date_plot = np.append(date_plot, dates.date2num(
            dt.datetime.strptime('{:4d}010100'.format(
            wanted_year+1), '%Y%m%d%H')))
        else:
            date_plot = np.append(date_plot, dates.date2num(
            days_dates[-1] + dt.timedelta(days = 1)))
        # for title and filename
        plot_type = 'daily {} {} - {}'.format(plot_kind,
            days_dates[0].strftime('%Y-%m-%d'),
            days_dates[-1].strftime('%Y-%m-%d'))
        plot_fn = 'daily_{}_{}_{}'.format(plot_kind_file,
            days_dates[0].strftime('%Y%m%d'),
            days_dates[-1].strftime('%Y%m%d'))
    # monthly means
    elif which_plot == 'monmeans':
        xplot = monthly_means
        # for the dates to be correct on the x-axis, another month has to be
        #  added at the end
        # End of the year:
        date_plot = dates.date2num([dt.datetime.strptime(
            '{:4d}{:02d}01'.format(wanted_year,x),
            '%Y%m%d') for x in which_months])
        if which_months[-1] == 12:
            date_plot = np.append(date_plot, dates.date2num(
                dt.datetime.strptime('{:4d}0101'.format(wanted_year+1),
                '%Y%m%d')))
        else:
            date_plot = np.append(date_plot, dates.date2num(
                dt.datetime.strptime('{:4d}{:02d}01'.format(wanted_year,
                which_months[-1]+1),'%Y%m%d')))
        # for title and filename
        plot_type = 'monthly means {:4d}{:02d} - {:4d}{:02d}'.format(
            wanted_year, which_months[0], wanted_year, which_months[-1])
        plot_fn = 'monthly_means_{:4d}{:02d}_{:4d}{:02d}'.format(
            wanted_year, which_months[0], wanted_year, which_months[-1])


# now comes the real plotting
plt.figure(figsize=(15, 7))
ax1 = plt.subplot(111)
v1 = ax1.imshow(xplot.transpose(), aspect='auto', origin = 'lower',
    extent = [date_plot[0], date_plot[-1], 0., (z_midmean[-1] + 
    0.5*(z_midmean[-1] - z_midmean[-2]))],
    cmap='RdBu_r', interpolation = 'none',
    vmin = args.minconc, vmax = args.maxconc)

if which_plot == 'monmeans':
    # tick at every 15th
    locat = dates.MonthLocator(bymonthday = 15)
    # matplotlib date format object
    hfmt = dates.DateFormatter('%b')
    rotat = 'horizontal'
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


# beautify x-labels
ax1.xaxis.set_major_formatter(hfmt)
ax1.xaxis.set_major_locator(locat)
#DS ha idea from http://stackoverflow.com/questions/14852821/aligning-rotated-xticklabels-with-their-respective-xticks
plt.xticks(rotation=rotat,ha='right')
# y-labels
upper_y = (z_midmean[-1] + 0.5*(z_midmean[-1] - z_midmean[-2]))
# to get the correctly matching z-labels on the y-axis (not to scale)
delta_y = upper_y / 20.
mid_y = [y*delta_y - 0.5*delta_y for y in range(1,21)]
ax1.set_ylim(0,upper_y)
plt.yticks(mid_y, ['{:d}'.format(int(round(x))) for x in z_midmean])
plt.ylabel('height in m')

# legend
cbar = plt.colorbar(v1, format='%.2e')
cbar.set_label('[{}] in {}'.format(wpoll, unit),size=14)

plt.title('{} for site {}, {}'.format(wpoll, wsite, plot_type))

plt.tight_layout()

out_plot = '{}_{}_sondes_{}.png'.format(wsite, wpoll, plot_fn)
plt.savefig(out_plot, bbox_inches = 'tight')
