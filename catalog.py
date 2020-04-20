#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Open Source EMEP/MSC-W model
simplified access to the source code, input data and benchmark results.
"""

import os
import sys
import hashlib
import tarfile
import shutil

_CONST = {
    'VERSION':"0.1.5",                          # script version
    'RELEASE': ['rv3', 'v201106', 'rv4_0', 'rv4_3', 'rv4_4', 'rv4_5', 'rv4_8',
               'rv4_10','rv4_15','rv4_17','rv4_32','rv4_33', 'rv4_34'],     # released model versions
    'METYEAR':[2005, 2008] + [year for year in range(2010, 2015+1)],   # released met-years
    'THREDDS':'http://thredds.met.no/thredds/fileServer/data/EMEP/OpenSource',
    'FTP':'ftp://ftp.met.no/projects/emep/OpenSource',
    'GIT':'https://github.com/metno/emep-ctm/',
    'DOC':'https://emep-ctm.readthedocs.io/',
    'RTD':'https://buildmedia.readthedocs.org/media/pdf/emep-ctm',
    'CSV':'/catalog.csv',                       # list all files from all releases
    'RAW':'https://raw.githubusercontent.com/metno/emep-ctm/'+
          'tools/catalog.csv',                  # catalog on the repo
    'TMPDIR':"./downloads",                     # temp path for downloads
    'DATADIR':'.'                               # base path for datasets
}

def parse_arguments(args):
    """Arguments from command line"""
    from optparse import OptionParser, OptionGroup, SUPPRESS_HELP

    usage = """usage: %prog [options]

Examples:

  Retrieve release dataset for revision REV ({REV})
    %prog -R REV

  Get Only the source code and user guide for revision REV
    %prog -R REV --source --docs

  Download meteorological input for YEAR ({MET})
    %prog -Y YEAR --meteo
""".format(
    REV="|".join(_CONST['RELEASE']),
    MET="|".join(["%d"%y for y in _CONST['METYEAR']]))

    parser = OptionParser(usage, version=_CONST['VERSION'])
    parser.set_defaults(verbose=1)
    parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose",
                      help="don't print status messages to stdout")
    parser.add_option("-v", "--verbose",
                      action="count", dest="verbose",
                      help="Increase verbosity")
    parser.add_option("--catalog", default=os.path.dirname(__file__)+_CONST['CSV'],
                      action="store", type="string", dest="catalog",
                      help="Override dataset cataloque path/file (default:%default)")

    group = OptionGroup(parser, "Release options", "Select release dataset")
    group.add_option("-R", "--revision",
                     type="string", metavar="REV",
                     action="append", dest="tag",
                     help="revision REV")
    group.add_option("-S", "--status",
                     type="int", metavar="YEAR",
                     action="append", dest="status",
                     help=SUPPRESS_HELP) # help="YEAR's status report")
    group.add_option("-Y", "--year",
                     type="int", metavar="YEAR",
                     action="append", dest="year",
                     help="Meteorological/run YEAR")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Dataset options", "Partial release dataset")
    group.add_option("-m", "--meteo", const="meteo",
                     action="append_const", dest="data",
                     help="get meteorology input")
    group.add_option("--met-domain", default=None,
                      action="store", type="string", dest="domain",
                      help="get only DOMAIN meteorology")
    group.add_option("-i", "--input", const="input",
                     action="append_const", dest="data",
                     help="get other input")
    group.add_option("-o", "--output", const="output",
                     action="append_const", dest="data",
                     help="get model benchmark")
    group.add_option("-s", "--source", const="source",
                     action="append_const", dest="data",
                     help="get source code for benchmark")
    group.add_option("-d", "--docs", const="docs",
                     action="append_const", dest="data",
                     help="get corresponding user guide")
    group.add_option("--extras",
                     action="store_true", dest="extras",
                     help="also get the extras, if any")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Download options", "")
    group.add_option("--yes", default=True,
                     action="store_false", dest="ask",
                     help="Don't ask before start downloading/unpacking")
    group.add_option("--outpath",
                     action="store", dest="outpath",
                     help="Override output (dataset) directory")
    group.add_option("--tmppath",
                     action="store", dest="tmppath",
                     help="Override temporary (download) directory")
    group.add_option("--cleanup", default=False,
                     action="store_true", dest="cleanup",
                     help="Remove ALL temporary (download) files")
    parser.add_option_group(group)

    opts, args = parser.parse_args(args)
    if all(getattr(opts, attr) is None for attr in ['tag', 'status', 'year']):
        opts.tag = [_CONST['RELEASE'][-1]]
    if opts.data is None:
        opts.data = ["meteo", "input", "output", "source", "docs"]
    if opts.extras:
        opts.data += ["extra"]
    if opts.outpath:
        _CONST['DATADIR'] = opts.outpath
    if opts.tmppath:
        _CONST['TMPDIR'] = opts.tmppath
    return opts, args

def user_consent(question, ask=True, default='yes'):
    """Ask user for confirmation"""

    # format question
    qq = {'y':"%s [Y]/n:", 'yes':"%s [Y]/n:", 'n':"%s y/[N]:", 'no':"%s y/[N]:"}
    try:
        question = qq[default.lower()]%question
    except KeyError:
        print("Unsupported option: user_consent(..,default='%s')"%default)

    # valid answers
    answer = {'y':True, 'yes':True, 'n':False, 'no':False}

    # ask until you get a valid answer
    while True:
        if not ask:                     # take the default
            print(question+default)
            data = default
        else:
            data = input(question)

        # use default when user replies ''
        data = data.lower() if len(data) else default.lower()
        if data in answer:
            return answer[data]

# Print iterations progress
# http://stackoverflow.com/a/34325723/2576368
def print_progress(iteration, total, prefix='', suffix='', decimals=2, length=100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : number of decimals in percent complete (Int)
        length   - Optional  : character length of bar (Int)
    """
    filled = int(round(length * iteration / float(total)))
    percents = round(100.00 * (iteration / float(total)), decimals)
    progress = '█' * filled + '-' * (length - filled)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, progress, percents, '%', suffix))
    sys.stdout.flush()
    if iteration == total:
        sys.stdout.write('\n')
        sys.stdout.flush()

def file_size(value):
    """
    Human readable file size (K/M/G/T bytes)
      derived from http://stackoverflow.com/a/1094933/2576368
    """
    num = max(float(value), 0)
    for unit in ['B', 'K', 'M', 'G']:
        if num < 1024.0:
            return "%3.1f%s"%(num, unit)
        num /= 1024.0
    return "%.1f%s"%(num, 'T')

class DataPoint(object):
    '''Info and retrieval/check methods'''

    def __init__(self, release, key, year, model, src,
                 dst='', byteSize=0, md5sum=None):
        '''Initialize object'''

        self.release = int(release)       # release date (YYYYMM)
        self.key = str(key)
        try:
            self.tag = {'other':'{KEY}{REL}',
                        'meteo':'{KEY}{YEAR}'}[self.key]
        except KeyError:
            self.tag = '{MOD}{KEY}'       # eg 'rv4_8source'
        try:
            self.year = int(year)         # met-year
        except ValueError:
            self.year = None
        try:
            self.model = str(model)       # model version (or meteo domain)
        except ValueError:
            self.model = None
        self.src = str(src)               # single source url/file
        self.dst = str(dst)               # path for uncompressed self.dst
        if self.dst == '':
            self.dst = '{TMPDIR}/{TAG}/'
        self.size = float(byteSize)       # self.src file size [bytes]
        self.md5sum = md5sum              # self.src checksum
        if self.md5sum:
            self.md5sum = str(self.md5sum)

        # replace keywords
        kwargs = _CONST
        kwargs.update(dict(REL=self.release, KEY=self.key,
                           YEAR=self.year, MOD=self.model))
        self.tag = self.tag.format_map(kwargs)
        kwargs.update(dict(TAG=self.tag))
        self.src = self.src.format_map(kwargs)
        self.dst = self.dst.format_map(kwargs)
        if os.path.basename(self.dst) == '':
            self.dst += os.path.basename(self.src)

    def __str__(self):
        return "%-14s %6s %s"%(self.tag, file_size(self.size), self.dst)
    def __repr__(self):
        return "%6s %s"%(file_size(self.size), os.path.basename(self.dst))

    def __hash__(self):
        '''find unique self.src occurrences'''
        return hash(repr(self.src))
    def __eq__(self, other):
        if isinstance(other, DataPoint):
            return self.src == other.src
        else:
            return False
    def __ne__(self, other):
        return not self.__eq__(other)

    def cleanup(self, verbose=True):
        """Remove (raw) downloads"""
        if not os.path.isfile(self.dst):
            return
        if verbose:
            print("%-8s %s"%('Cleanup', self))

        try:
            os.remove(self.dst)
            dname = os.path.dirname(self.dst)
            if dname != '':
                os.rmdir(dname)
        except OSError as error:
            if error.errno == 1:  # opperation not permitted (permissions?)
                pass
            elif error.errno == 39: # directory was not empty
                pass
            else:
                raise error
                

    def check(self, verbose=True, cleanup=False):
        """Check download against md5sum"""
        if not os.path.isfile(self.dst):
            return False
        if verbose:
            print("%-8s %s"%('Check', self))

        with open(self.dst, 'rb') as infile:
            if self.md5sum != hashlib.md5(infile.read()).hexdigest():
                if cleanup:    # remove broken file
                    self.cleanup(verbose)
                if verbose:
                    print("%-8s %s"%('  md5/=', self.md5sum))
                return False
            return True

    def download(self, verbose=True):
        """derived from http://stackoverflow.com/a/22776/2576368"""
        from urllib.request import urlopen
        # check if file exists/md5sum, remove file if fail md5sum
        if self.check(verbose > 2, cleanup=True):
            return

        if verbose:
            print("%-8s %s"%('Download', self))

        dname = os.path.dirname(self.dst)
        if (not os.path.isdir(dname))and(dname != ''):
            os.makedirs(dname)

        url = urlopen(self.src)
        with open(self.dst, 'wb') as outfile:
            block = 1024*8          #   8K
            big_file = self.size > 1024*1024*128 # 128M
            if verbose and big_file:
                n = 0
                ntot = self.size/block
                print_progress(n, ntot, length=50)
            while True:
                buff = url.read(block)
                if not buff:
                    break
                outfile.write(buff)
                if verbose and big_file:
                    n += 1
                    if n%128 == 0:  # 8K*128=1M
                        print_progress(n, ntot, length=50)
            if verbose and big_file and ntot%128 != 0:
                print_progress(ntot, ntot, length=50)

    def unpack(self, verbose=True, inspect=False):
        """Unpack download"""
        if not self.check(verbose > 2):
            return

        if tarfile.is_tarfile(self.dst):
            print("%-8s %s"%('Untar', self))
            infile = tarfile.open(self.dst, 'r')
            if user_consent('  See the contents first?', inspect, 'no'):
                print(infile.list(verbose=(verbose > 1)))
                if not user_consent('  Do you wish to proceed?', inspect):
                    print("OK, skipping file")
                    return
            dname = "%s/"%_CONST['DATADIR']
            if (not os.path.isdir(dname))and(dname != ''):
                os.makedirs(dname)
            try:
                infile.extractall(dname)
            except EOFError as error:
                print("  Failed unpack '%s':\n    %s."%(self.dst, error))
                if not user_consent('    Do you wish to continue?', inspect, 'yes'):
                    sys.exit(-1)
            infile.close()

        else:
            outfile = "%s/%s/%s"%(_CONST['DATADIR'], self.key, os.path.basename(self.dst))
            if verbose > 1:
                print("%-8s %s"%('Copy', outfile))
            dname = os.path.dirname(outfile)
            if (not os.path.isdir(dname))and(dname != ''):
                os.makedirs(dname)
            shutil.copyfile(self.dst, outfile)

class DataSet(object):
    '''Info and retrieval/check methods'''

    def __init__(self, tag, release, year, status, dataset=None, byteSize=None):
        '''Initialize object'''
        self.tag = str(tag)             # revision tag
        self.release = int(release)     # release date (YYYYMM)
        self.year = year                # met-year(meteo)/status-year(model)
        self.status = status            # Status report
        self.dataset = dataset          # {'input':input,'meteo':meteo,..}
        self.size = 0
        if byteSize:
            self.size = float(byteSize) or 0
        else:
            self.size = sum([x.size for x in self.dataset.values() if hasattr(x, "size")])

    def __str__(self):
        return "%-8s (meteo:%s, status:%s)"%(self.tag, self.year, self.status)
    def __repr__(self):
        return "%s: %s"%(self, {k: str(v) for k, v in self.dataset.items() if v})
#       return "%s: %s"%(self,self.dataset)

def read_catalog(filename, verbose=1):
    '''
    Returns releases read from catalog csv-file

    Definitions
      DataPoint(class): single remote file/tarfile
      DataSet(class):   dataPoints from a model release,
                        sorted into meteo|input|output|source|docs
      catalog(list):    all dataPoints
      index(dict):      catalog sorted into meteo|input|output|source|docs
      archive(list):    all DataSet
    '''
    from csv import reader as csvreader

    # download catalog if file not found
    if not os.path.isfile(filename):
        DataPoint(0, 'catalog', 0, '', _CONST['RAW'], filename).download(verbose)

    # catalog file should be present at this point
    with open(filename) as infile:
        reader = csvreader(infile, delimiter=',')
        next(reader)              # skip header
        catalog = []              # list all src files (1 file per row)
        rels, keys = set(), set() # unique DataPoint.release/.keys for indexing

        if verbose > 2:
            print("Reading %s"%filename)
        for row in reader:
            if row and "".join(row).strip():  # skip empty lines
                try:
                    catalog += [DataPoint(*row)]
                    rels.add(catalog[-1].release) # unique releases
                    keys.add(catalog[-1].key)     # unique keys (meteo,source,&c)
                    if verbose > 2:
                        print("  %s"%catalog[-1])
                except:
                    print("Failed to parse (%s): %r"%(filename, row))
                    print(DataPoint(*row))
                    raise
        if verbose > 1:
            print("%s read(srcs:%d)"%(filename, len(catalog)))

    index = dict.fromkeys(rels) # index[releases][keys:meteo,source,&c]
    if verbose > 2:
        print("Indexing")
    for r in rels:            # r:release
        index[r] = dict.fromkeys(keys)
        for k in keys:          # k:meteo,source,&c
            index[r][k] = [x for x in catalog if (x.release == r)and(x.key == k)]
            if verbose > 2:
                print("  index[%s][%s](srcs:%d)"%(r, k, len(index[r][k])))
    if verbose > 1:
        print("%s index[release:%d][sets:%d]"%(filename, len(rels), len(keys)))

    # repack index from dict({key:DataPoint,}) to dict(DataSet)
    if verbose > 2:
        print("Compiling releases")
    for r, v in index.items():
        # pack index[r] into a DataSet object
        #   before index[r]:{meteo:DataPoint,source:DataPoint,..}
        #   after  index[r]:DataSet (with search metadata)
        index[r] = DataSet(v['source'][0].model, r, v['meteo'][0].year if v['meteo'] else None,
                           v['source'][0].year, v)
        if verbose > 2:
            print("  index[%s]: %s"%(r, index[r]))
    if verbose > 1:
        print("%s index(release:%d)"%(filename, len(rels)))

    return index

def main(opts):
    """Command line function"""
    def get_datasets(catalog, attr, target):
        """Search catalog for tag|status|year DataSet"""
        if opts.verbose > 1:
            print("Searching %s(s):%s"%(attr, target))
        try:
            dataset = [v for _, v in catalog.items() if getattr(v, attr) in target]
            if len(dataset) == 0:
                raise IndexError
        except IndexError:
            print("No datasets found for --%s=%s"%(attr, target))
            sys.exit(-1)
        if opts.verbose > 1:
            for x in dataset:
                print("  Found %s"%x)
        return dataset
    def get_downloads(catalog, attr, target):
        """List downloads for tag|status|year DataSet"""
        downloads = []
        if opts.verbose > 1:
            print("Searching datasets:%s"%(opts.data))
        for ds in get_datasets(catalog, attr, target):
            try:
                ds = {key: ds.dataset[key] for key in opts.data}
                # only download meteo with matching --met-domain option
                if 'meteo' in ds and opts.domain:
                    ds['meteo'] = [ x for x in ds['meteo'] if x.model == opts.domain ]
            except KeyError:
                print("No datasets found for --%s=%s"%(attr, target))
                sys.exit(-1)
            if opts.verbose > 1:
                for key in ds:
                    print("  Found %-6s:%s"%(key, ds[key]))

            for key in ds:
                # do not try to download 0-size files
                aux = [x for x in ds[key] if x.size > 0]
                total = sum([x.size for x in aux])
                if(total>0):
                    downloads += aux
                if opts.verbose and aux:
                    print("Queue download: %6s %s"%(file_size(total), aux[0].tag))
        return downloads

    # list files to download
    catalog = read_catalog(opts.catalog, opts.verbose)
    downloads = [] # files to download
    for attr in ['tag', 'status', 'year']:
        target = getattr(opts, attr)
        if target is None:
            continue
        downloads += get_downloads(catalog, attr, target)

    # list file_sizes to download
    downloads = list(set(downloads))  # unique files, for single download and total size
    total = sum([x.size for x in downloads])
    print("Queue download: %6s %s"%(file_size(total), 'Total'))
    if not user_consent('Do you wish to proceed?', opts.ask):
        print("OK, bye")
        sys.exit(0)

    # download files
    for x in downloads:
        x.download(opts.verbose)
        x.unpack(opts.verbose, opts.ask)
        if opts.cleanup:
            x.cleanup(opts.verbose)

if __name__ == "__main__":
    main(parse_arguments(sys.argv[1:] or ['-h'])[0])
