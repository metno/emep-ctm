#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Open Source EMEP/MSC-W model
simplified access to the source code, input data and benchmark results.
"""

_CONST={
  'VERSION':"0.0.1",      # script version
  'LASTREL':"rv4.8",      # last model release
  'LASTMET':2013,         # last met-year released
  'FTP':"ftp://ftp.met.no/projects/emep/OpenSource",
  'CSV':'./catalog.csv',  # list all files from all releases
  'TMPDIR':"./downloads", # temp path for downloads
  'DATADIR':'.'           # base path for datasets
}

# Print iterations progress
# http://stackoverflow.com/a/34325723/2576368
def printProgress(iteration,total,prefix='',suffix='',decimals=2,barLength=100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    from sys import stdout
    filledLength    = int(round(barLength * iteration / float(total)))
    percents        = round(100.00 * (iteration / float(total)), decimals)
    bar             = '█' * filledLength + '-' * (barLength - filledLength)
    stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
    stdout.flush()
    if iteration == total:
      stdout.write('\n')
      stdout.flush()

class fileSize(int):
  """Human readable file size"""
  def __new__(cls, value):
    return int.__new__(cls,max(value,0))  # only positive values make sence
  def __repr__(self):
    return "%d%s"%(self*4,'b')            # show in bites
  def __str__(self):                      # show K/M/G/T bytes
    """ derived from http://stackoverflow.com/a/1094933/2576368"""
    num=float(self)
    for unit in ['B','K','M','G']:
      if num < 1024.0:
        return "%3.1f%s"%(num, unit)
      num /= 1024.0
    return "%.1f%s"%(num, 'T')

class dataPoint(object):
  '''Info and retrieval/check methods'''

  def __init__(self,release,key,year,model,src,dst='',byteSize=0,md5sum=None):
    '''Initialize object'''
    from os.path import basename
    
    self.release= int(release)        # release date (YYYYMM)
    self.key    = str(key)
    try: 
      self.tag  = {'other':'{KEY}{REL}',
                   'meteo':'{KEY}{YEAR}'}[self.key]
    except:
      self.tag  = '{MOD}{KEY}'      # eg 'rv4_8source'
    try:
      self.year = int(year)           # met-year
    except:
      self.year = None
    try:
      self.model= str(model)          # model version
    except:
      self.model= None
    self.src    = str(src)            # single source url/file
    self.dst    = str(dst)            # path for uncompressed self.dst
    if self.dst=='':
      self.dst= '{TMPDIR}/{TAG}/'
    self.size   = fileSize(byteSize)  # self.src file size [bytes]
    self.md5sum  = md5sum             # self.src checksum
    if self.md5sum:
      self.md5sum= str(self.md5sum)

    # replace keywords
    self.tag = self.tag.format(REL=self.release,KEY=self.key,YEAR=self.year,MOD=self.model)
    self.src = self.src.format(REL=self.release,TAG=self.tag,**_CONST)
    self.dst = self.dst.format(REL=self.release,TAG=self.tag,**_CONST)
    if basename(self.dst)=='':
      self.dst+=basename(self.src)

  def __str__(self):
    return "%-14s %6s %s"%(self.tag,self.size,self.dst)
#   return "%s (%s) --> %s"%(self.tag,self.size,self.dst)
#   return "%s %s --> %s"%(self,self.src,self.dst)
  def __repr__(self):
    from os.path import basename
    return "%6s %s"%(self.size,basename(self.dst))
#   return "%6s %s"%(self.size,self.tag)
#   return "%s (%s)"%(self.tag,self.size)

  def __hash__(self):
    '''find unique self.src occurences'''
    return hash(repr(self.src))
  def __eq__(self, other):
    if isinstance(other, dataPoint):
      return (self.src == other.src)
    else:
      return False
  def __ne__(self, other):
    return (not self.__eq__(other))

  def cleanup(self,verbose=True):
    from os.path import isfile,dirname
    from os import remove,rmdir
    if not isfile(self.dst):
      return

    if(verbose):
      print("%-8s %s"%('Cleanup',self))

    try:
      remove(self.dst)
      d=dirname(self.dst)
      if d!='':
        rmdir(d)
    except OSError as e:
      if e.errno==1:  # opperation not permited (permissions?)
        pass
      if e.errno==39: # directory was not empty
        pass
      else:
        print(e)
        raise
      
  def download(self,verbose=True):
    """derived from http://stackoverflow.com/a/22776/2576368"""
    from os.path import isfile,dirname,isdir
    from os import makedirs
    import urllib2
  
    if(isfile(self.dst)):
      if(verbose):
        print("%-8s %s"%('Found',self))
    else:
      if(verbose):
        print("%-8s %s"%('Download',self))

      d=dirname(self.dst)
      if (not isdir(d))and(d!=''):
        makedirs(d)
      
      u = urllib2.urlopen(self.src)
      with open(self.dst,'wb') as f:
        block = 1024*8        #   8K
        bigFile=1024*1024*128 # 128M
        bigFile=(self.size>bigFile)
        if verbose and bigFile:
          n=0
          ntot = self.size/block
          printProgress(n,ntot,barLength=50)
        while True:
          buffer = u.read(block)
          if not buffer:
            break
          f.write(buffer)
          if verbose and bigFile:
            n+=1
            if(n%128==0):     # 8K*128=1M
              printProgress(n,ntot,barLength=50)
        if verbose and bigFile and ntot%128!=0:
          printProgress(ntot,ntot,barLength=50)

class dataSet(object):
  '''Info and retrieval/check methods'''

  def __init__(self,tag,release,year,status,dataset=None,byteSize=None):
    '''Initialize object'''
    self.tag    = str(tag)            # revision tag
    self.release= int(release)        # release date (YYYYMM)
    self.year   = year                # met-year(meteo)/status-year(model)
    self.status = status              # Status report
    self.dataset= dataset             # {'input':input,'meteo':meteo,..}
    self.size   = 0
    if byteSize:
      self.size = byteSize or 0
    else:
      self.size = sum([x.size for _,x in self.dataset.items() if hasattr(x,"size")])
    self.size = fileSize(self.size)

  def __str__(self):
    return "%-8s (meteo:%s, status:%s)"%(self.tag,self.year,self.status)
  def __repr__(self):
    return "%s: %s"%(self,{k:"%s"%v for k,v in self.dataset.items() if v})
#   return "%s: %s"%(self,self.dataset)

def readCatalog(filename,verbose=1):
  '''
  Returns releases read from catalog csv-file

  Definitions
    dataPoint(class): single remote file/tarfile
    dataSet(class):   dataPoints from a model release,
                      sorted into meteo|input|output|source|docs
    catalog(list):    all dataPoints
    index(dict):      catalog sorted into meteo|input|output|source|docs
    archive(list):    all dataSet
  '''
  from csv import reader as csvreader
  
  try:
    f=open(filename)
  except IOError as e:
    print "Failed to open '%s'.\n%s."%(filename,e.strerror)
    sys.exit(-1)

  reader = csvreader(f,delimiter=',')
  next(reader)              # skip header
  catalog = []              # list all src files (1 file per row)
  rels,keys = set(),set()   # unique dataPoint.release/.keys for indexing

  if verbose>2:
    print("Reading %s"%filename)
  for row in reader:
    if row and "".join(row).strip():  # skip empty lines
      try:
        catalog+=[dataPoint(*row)]
        rels.add(catalog[-1].release) # unique releases
        keys.add(catalog[-1].key)     # unique keys (meteo,source,&c)
        if verbose>2:
          print("  %s"%catalog[-1])
      except:
        print("Failed to parse (%s): %r"%(filename,row))
        print(dataPoint(*row))
        raise
  if verbose>1:
    print("%s read(srcs:%d)"%(filename,len(catalog)))
  f.close()
    
  index=dict.fromkeys(rels) # index[releases][keys:meteo,source,&c]
  if verbose>2:
    print("Indexing")
  for r in rels:            # r:release
    index[r]=dict.fromkeys(keys)
    for k in keys:          # k:meteo,source,&c
      index[r][k]=[x for x in catalog if (x.release==r)and(x.key==k)]
      if verbose>2:
        print("  index[%s][%s](srcs:%d)"%(r,k,len(index[r][k])))
  if verbose>1:
    print("%s index[release:%d][sets:%d]"%(filename,len(rels),len(keys)))

  # repack index from dict({key:dataPoint,}) to dict(dataSet)
  if verbose>2:
    print("Compiling releases")
  for r,v in index.items(): 
    # pack index[r] into a dataSet object
    #   before index[r]:{meteo:dataPoint,source:dataPoint,..}
    #   after  index[r]:dataSet (with search metadata)
    index[r]=dataSet(v['source'][0].model,r,v['meteo'][0].year,v['source'][0].year,v)
    if verbose>2:
      print("  index[%s]: %s"%(r,index[r]))
  if verbose>1:
    print("%s index(release:%d)"%(filename,len(rels)))

  return index

def parse_arguments():
  from optparse import OptionParser,OptionGroup
  from sys      import argv as args

  usage = """usage: %prog [options]

Examples:

  Retrieve release dataset for revision REV ({REV})
    %prog -R REV          

  Get Only the source code and user guide for revision REV
    %prog -R REV -sd

  Download meteorological input for YEAR ({MET})
    %prog -Y YEAR -m
""".format(
  REV="|".join(['rv3','v201106','rv4_0','rv4_3','rv4_4','rv4_5','rv4_8']),
  MET="|".join(['2005','2008','2010..%d'%_CONST['LASTMET']]))

  parser = OptionParser(usage,version=_CONST['VERSION'])
  parser.set_defaults(verbose=1)
  parser.add_option("-q", "--quiet",
    action="store_false", dest="verbose",
    help="don't print status messages to stdout")
  parser.add_option("-v", "--verbose",
    action="count", dest="verbose",
    help="Increase verbosity")
  parser.add_option("--catalog", default=_CONST['CSV'],
    action="store", type="string", dest="catalog",
    help="Override dataset cataloque path/file (default:%default)")

  group = OptionGroup(parser, "Release options",
    "Select a release dataset")
  group.add_option("-R","--revision",
    type="string", metavar="REV",
    action="store", dest="tag",
    help="revision REV")
  group.add_option("-S","--status",
    type="int", metavar="YEAR",
    action="store", dest="status",
    help="YEAR's status report")
  group.add_option("-Y","--year",
    type="int", metavar="YEAR",
    action="store", dest="year",
    help="Meteorological/run YEAR")
  parser.add_option_group(group)

  group = OptionGroup(parser, "Dataset options",
    "Get parts of a release dataset")
  group.add_option("-m", "--meteo", const="meteo",
    action="append_const", dest="data",
    help="get meteorology input")
  group.add_option("-i", "--input", const="input",
    action="append_const",dest="data",
    help="get other input")
  group.add_option("-o", "--output",const="output",
    action="append_const", dest="data",
    help="get model benckmark")
  group.add_option("-s", "--source",const="source",
    action="append_const", dest="data",
    help="get source code for benckmark")
  group.add_option("-d", "--docs",const="docs",
    action="append_const", dest="data",
    help="get corresponding user guide")
  parser.add_option_group(group)

  group = OptionGroup(parser, "Download options","")
  group.add_option("--yes", default=True,
    action="store_false", dest="ask",
    help="Don't ask before start downloading")
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
 
  opts,args = parser.parse_args(args[1:])
  if(all(getattr(opts,attr) is None for attr in ['tag','status','year'])):
    opts.tag=_CONST['LASTREL']
  if opts.data==None :
    opts.data=["meteo","input","output","source","docs"]
  if opts.outpath:
    _CONST['DATADIR']=opts.outpath
  if opts.tmppath:
    _CONST['TMPDIR']=opts.tmppath
  
  return opts,args

if __name__ == "__main__":
  import sys
  opts,args = parse_arguments()
 
  catalog=readCatalog(opts.catalog,opts.verbose)

  down=[] # files to download
  for attr in ['tag','status','year']:
    target=getattr(opts,attr)
    if(target is None):
      continue

    if opts.verbose>1:
      print("Searching %s:%s"%(attr,target))
    try:
      dataSet=[v for _,v in catalog.items() if getattr(v,attr)==target]
      if len(dataSet)==0:
        raise
    except:
      print("No datasets found for --%s=%s"%(attr,target))
      sys.exit(-1)
    if opts.verbose>1:
      for x in dataSet:
        print("  Found %s"%x)

    if opts.verbose>1:
      print("Searching datasets:%s"%(opts.data))
    for ds in dataSet:
      try:
        ds={key:ds.dataset[key] for key in opts.data}
      except:
        print("No datasets found for --%s=%s"%(attr,target))
        sys.exit(-1)
      if opts.verbose>1:
        for key in ds:
          print("  Found %-6s:%s"%(key,ds[key]))

      for key in ds:
        down+=ds[key]
        if opts.verbose:
          print("Queue download: %6s %s"%(
            fileSize(sum([x.size for x in ds[key]])),ds[key][0].tag))

  down=list(set(down))  # unique files, for single download and total size
  if opts.verbose:
    print("Queue download: %6s %s"%(
      fileSize(sum([x.size for x in down])),'Total'))

  while opts.ask:
    if sys.version_info[0]>2: # Python3.x
      data = input("Do you wish to proceed? [Y]/n:")
    else:                     # Python2.x
      data = raw_input("Do you wish to proceed? [Y]/n:")
    if data.lower() in ('y','yes',''):
      break
    elif data.lower() in ('n','no'):
      print("OK, bye")
      sys.exit(0)
    
  for x in down:
    x.download(opts.verbose)

    if opts.cleanup:
      x.cleanup(opts.verbose)

