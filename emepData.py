#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Open Source EMEP/MSC-W model
simplified access to the source code, input data and benchmark results.
"""
VERSION="0.0.1"
LASTREL="rv4.8"

class unimod(object):
  '''Model release info'''

  def __init__(self,tag,year,report,dataset=None):
    '''Initialize object using the data provided by one csv row.'''
    self.tag    = tag                 # name/revision/tag
    self.year   = year                # met-year
    self.report = report              # Status report
    self.dataset= dataset

  def __repr__(self):
    return "%s(%s): Status %s"%(self.name,self.year,self.report)

OpenSource=[
  unimod("rv4.8",2013,2015,"ftp://ftp.met.no/projects/emep/OpenSource/201510/"),
  unimod("rv4.5",2012,2014,"ftp://ftp.met.no/projects/emep/OpenSource/201409/"),
  unimod("rv4.4",2011,2013,"ftp://ftp.met.no/projects/emep/OpenSource/201309/"),
  unimod("rv4.3",None,None,"ftp://ftp.met.no/projects/emep/OpenSource/201304/"),
  unimod("rv4.0",2010,2012,"ftp://ftp.met.no/projects/emep/OpenSource/201209/"),
  unimod("v.2011-06",2009,2011,"ftp://ftp.met.no/projects/emep/OpenSource/201108/"),
  unimod("rv3"  ,2006,2008,"ftp://ftp.met.no/projects/emep/OpenSource/200802/"),
]

def parse_arguments():
  usage = """usage: %prog [options]"""
  from optparse import OptionParser,OptionGroup
  from sys      import argv as args

  parser = OptionParser(usage,version=VERSION)
  parser.add_option("-q", "--quiet",
    action="store_false", dest="verbose", default=True,
    help="don't print status messages to stdout")

  group = OptionGroup(parser, "Release options",
    "Select a release dataset")
  group.add_option("-R","--release",
    type="string", metavar="rvX.Y",
    action="store", dest="tag",
    help="revision rvX.Y")
  group.add_option("-S","--status",
    type="int", metavar="YEAR",
    action="store", dest="status",
    help="YEAR's status report")
  group.add_option("-Y","--year",
    type="int", metavar="YEAR",
    action="store", dest="year",
    help="YEAR's benckmark")
  parser.add_option_group(group)


  group = OptionGroup(parser, "Data-set options",
    "Get parts of a release dataset")
  group.add_option("-m", "--meteo", const="meteo",
    action="append_const", dest="data",
    help="get meteorology input")
  group.add_option("-i", "--input", const="other",
    action="append_const",dest="data",
    help="get other input")
  group.add_option("-o", "--output",const="output",
    action="append_const", dest="data",
    help="get model benckmark")
  group.add_option("-s", "--source",const="source",
    action="append_const", dest="data",
    help="get source code for benckmark")
  parser.add_option_group(group)
 
  opts,args = parser.parse_args(args[1:])
  if(all(getattr(opts,attr)==None for attr in ['tag','status','year'])):
    opts.tag=LASTREL
  if(opts.data==None):
    opts.data=["meteo","other","output","source"]
  
  return opts,args

if __name__ == "__main__":
  from operator import attrgetter
  opts,args = parse_arguments()
  for attr in ['tag','status','year']:
    g=attrgetter(attr)
    if(g(opts)==None):
      continue
    print "Search %s = %s" % (attr, g(opts))
    print [g(rel)==g(opts) for rel in OpenSource]

