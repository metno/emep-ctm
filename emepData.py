#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Open Source EMEP/MSC-W model
simplified access to the source code, input data and benchmark results.
"""

_CONST={'VERSION':"0.0.1",'LASTREL':"rv4.8",
  'FTP':"ftp://ftp.met.no/projects/emep/OpenSource",
  'TMPDIR':"./tmp",'DATADIR':'.'}

class dataset(object):
  '''Dataset info and retrieval/check methods'''

  def __init__(self,name,src,dst='{DATADIR}/{NAME}',size=0,cksum=None):
    '''Initialize object'''
    self.name   = str(name)           # name/revision/tag
    self.src    = str(src)            # list source url/file(s)
    self.dst    = str(dst)            # path for uncompressed files
    self.size   = int(size)           # disk space required in bytes
    self.cksum  = str(cksum)          # checksum url/file for uncompressed files

    # replace _CONST keywords
    self.src    = self.src.format(NAME=self.name,**_CONST)
    self.dst    = self.dst.format(NAME=self.name,**_CONST)
    self.cksum  = self.cksum.format(NAME=self.name,**_CONST)

  def __repr__(self):
#   return "%s(%d bytes): %s"%(self.name,self.size,self.dst)
    return "%s (%d bytes)"%(self.dst,self.size)

archive={
  '201510':{
    'source':dataset('rv4_8',
      ['{FTP}/201510/model_code/EMEP_MSC-W_model.rv4.8.OpenSource.tar.gz'],
      '{DATADIR}/source/Unimod.{NAME}'),
    'docs':dataset('userguide_rv4_8',
      ['{FTP}/201510/userguide_rv4_8.pdf'],
      '{DATADIR}/docs'),
    'meteo':dataset('meteo2013',
      ['{FTP}/201510/input_data/meteo/meteo%06d.tar.bz2'%i for i in range(201301,201312)],
      '{DATADIR}/input/meteo/2013'),
    'input':dataset('input2013',
      ['{FTP}/201510/input_data/other_input_files.tar.bz2'],
      '{DATADIR}/input/other/2013'),
    'output':dataset('output2013',
      ['{FTP}/201510/model_results_rv4_8_2013/Base_model_rv4_8_results_NoFF.tar.bz2',
       '{FTP}/201510/model_results_rv4_8_2013/Base_model_rv4_8_results_withFF.tar.bz2',
       '{FTP}/201510/model_results_rv4_8_2013/README'],
      '{DATADIR}/source/Unimod.{NAME}'),
    },
  '201409':'ftp://ftp.met.no/projects/emep/OpenSource/201409/',
  '201309':'ftp://ftp.met.no/projects/emep/OpenSource/201309/',
  '201304':'ftp://ftp.met.no/projects/emep/OpenSource/201304/',
  '201209':'ftp://ftp.met.no/projects/emep/OpenSource/201209/',
  '201108':'ftp://ftp.met.no/projects/emep/OpenSource/201108/',
  '200802':'ftp://ftp.met.no/projects/emep/OpenSource/200802/',
}


class release(object):
  '''Model release info'''

  def __init__(self,tag,year,status,dataset=None):
    '''Initialize object'''
    self.tag    = tag                 # name/revision/tag
#   self.date   = date                # release date
    self.year   = year                # met-year
    self.status = status              # Status report
    self.dataset= dataset

  def __repr__(self):
    return "%s (meteo:%s, status:%s)"%(self.tag,self.year,self.status)
#   return "%s (meteo:%s, status:%s): %s"%(self.tag,self.year,self.status,self.dataset)

OpenSource=[
  release('rv4.8'    ,2013,2015,archive['201510']),
  release('rv4.5'    ,2012,2014,archive['201409']),
  release('rv4.4'    ,2011,2013,archive['201309']),
  release('rv4.3'    ,None,None,archive['201304']),
  release('rv4.0'    ,2010,2012,archive['201209']),
  release('v.2011-06',2009,2011,archive['201108']),
  release('rv3'      ,2006,2008,archive['200802']),
]

def parse_arguments():
  usage = """usage: %prog [options]"""
  from optparse import OptionParser,OptionGroup
  from sys      import argv as args

  parser = OptionParser(usage,version=_CONST['VERSION'])
  parser.set_defaults(verbose=1)
  parser.add_option("-q", "--quiet",
    action="store_false", dest="verbose",
    help="don't print status messages to stdout")
  parser.add_option("-v", "--verbose",
    action="count", dest="verbose",
    help="Increase verbosity")

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
 
  opts,args = parser.parse_args(args[1:])
  if(all(getattr(opts,attr) is None for attr in ['tag','status','year'])):
    opts.tag=_CONST['LASTREL']
  if(opts.data==None):
    opts.data=["meteo","input","output","source","docs"]
  
  return opts,args

if __name__ == "__main__":
  opts,args = parse_arguments()
 
  for attr in ['tag','status','year']:
    target=getattr(opts,attr)
    if(target is None):
      continue

    if opts.verbose>1:
      print("Searching %s:%s"%(attr,target))
    try:
      dataSet=[getattr(rel,attr)==target for rel in OpenSource].index(True)
      dataSet=OpenSource[dataSet]
    except:
      print("Release not found")
      raise
    if opts.verbose>1:
      print("  Found %s"%dataSet)

    if opts.verbose>1:
      print("Searching datasets:%s"%(opts.data))
    try:
      dataSet=dataSet.dataset
      dataSet={key:dataSet[key] for key in opts.data}
    except:
      print("Datasets not found")
      raise
    if opts.verbose>1:
      for key in dataSet:
        print("  Found %-6s:%s"%(key,dataSet[key].name))

    if opts.verbose:
      print("From %s:%s will retrieve:"%(attr,target))
      for key in dataSet:
        print("  %-6s:%s"%(key,dataSet[key]))
