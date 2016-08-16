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

  def __init__(self,tag,src,dst='{DATADIR}/{TAG}',size=0,cksum=None):
    '''Initialize object'''
    self.tag    = str(tag)            # revision tag
    self.src    = str(src)            # list source url/file(s)
    self.dst    = str(dst)            # path for uncompressed files
    self.size   = int(size)           # disk space required in bytes
    self.cksum  = str(cksum)          # checksum url/file for uncompressed files

    # replace _CONST keywords
    self.src    = self.src.format(TAG=self.tag,**_CONST)
    self.dst    = self.dst.format(TAG=self.tag,**_CONST)
    self.cksum  = self.cksum.format(TAG=self.tag,**_CONST)

  def __repr__(self):
#   return "%s(%d bytes): %s"%(self.tag,self.size,self.dst)
    return "%s (%d bytes)"%(self.dst,self.size)

class release(object):
  '''Model release info'''

  def __init__(self,tag,year,status,dataset=None):
    '''Initialize object'''
    self.tag    = tag                 # revision tag
#   self.date   = date                # release date
    self.year   = year                # met-year
    self.status = status              # Status report
    self.dataset= dataset

  def __repr__(self):
    return "%s (meteo:%s, status:%s)"%(self.tag,self.year,self.status)
#   return "%s (meteo:%s, status:%s): %s"%(self.tag,self.year,self.status,self.dataset)

archive={
  '201510':{
    'source':dataset('rv4_8',
      ['{FTP}/201510/model_code/EMEP_MSC-W_model.rv4.8.OpenSource.tar.gz'],
      '{DATADIR}/source/Unimod.{TAG}'),
    'docs':dataset('rv4_8',
      ['{FTP}/201510/userguide_{TAG}.pdf'],
      '{DATADIR}/docs/'),
    'meteo':dataset('met2013',
      ['{FTP}/201510/input_data/meteo/meteo%d.tar.bz2'%i for i in range(201301,201312)],
      '{DATADIR}/input/{TAG}'),
    'input':dataset('rv4_8',
      ['{FTP}/201510/input_data/other_input_files.tar.bz2'],
      '{DATADIR}/input/{TAG}'),
    'output':dataset('rv4_8_2013',
      ['{FTP}/201510/model_results_rv4_8_2013/Base_model_rv4_8_results_NoFF.tar.bz2',
       '{FTP}/201510/model_results_rv4_8_2013/Base_model_rv4_8_results_withFF.tar.bz2',
       '{FTP}/201510/model_results_rv4_8_2013/README'],
      '{DATADIR}/output/{TAG}'),
    },
  '201409':{
    'source':dataset('rv4_5',
      ['{FTP}/201409/model_code/EMEP_MSC-W_model.rv4_5.OpenSource.tar.gz'],
      '{DATADIR}/source/Unimod.{TAG}'),
    'docs':dataset('rv4_5',
      ['{FTP}/201409/userguide_{TAG}.pdf'],
      '{DATADIR}/docs/'),
    'meteo':dataset('met2012',
      ['{FTP}/201409/input_data_2012/meteo/meteo%d.tar.bz2'%i for i in range(201201,201212)],
      '{DATADIR}/input/{TAG}'),
    'input':dataset('rv4_5',
      ['{FTP}/201409/input_data_2012/other_input_files.tar.bz2'],
      '{DATADIR}/input/{TAG}'),
    'output':dataset('rv4_5_2012',
      ['{FTP}/201409/model_results_rv4_5_2012/Base_model_rv4_5_results.tar.bz2'],
      '{DATADIR}/output/{TAG}'),
    },
  '201309':{
    'source':dataset('rv4_4',
      ['{FTP}/201309/model_code/EMEP_MSC-W_model.rv4_4.OpenSource.tar.gz'],
      '{DATADIR}/source/Unimod.{TAG}'),
    'docs':dataset('rv4_4',
      ['{FTP}/201309/userguide_{TAG}.pdf'],
      '{DATADIR}/docs/'),
    'meteo':dataset('met2011',
      ['{FTP}/201309/input_data_2011/meteo/meteo%d.tar.bz2'%i for i in range(201101,201112)],
      '{DATADIR}/input/{TAG}'),
    'input':dataset('rv4_4',
      ['{FTP}/201309/input_data_2011/other_input_files.tar.bz2'],
      '{DATADIR}/input/{TAG}'),
    'output':dataset('rv4_4_2011',
      ['{FTP}/201309/model_results_rv4_4_2011/Base_model_rv4_4_results.tar.bz2'],
      '{DATADIR}/output/{TAG}'),
#   'output':dataset('rv4_4_2010',
#     ['{FTP}/201309/model_results_rv4_4_2010/Base_model_rv4_4_results.tar.bz2'],
#     '{DATADIR}/output/{TAG}'),
    },
  '201304':{
    'source':dataset('rv4_3',
      ['{FTP}/201304/model_code/EMEP_MSC-W_model.rv4_3.OpenSource.tar.gz',
       '{FTP}/201304/README_rv4_3special.txt'],
      '{DATADIR}/source/Unimod.{TAG}'),
    'docs':dataset('rv4_3',
      ['{FTP}/201304/userguide_{TAG}.pdf'],
      '{DATADIR}/docs/'),
    'meteo':dataset('met2010',
      ['{FTP}/met2010/meteo%d.tar.bz2'%i for i in range(201001,201012)],
      '{DATADIR}/input/{TAG}'),
    'input':dataset('rv4_3',
      ['{FTP}/201304/input_data/other_input_files.tar.bz2'],
      '{DATADIR}/input/{TAG}'),
    'output':dataset('rv4_3_2010',
      ['{FTP}/201304/model_results_rv4_3_2010/Base_model_rv4_3_results.tar.bz2'],
      '{DATADIR}/output/{TAG}'),
#   'tools':dataset('NCL_LOCAL',
#     ['201304/mscw-osrc_ncl.tar.gz'],
#     '{DATADIR}/tools/{TAG}'),
    },
  '201209':{
    'source':dataset('rv4_0',
      ['{FTP}/201209/model_code/Unified_EMEP_model.tar.bz2'],
      '{DATADIR}/source/Unimod.{TAG}'),
    'docs':dataset('rv4_0',
      ['{FTP}/201209/userguide092012.pdf'],
      '{DATADIR}/docs/userguide_{TAG}.pdf'),
    'meteo':dataset('met2010',
      ['{FTP}/met2010/meteo%d.tar.bz2'%i for i in range(201001,201012)],
      '{DATADIR}/input/{TAG}'),
    'input':dataset('rv4_0',
      ['{FTP}/201209/input_data/other_input_files.tar.bz2'],
      '{DATADIR}/input/{TAG}'),
    'output':dataset('rv4_0_2010',
      ['{FTP}/201209/model_results_2010/Base_model_results.tar.bz2'],
#      '{FTP}/201209/model_results_2010/April_model_results.tar.bz2'],
      '{DATADIR}/output/{TAG}'),
    },
  '201108':{
    'source':dataset('v201106',
      ['{FTP}/201108/model_code/MSC-W_EMEP_ModelCode.tar.bz2'],
      '{DATADIR}/source/Unimod.{TAG}'),
    'docs':dataset('v201106',
      ['{FTP}/201108/userguide_062011.pdf'],
      '{DATADIR}/docs/userguide_{TAG}.pdf'),
    'meteo':dataset('met2008',
      ['{FTP}/201108/input_data/meteo%d.tar.bz2'%i for i in range(200801,200812)],
      '{DATADIR}/input/{TAG}'),
    'input':dataset('v201106',
      ['{FTP}/201108/input_data/other_input_files.tar.bz2'],
      '{DATADIR}/input/{TAG}'),
    'output':dataset('v201106_2008',
      ['{FTP}/201108/model_results_2008/model_results_Base.tar.bz2'],
#      '{FTP}/201108/model_results_2008/model_results_Jan_01_02.tar.bz2'],
      '{DATADIR}/output/{TAG}'),
#   'tools':dataset('tools',
#     ['201108/model_tools/tools.tar.bz2'],
#     '{DATADIR}/tools/{TAG}'),
    },
  '200802':{
    'source':dataset('rv3',
      ['{FTP}/200802/model_code/Unified_EMEP.tar'],
      '{DATADIR}/source/Unimod.{TAG}'),
    'docs':dataset('rv3',
      ['{FTP}/200802/User_Guide_Unified_EMEP_rv3.pdf'],
      '{DATADIR}/docs/userguide_{TAG}.pdf'),
    'meteo':dataset('met2005',
      ['{FTP}/201108/input_data/meteo%d.tar.bz2'%i for i in range(200501,200512)]
     +['{FTP}/200802/input_data/meteo/meteo20060101.nc.bz2'],
      '{DATADIR}/input/{TAG}'),
    'input':dataset('rv3',
      ['{FTP}/200802/input_data/other/other_input_files.tar.bz2'],
      '{DATADIR}/input/{TAG}'),
    'output':dataset('rv3_2005',
      ['200802/model_results_2005/Results_2005.tar.bz2'],
      '{DATADIR}/output/{TAG}'),
    },
}

OpenSource=[
  release('rv4.8'    ,2013,2015,archive['201510']),
  release('rv4.5'    ,2012,2014,archive['201409']),
  release('rv4.4'    ,2011,2013,archive['201309']),
  release('rv4.3'    ,None,None,archive['201304']),
  release('rv4.0'    ,2010,2012,archive['201209']),
  release('v.2011-06',2008,None,archive['201108']),
  release('rv3'      ,2005,None,archive['200802']),
]

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
  REV="|".join([    x.tag     for x in OpenSource]),
  MET="|".join([str(x.year)   for x in OpenSource if x.year]),
  REP="|".join([str(x.status) for x in OpenSource if x.status]))

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
 
  opts,args = parser.parse_args(args[1:])
  if(all(getattr(opts,attr) is None for attr in ['tag','status','year'])):
    opts.tag=_CONST['LASTREL']
  if(opts.data==None):
    opts.data=["meteo","input","output","source","docs"]
  
  return opts,args

if __name__ == "__main__":
  import sys
  opts,args = parse_arguments()
 
  for attr in ['tag','status','year']:
    target=getattr(opts,attr)
    if(target is None):
      continue

    if opts.verbose>1:
      print("Searching %s:%s"%(attr,target))
    try:
      dataSet=[getattr(x,attr) for x in OpenSource].index(target)
      dataSet=OpenSource[dataSet]
    except:
      print("No datasets found for --%s=%s"%(attr,target))
      sys.exit(-1)
    if opts.verbose>1:
      print("  Found %s"%dataSet)

    if opts.verbose>1:
      print("Searching datasets:%s"%(opts.data))
    try:
      dataSet=dataSet.dataset
      dataSet={key:dataSet[key] for key in opts.data}
    except:
      print("No datasets found for --%s=%s"%(attr,target))
      sys.exit(-1)
    if opts.verbose>1:
      for key in dataSet:
        print("  Found %-6s:%s"%(key,dataSet[key].tag))

    if opts.verbose:
      print("From %s:%s will retrieve:"%(attr,target))
      for key in dataSet:
        print("  %-6s:%s"%(key,dataSet[key]))
