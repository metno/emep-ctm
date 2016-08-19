# Open Source EMEP/MSC-W model
Simplified access to the source code, input data and benchmark results

## Get the latest release dataset
```bash
# download the catalog tool/files
wget https://github.com/metno/emep-ctm/blob/master/tools/catalog.py
wget https://github.com/metno/emep-ctm/blob/master/tools/catalog.csv

# make it executable and run it
chmod +x catalog.py
catalog.py
```

## Get partial/older datasets
```
Usage: catalog.py [options]

Examples:

  Retrieve release dataset for revision REV (rv3|v201106|rv4_0|rv4_3|rv4_4|rv4_5|rv4_8)
    catalog.py -R REV          

  Get Only the source code and user guide for revision REV
    catalog.py -R REV -sd

  Download meteorological input for YEAR (2005|2008|2010..2013)
    catalog.py -Y YEAR -m


Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -q, --quiet           don't print status messages to stdout
  -v, --verbose         Increase verbosity
  --catalog=CATALOG     Override dataset cataloque path/file
                        (default:./catalog.csv)

  Release options:
    Select a release dataset

    -R REV, --revision=REV
                        revision REV
    -S YEAR, --status=YEAR
                        YEAR's status report
    -Y YEAR, --year=YEAR
                        Meteorological/run YEAR

  Dataset options:
    Get parts of a release dataset

    -m, --meteo         get meteorology input
    -i, --input         get other input
    -o, --output        get model benckmark
    -s, --source        get source code for benckmark
    -d, --docs          get corresponding user guide

  Download options:
    --yes               Don't ask before start downloading/unpacking
    --outpath=OUTPATH   Override output (dataset) directory
    --tmppath=TMPPATH   Override temporary (download) directory
    --cleanup           Remove ALL temporary (download) files
```
