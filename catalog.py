#!/usr/bin/env python3
"""
Open Source EMEP/MSC-W model
simplified access to the source code, input data and benchmark results.
"""

import hashlib
import json
import shutil
import sys
import tarfile
from contextlib import closing
from pathlib import Path
from textwrap import dedent

assert sys.version_info >= (3, 6), "This script requires python3.6 or better"

_CONST = {
    "VERSION": "0.3.0",  # script version
    "RELEASE": [  # released model versions
        "rv3",
        "v201106",
        "rv4_0",
        "rv4_3",
        "rv4_4",
        "rv4_5",
        "rv4_8",
        "rv4_10",
        "rv4_15",
        "rv4_17",
        "rv4_32",
        "rv4_33",
        "rv4_34",
        "rv4_36",
        "rv4_45",
        "5.0"
    ],
    "METYEAR": [  # released met-years
        2005,
        2008,
        2010,
        2011,
        2012,
        2013,
        2014,
        2015,
        2017,
    ],
    "THREDDS": "https://projects.met.no/emep/Thredds_Meteo",
    "FTP": "https://projects.met.no/emep",
    "GIT": "https://github.com/metno/emep-ctm/",
    "DOC": "https://emep-ctm.readthedocs.io/",
    "RTD": "https://emep-ctm.readthedocs.io/_/downloads/en",
    "CSV": Path(__file__).parent / "catalog.csv",  # list all files from all releases
    "RAW": "https://raw.githubusercontent.com/metno/emep-ctm/tools/catalog.csv",  # catalog on the repo
    "TMPDIR": "./downloads",  # temp path for downloads
    "DATADIR": Path("."),  # base path for datasets
}


def parse_arguments(args):
    """Arguments from command line"""
    from optparse import SUPPRESS_HELP, OptionGroup, OptionParser

    usage = f"""
            usage: %prog [options]

            Examples:

            Retrieve release dataset for revision REV ({"|".join(_CONST["RELEASE"])})
                %prog -R REV

            Get Only the source code and user guide for revision REV
                %prog -R REV --source --docs

            Download meteorological input for YEAR ({"|".join(str(y) for y in _CONST["METYEAR"])})
                %prog -Y YEAR --meteo
            """

    parser = OptionParser(dedent(usage), version=_CONST["VERSION"])
    parser.set_defaults(verbose=1)
    parser.add_option(
        "-q",
        "--quiet",
        action="store_false",
        dest="verbose",
        help="don't print status messages to stdout",
    )
    parser.add_option(
        "-v", "--verbose", action="count", dest="verbose", help="Increase verbosity"
    )
    parser.add_option(
        "--catalog",
        default=_CONST["CSV"],
        action="store",
        type="string",
        dest="catalog",
        help="Override dataset cataloque path/file (default:%default)",
    )

    group = OptionGroup(parser, "Release options", "Select release dataset")
    group.add_option(
        "-R",
        "--revision",
        type="string",
        metavar="REV",
        action="append",
        dest="tag",
        help="revision REV",
    )
    group.add_option(
        "-S",
        "--status",
        type="int",
        metavar="YEAR",
        action="append",
        dest="status",
        help=SUPPRESS_HELP,
    )
    group.add_option(
        "-Y",
        "--year",
        type="int",
        metavar="YEAR",
        action="append",
        dest="year",
        help="Meteorological/run YEAR",
    )
    parser.add_option_group(group)

    group = OptionGroup(parser, "Dataset options", "Partial release dataset")
    group.add_option(
        "-m",
        "--meteo",
        const="meteo",
        action="append_const",
        dest="data",
        help="get meteorology input",
    )
    group.add_option(
        "--met-domain",
        default=None,
        action="store",
        type="string",
        dest="domain",
        help="get only DOMAIN meteorology",
    )
    group.add_option(
        "-i",
        "--input",
        const="input",
        action="append_const",
        dest="data",
        help="get other input",
    )
    group.add_option(
        "-o",
        "--output",
        const="output",
        action="append_const",
        dest="data",
        help="get model benchmark",
    )
    group.add_option(
        "-s",
        "--source",
        const="source",
        action="append_const",
        dest="data",
        help="get source code for benchmark",
    )
    group.add_option(
        "-d",
        "--docs",
        const="docs",
        action="append_const",
        dest="data",
        help="get corresponding user guide",
    )
    group.add_option(
        "--extras",
        action="store_true",
        dest="extras",
        help="also get the extras, if any",
    )
    parser.add_option_group(group)

    group = OptionGroup(parser, "Download options", "")
    group.add_option(
        "--yes",
        default=True,  # help="YEAR's status report")
        action="store_false",
        dest="ask",
        help="Don't ask before start downloading/unpacking",
    )
    group.add_option(
        "--outpath",
        action="store",
        dest="outpath",
        help="Override output (dataset) directory",
    )
    group.add_option(
        "--tmppath",
        action="store",
        dest="tmppath",
        help="Override temporary (download) directory",
    )
    group.add_option(
        "--cleanup",
        default=False,
        action="store_true",
        dest="cleanup",
        help="Remove ALL temporary (download) files",
    )
    parser.add_option_group(group)

    opts, args = parser.parse_args(args)
    if all(getattr(opts, attr) is None for attr in ["tag", "status", "year"]):
        opts.tag = [_CONST["RELEASE"][-1]]
    if opts.data is None:
        opts.data = ["meteo", "input", "output", "source", "docs"]
    if opts.extras:
        opts.data += ["extra"]
    if opts.outpath:
        _CONST["DATADIR"] = Path(opts.outpath)
    if opts.tmppath:
        _CONST["TMPDIR"] = opts.tmppath
    return opts, args


def user_consent(question, ask=True, default="yes"):
    """Ask user for confirmation"""

    # format question
    qq = {"y": "%s [Y]/n:", "yes": "%s [Y]/n:", "n": "%s y/[N]:", "no": "%s y/[N]:"}
    try:
        question = qq[default.lower()] % question
    except KeyError:
        print(f"Unsupported option: user_consent(..,default='{default}')")

    # valid answers
    answer = {"y": True, "yes": True, "n": False, "no": False}

    # ask until you get a valid answer
    while True:
        if not ask:  # take the default
            print(question + default)
            data = default
        else:
            data = input(question)

        # use default when user replies ''
        data = data.lower() if len(data) else default.lower()
        if data in answer:
            return answer[data]


# Print iterations progress
# http://stackoverflow.com/a/34325723/2576368
def print_progress(iteration, total, prefix="", suffix="", decimals=2, length=100):
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
    progress = "█" * filled + "-" * (length - filled)
    sys.stdout.write(f"\r{prefix} |{progress}| {percents}% {suffix}")
    sys.stdout.flush()
    if iteration == total:
        sys.stdout.write("\n")
        sys.stdout.flush()


def file_size(value):
    """
    Human readable file size (K/M/G/T bytes)
      derived from http://stackoverflow.com/a/1094933/2576368
    """
    num = max(float(value), 0)
    for unit in ["B", "K", "M", "G"]:
        if num < 1024.0:
            return f"{num:.1f}{unit}"
        num /= 1024.0
    return f"{num:.1f}T"


class DataPoint:
    """Info and retrieval/check methods"""

    def __init__(self, release, key, year, model, src, dst="", byteSize=0, md5sum=None):
        """Initialize object"""

        self.release = int(release)  # release date (YYYYMM)
        self.key = str(key)
        self.tag = {  # eg 'rv4_8source'
            "other": "{KEY}{REL}",
            "meteo": "{KEY}{YEAR}",
        }.get(self.key, "{MOD}{KEY}")
        try:
            self.year = int(year)  # met-year
        except ValueError:
            self.year = None
        try:
            self.model = str(model)  # model version (or meteo domain)
        except ValueError:
            self.model = None
        self.src: str  # single source url/file
        self.dst: Path  # path for uncompressed self.dst
        self.size = float(byteSize)  # self.src file size [bytes]
        self.md5sum = md5sum  # self.src checksum

        # replace keywords
        kwargs = _CONST
        kwargs.update(REL=self.release, KEY=self.key, YEAR=self.year, MOD=self.model)
        self.tag = self.tag.format_map(kwargs)
        kwargs.update(TAG=self.tag)
        self.src = src.format_map(kwargs)
        if isinstance(dst, str):
            dst = dst.format_map(kwargs)
        self.dst = Path(dst)
        if self.dst.name == "":
            self.dst /= Path(self.src).name

    def __str__(self) -> str:
        return f"{self.tag:>14} {self.file_size:>6} {self.dst}"

    def __repr__(self) -> str:
        return f"{self.file_size:>6} {self.dst.name}"

    def __hash__(self) -> int:
        """find unique self.src occurrences"""
        return hash(repr(self.src))

    def __eq__(self, other) -> bool:
        if isinstance(other, DataPoint):
            return self.src == other.src
        else:
            return False

    def __ne__(self, other) -> bool:
        return not self.__eq__(other)

    @property
    def file_size(self) -> str:
        return file_size(self.size)

    def cleanup(self, verbose=True):
        """Remove (raw) downloads"""
        if not self.dst.is_file():
            return
        if verbose:
            print(f"Cleanup  {self}")

        try:
            self.dst.unlink()
            if self.dst.parent != "":
                self.dst.parent.rmdir()
        except OSError as error:
            if error.errno == 1:  # opperation not permitted (permissions?)
                pass
            elif error.errno == 39:  # directory was not empty
                pass
            else:
                raise error

    def check(self, verbose=True, cleanup=False):
        """Check download against md5sum"""
        if not self.dst.is_file():
            return False
        if verbose:
            print(f"Check    {self}")

        if self.md5sum != hashlib.md5(self.dst.read_bytes()).hexdigest():
            if verbose:
                print(f"  md5 /= {self.md5sum}")
            if cleanup:  # remove broken file
                self.cleanup(verbose)
            return False
        return True

    def download(self, verbose=True):
        """derived from http://stackoverflow.com/a/22776/2576368"""
        from urllib.request import HTTPError, urlopen

        # check if file exists/md5sum, remove file if fail md5sum
        if self.check(verbose > 2, cleanup=True):
            return

        if verbose:
            print(f"Download {self}")

        try:
            url = urlopen(self.src)
        except HTTPError:
            print(f"Could not download {self.src}")
            sys.exit(-1)

        self.dst.parent.mkdir(parents=True, exist_ok=True)
        with url, self.dst.open("wb") as file:
            block = 1_024 * 8  #   8K
            big_file = self.size > 1_048_576 * 128  # 128M
            total = self.size if verbose and big_file else 0
            n = 0
            while True:
                buff = url.read(block)
                if not buff:
                    break
                file.write(buff)
                n += len(buff)
                if total and (n == total or n % 1_048_576 == 0):  # 1M
                    print_progress(n, total, length=50)

    def unpack(self, verbose=True, inspect=False):
        """Unpack download"""
        if not self.check(verbose > 2):
            return

        if tarfile.is_tarfile(self.dst):
            print(f"Untar    {self}")
            with closing(tarfile.open(self.dst, "r")) as file:
                if user_consent("  See the contents first?", inspect, "no"):
                    print(file.list(verbose=(verbose > 1)))
                    if not user_consent("  Do you wish to proceed?", inspect):
                        print("OK, skipping file")
                        return

                _CONST["DATADIR"].parent.mkdir(parents=True, exist_ok=True)
                try:
                    file.extractall(_CONST["DATADIR"])
                except EOFError as error:
                    print(f"  Failed unpack '{self.dst}':\n    {error}.")
                    if not user_consent("    Do you wish to continue?", inspect, "yes"):
                        sys.exit(-1)

        elif self.dst.name == "pdf":
            outfile = _CONST["DATADIR"] / self.key / self.model / "emep-ctm.pdf"
            outfile.parent.mkdir(parents=True, exist_ok=True)
            self.dst.rename(outfile)
        else:
            outfile = _CONST["DATADIR"] / self.key / self.dst.name
            if verbose > 1:
                print(f"Copy     {outfile}")
            outfile.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(self.dst, outfile)


class DataSet:
    """Info and retrieval/check methods"""

    def __init__(self, tag, release, year, status, dataset=None, byteSize=None):
        """Initialize object"""
        self.tag = str(tag)  # revision tag
        self.release = int(release)  # release date (YYYYMM)
        self.year = year  # met-year(meteo)/status-year(model)
        self.status = status  # Status report
        self.dataset = dataset  # {'input':input,'meteo':meteo,..}
        self.size = 0
        if byteSize:
            self.size = float(byteSize) or 0
        else:
            self.size = sum(x.size for x in self.dataset.values() if hasattr(x, "size"))

    def __str__(self):
        return f"{self.tag:>8} (meteo:{self.year}, status:{self.status})"

    def __repr__(self):
        return f"{self}: {json.dumps(self.dataset)}"


def read_catalog(filename: Path, verbose: int = 1):
    """
    Returns releases read from catalog csv-file

    Definitions
      DataPoint(class): single remote file/tarfile
      DataSet(class):   dataPoints from a model release,
                        sorted into meteo|input|output|source|docs
      catalog(list):    all dataPoints
      index(dict):      catalog sorted into meteo|input|output|source|docs
      archive(list):    all DataSet
    """
    from csv import reader

    # download catalog if file not found
    if not filename.is_file():
        DataPoint(0, "catalog", 0, "", _CONST["RAW"], filename).download(verbose)

    # catalog file should be present at this point
    with open(filename) as file:
        csv = reader(file, delimiter=",")
        next(csv)  # skip header
        catalog: list = []  # list all src files (1 file per row)
        rels, keys = set(), set()  # unique DataPoint.release/.keys for indexing

        if verbose > 2:
            print(f"Reading {filename}")
        for row in csv:
            if not row or "".join(row).strip() == "":
                continue  # skip empty lines
            try:
                catalog += [DataPoint(*row)]
                rels.add(catalog[-1].release)  # unique releases
                keys.add(catalog[-1].key)  # unique keys (meteo,source,&c)
                if verbose > 2:
                    print(f"  {catalog[-1]}")
            except:
                print("Failed to parse ({filename}): {row}")
                print(DataPoint(*row))
                raise
        if verbose > 1:
            print(f"{filename} read(srcs:{len(catalog)})")

    index = dict.fromkeys(rels)  # index[releases][keys:meteo,source,&c]
    if verbose > 2:
        print("Indexing")
    for r in rels:  # r:release
        index[r] = dict.fromkeys(keys)
        for k in keys:  # k:meteo,source,&c
            index[r][k] = [x for x in catalog if (x.release == r) and (x.key == k)]
            if verbose > 2:
                print(f"  index[{r}][{k}](srcs:{len(index[r][k])})")
    if verbose > 1:
        print(f"{filename} index[release:{len(rels)}][sets:{len(keys)}]")

    # repack index from dict({key:DataPoint,}) to dict(DataSet)
    if verbose > 2:
        print("Compiling releases")
    for r, v in index.items():
        # pack index[r] into a DataSet object
        #   before index[r]:{meteo:DataPoint,source:DataPoint,..}
        #   after  index[r]:DataSet (with search metadata)
        index[r] = DataSet(
            v["source"][0].model,
            r,
            v["meteo"][0].year if v["meteo"] else None,
            v["source"][0].year,
            v,
        )
        if verbose > 2:
            print(f"  index[{r}]: {index[r]}")
    if verbose > 1:
        print(f"{filename} index(release:{len(rels)})")

    return index


def main(opts):
    """Command line function"""

    def get_datasets(catalog, attr, target):
        """Search catalog for tag|status|year DataSet"""
        if opts.verbose > 1:
            print(f"Searching {attr}(s):{target}")
        try:
            dataset = [v for v in catalog.values() if getattr(v, attr) in target]
            if len(dataset) == 0:
                raise IndexError
        except IndexError:
            print(f"No datasets found for --{attr}={target}")
            sys.exit(-1)
        if opts.verbose > 1:
            for x in dataset:
                print(f"  Found {x}")
        return dataset

    def get_downloads(catalog, attr, target):
        """List downloads for tag|status|year DataSet"""
        downloads = []
        if opts.verbose > 1:
            print(f"Searching datasets:{opts.data}")
        for ds in get_datasets(catalog, attr, target):
            try:
                ds = {key: ds.dataset[key] for key in opts.data}
                # only download meteo with matching --met-domain option
                if "meteo" in ds and opts.domain:
                    ds["meteo"] = [x for x in ds["meteo"] if x.model == opts.domain]
            except KeyError:
                print(f"No datasets found for --{attr}={target}")
                sys.exit(-1)
            if opts.verbose > 1:
                for key in ds:
                    print(f"  Found {key:>6}:{ds[key]}")

            for key in ds:
                # do not try to download 0-size files
                aux = [x for x in ds[key] if x.size > 0]
                total = sum(x.size for x in aux)
                if total > 0:
                    downloads += aux
                if opts.verbose and aux:
                    print(f"Queue download: {file_size(total):>6} {aux[0].tag}")
        return downloads

    # list files to download
    catalog = read_catalog(opts.catalog, opts.verbose)
    downloads = []  # files to download
    for attr in ["tag", "status", "year"]:
        target = getattr(opts, attr)
        if target is None:
            continue
        downloads += get_downloads(catalog, attr, target)

    # list file_sizes to download
    downloads = list(set(downloads))  # unique files, for single download and total size
    total = file_size(sum(x.size for x in downloads))
    print(f"Queue download: {total:>6} Total")
    if not user_consent("Do you wish to proceed?", opts.ask):
        print("OK, bye")
        sys.exit(0)

    # download files
    for x in downloads:
        x.download(opts.verbose)
        x.unpack(opts.verbose, opts.ask)
        if opts.cleanup:
            x.cleanup(opts.verbose)


if __name__ == "__main__":
    main(parse_arguments(sys.argv[1:] or ["-h"])[0])
