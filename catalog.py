#!/usr/bin/env python3
"""
Open Source EMEP/MSC-W model
simplified access to the source code, input data and benchmark results.
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import shutil
import sys
import tarfile
import warnings
from collections import defaultdict
from collections.abc import Collection, Iterator
from contextlib import closing
from pathlib import Path
from textwrap import dedent
from types import SimpleNamespace
from typing import NamedTuple

assert sys.version_info >= (3, 10), "This script requires python3.10 or better"
warnings.filterwarnings("ignore", r".*CVE-2007-4559", RuntimeWarning, "tarfile")

DEFAULT = SimpleNamespace(
    # script version
    version="0.5.0",
    # released model versions
    releases=(
        "rv3,v201106,rv4_0,rv4_3,rv4_4,rv4_5,rv4_8,rv4_10,"
        "rv4_15,rv4_17,rv4_32,rv4_33,rv4_34,rv4_36,rv4_45,"
        "5.0,5.5,5.6"
    ),
    # released met-years
    years="2005,2008,2010,2011,2012,2013,2014,2015,2017,2018",
    # list all files from all releases
    csv=Path(__file__).with_name("catalog.csv"),
    # catalog on the repo
    remote="https://raw.githubusercontent.com/metno/emep-ctm/tools/catalog.csv",
    # temp path for downloads
    downloads=Path("./downloads/"),
    # base path for datasets
    output=Path("."),
)

_CONST = {
    "THREDDS": "https://projects.met.no/emep/Thredds_Meteo",
    "FTP": "https://projects.met.no/emep",
    "GIT": "https://github.com/metno/emep-ctm/",
    "DOC": "https://emep-ctm.readthedocs.io/",
    "RTD": "https://emep-ctm.readthedocs.io/_/downloads/en",
}


def parse_arguments(args: list[str]) -> argparse.Namespace:
    """Arguments from command line"""

    usage = dedent(
        f"""
        %(prog)s [options]

        Examples:

        Retrieve release dataset for revision REV ({DEFAULT.releases})
            %(prog)s -R REV

        Get Only the source code and user guide for revision REV
            %(prog)s -R REV --source --docs

        Download meteorological input for YEAR ({DEFAULT.years})
            %(prog)s -Y YEAR --meteo
        """
    ).strip()
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("-V", "--version", action="version", version=f"%(prog)s {DEFAULT.version}")

    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument(
        "-q",
        "--quiet",
        action="store_const",
        dest="log_level",
        const=logging.ERROR,
        help="hide info messages",
    )
    verbosity.add_argument(
        "-v",
        "--verbose",
        action="store_const",
        dest="log_level",
        const=logging.DEBUG,
        help="show debug messages",
    )
    parser.set_defaults(log_level=logging.INFO)

    parser.add_argument(
        "--catalog",
        default=DEFAULT.csv,
        action="store",
        type=Path,
        dest="catalog",
        help="Override dataset cataloque path/file (default:%(default)s)",
    )

    group = parser.add_argument_group("Release options", "Select release dataset")
    group.add_argument(
        "-R",
        "--revision",
        type=str,
        choices=DEFAULT.releases.split(","),
        metavar="REV",
        action="store",
        dest="tag",
        help="revision REV",
    )
    group.add_argument(
        "-Y",
        "--year",
        type=int,
        choices=[int(year) for year in DEFAULT.years.split(",")],
        metavar="YEAR",
        action="store",
        dest="year",
        help="Meteorological/run YEAR",
    )

    group = parser.add_argument_group("Dataset options", "Partial release dataset")
    group.add_argument(
        "-m",
        "--meteo",
        const="meteo",
        action="append_const",
        dest="data",
        help="get meteorology input",
    )
    group.add_argument(
        "--met-domain",
        type=str,
        choices=["EMEP", "EECCA", "EMEP01", "MACC14", "EMEP0302"],
        default=None,
        action="store",
        dest="domain",
        help="get only DOMAIN meteorology",
    )
    group.add_argument(
        "-i",
        "--input",
        const="input",
        action="append_const",
        dest="data",
        help="get other input",
    )
    group.add_argument(
        "-o",
        "--output",
        const="output",
        action="append_const",
        dest="data",
        help="get model benchmark",
    )
    group.add_argument(
        "-s",
        "--source",
        const="source",
        action="append_const",
        dest="data",
        help="get source code for benchmark",
    )
    group.add_argument(
        "-d",
        "--docs",
        const="docs",
        action="append_const",
        dest="data",
        help="get corresponding user guide",
    )
    group.add_argument(
        "--extras",
        action="store_true",
        dest="extras",
        help="also get the extras, if any",
    )

    group = parser.add_argument_group("Download options")
    group.add_argument(
        "--yes",
        default=True,
        action="store_false",
        dest="ask",
        help="Don't ask before start downloading/unpacking",
    )
    group.add_argument(
        "--outpath",
        action="store",
        type=Path,
        dest="outpath",
        help="Override output (dataset) directory",
    )
    group.add_argument(
        "--tmppath",
        action="store",
        type=Path,
        dest="tmppath",
        help="Override temporary (download) directory",
    )
    group.add_argument(
        "--cleanup",
        default=False,
        action="store_true",
        dest="cleanup",
        help="Remove ALL temporary (download) files",
    )

    opts = parser.parse_args(args)
    logging.basicConfig(level=opts.log_level, format="%(message)s")

    if opts.tag is None and opts.year is None:
        opts.tag = DEFAULT.releases.split(",")[-1]
    if opts.data is None:
        opts.data = ["meteo", "input", "output", "source", "docs"]
    if opts.extras:
        opts.data.append("extra")
    if opts.outpath:
        assert isinstance(opts.outpath, Path)
        DEFAULT.output = opts.outpath
    if opts.tmppath:
        assert isinstance(opts.tmppath, Path)
        DEFAULT.downloads = opts.tmppath

    return opts


def user_consent(question, ask: bool = True, default: str = "yes"):
    """Ask user for confirmation"""

    # format question
    qq = {"y": "%s [Y]/n:", "yes": "%s [Y]/n:", "n": "%s y/[N]:", "no": "%s y/[N]:"}
    try:
        question = qq[default.lower()] % question
    except KeyError:
        logging.error(f"Unsupported option: user_consent(..,default='{default}')")

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
def print_progress(
    iteration, total: int, prefix: str = "", suffix: str = "", decimals: int = 2, length: int = 100
):
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


def file_size(value: float):
    """
    Human readable file size (K/M/G/T bytes)
      derived from http://stackoverflow.com/a/1094933/2576368
    """
    num = max(value, 0)
    for unit in "B", "K", "M", "G":
        if num < 1024.0:
            return f"{num:.1f}{unit}"
        num /= 1024.0
    return f"{num:.1f}T"


class DataPoint:
    """Info and retrieval/check methods"""

    def __init__(
        self,
        release: str | int,
        key: str,
        year: str | int,
        model: str,
        src: str,
        dst: str | Path = "",
        byteSize: str | int = 0,
        md5sum: str | None = None,
    ):
        """Initialize object"""

        self.release: int = int(release)  # release date (YYYYMM)
        self.key: str = key
        self.tag: str = {  # eg 'rv4_8source'
            "other": "{KEY}{REL}",
            "meteo": "{KEY}{YEAR}",
        }.get(self.key, "{MOD}{KEY}")

        self.year: int | None
        try:
            self.year = int(year)  # met-year
        except ValueError:
            self.year = None

        self.model: str | None
        try:
            self.model = str(model)  # model version (or meteo domain)
        except ValueError:
            self.model = None

        self.src: str  # single source url/file
        self.dst: Path  # path for uncompressed self.dst
        self.size: int = int(byteSize)  # self.src file size [bytes]
        self.md5sum = md5sum  # self.src checksum

        # replace keywords
        kwargs: dict[str, str | int | None] = _CONST.copy()  # type:ignore[assignment]
        kwargs.update(REL=self.release, KEY=self.key, YEAR=self.year, MOD=self.model)
        self.tag = self.tag.format_map(kwargs)
        kwargs.update(TAG=self.tag)
        self.src = src.format_map(kwargs)
        if isinstance(dst, Path):
            self.dst = dst
        elif not dst:
            assert isinstance(DEFAULT.downloads, Path)
            self.dst = DEFAULT.downloads / self.tag / Path(self.src).name
        else:
            self.dst = Path(dst.format_map(kwargs))
        if not self.dst.name:
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

    def cleanup(self):
        """Remove (raw) downloads"""
        if not self.dst.is_file():
            return
        logging.info(f"Cleanup  {self}")

        try:
            self.dst.unlink()
            if self.dst.parent != "":
                self.dst.parent.rmdir()
        except OSError as error:
            if error.errno == 1:  # operation not permitted (permissions?)
                pass
            elif error.errno == 39:  # directory was not empty
                pass
            else:
                raise error

    def check(self, *, cleanup: bool = False, quiet: bool = False):
        """Check download against md5sum"""
        if not self.dst.is_file():
            return False
        if not quiet:
            logging.info(f"Check    {self}")

        if self.md5sum != hashlib.md5(self.dst.read_bytes()).hexdigest():
            logging.info(f"  md5 /= {self.md5sum}")
            if cleanup:  # remove broken file
                self.cleanup()
            return False
        return True

    def download(self):
        """derived from http://stackoverflow.com/a/22776/2576368"""
        from urllib.request import HTTPError, urlopen

        # check if file exists/md5sum, remove file if fail md5sum
        if self.check(cleanup=True, quiet=True):
            return

        logging.info(f"Download {self}")

        try:
            url = urlopen(self.src)
        except HTTPError as e:
            logging.debug(e)
            logging.error(f"Could not download {self.src}")
            sys.exit(-1)

        self.dst.parent.mkdir(parents=True, exist_ok=True)
        with url, self.dst.open("wb") as file:
            block = 1_024 * 8  #   8K
            big_file = self.size > 1_048_576 * 128  # 128M
            n = 0
            while True:
                buff = url.read(block)
                if not buff:
                    break
                file.write(buff)
                n += len(buff)
                if big_file and (n == self.size or n % 1_048_576 == 0):  # 1M
                    print_progress(n, self.size, length=50)

    def unpack(self, inspect: bool = False, output: Path = DEFAULT.output):
        """Unpack download"""
        if not self.check():
            return

        if tarfile.is_tarfile(self.dst):
            logging.info(f"Untar    {self}")
            with closing(tarfile.open(self.dst, "r")) as file:
                if user_consent("  See the contents first?", inspect, "no"):
                    file.list(verbose=logging.getLogger().level <= logging.INFO)
                    if not user_consent("  Do you wish to proceed?", inspect):
                        logging.info("OK, skipping file")
                        return

                output.parent.mkdir(parents=True, exist_ok=True)
                try:
                    file.extractall(output)
                except EOFError as error:
                    logging.error(f"  Failed unpack '{self.dst}':\n    {error}.")
                    if not user_consent("    Do you wish to continue?", inspect, "yes"):
                        sys.exit(-1)

        elif self.dst.name == "pdf":
            assert self.model is not None
            outfile = output / self.key / self.model / "emep-ctm.pdf"
            outfile.parent.mkdir(parents=True, exist_ok=True)
            self.dst.rename(outfile)
        else:
            outfile = output / self.key / self.dst.name
            logging.info(f"Copy     {outfile}")
            outfile.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(self.dst, outfile)


class DataSet(NamedTuple):
    """Info and retrieval/check methods"""

    """release date (YYYYMM)"""
    release: int

    meteo: set[DataPoint]
    input: set[DataPoint]
    output: set[DataPoint]
    source: set[DataPoint]
    docs: set[DataPoint]
    extra: set[DataPoint]

    @property
    def tag(self) -> str:
        """revision tag"""
        for src in self.source:
            if src.model is not None:
                return src.model
        else:
            raise ValueError("no sources in DataSet")

    @property
    def year(self) -> int:
        """ "met year"""
        for met in self.meteo:
            if met.year is not None:
                return met.year
        else:
            raise ValueError("no meteo in DataSet")

    def __str__(self):
        return f"{self.tag:>8} (release:{self.release}, meteo:{self.year})"


def read_catalog(filename: Path) -> Iterator[DataSet]:
    """releases read from catalog csv-file"""
    from csv import reader

    # download catalog if file not found
    if not filename.is_file():
        DataPoint(0, "catalog", 0, "", DEFAULT.remote, filename).download()

    # catalog file should be present at this point
    with open(filename) as file:
        csv = reader(file, delimiter=",")
        next(csv)  # skip header
        catalog: list[DataPoint] = []  # list all src files (1 file per row)

        logging.debug(f"Reading {filename}")
        for row in csv:
            if not row or not any(row):
                continue  # skip empty lines
            try:
                catalog.append(DataPoint(*row))
                logging.debug(f"  {catalog[-1]}")
            except Exception:
                logging.error(f"Failed to parse ({filename}): {row}")
                logging.error(DataPoint(*row))
                raise
        logging.debug(f"{filename} read(srcs:{len(catalog)})")

    logging.debug(f"Indexing {filename}")

    dataset: defaultdict[str, set[DataPoint]]
    for rel in {c.release for c in catalog}:
        dataset = defaultdict(set)
        for c in catalog:
            if c.release == rel:
                dataset[c.key].add(c)

        ds = DataSet(
            rel,
            meteo=dataset["meteo"],
            input=dataset["input"],
            output=dataset["output"],
            source=dataset["source"],
            docs=dataset["docs"],
            extra=dataset["extra"],
        )
        logging.debug(f"  {ds}")
        yield ds


def main(
    path: Path,
    *,
    tag: str | None,
    year: str | None,
    data: list[str],
    domain: str | None,
    ask: bool,
    cleanup: bool,
) -> None:
    """Command line function"""
    if not set(data).issubset(DataSet._fields):
        raise ValueError(f"unsupported {set(data) - set(DataSet._fields)}")

    def get_datasets(
        catalog: Collection[DataSet], attr: str, target: str | int
    ) -> Iterator[DataSet]:
        """search catalog for tag|status|year DataSet"""
        logging.debug(f"Searching {attr}(s):{target}")

        for v in catalog:
            if getattr(v, attr) == target:
                yield v

    def get_downloads(catalog: Collection[DataSet], attr: str, target: str) -> Iterator[DataPoint]:
        """downloads for tag|status|year DataSet"""
        logging.debug(f"Searching datasets:{data}")

        for ds in get_datasets(catalog, attr, target):
            for key in data:
                # do not try to download 0-size files
                dataset = (x for x in getattr(ds, key) if x.size > 0)
                if key == "meteo" and domain is not None:
                    dataset = (x for x in dataset if x.model == domain)

                total = 0
                for x in dataset:
                    yield x
                    total += x.size
                if total:
                    logging.info(f"Queue download: {file_size(total):>6} {key:>6} {ds.tag}")

    # list files to download
    catalog = tuple(read_catalog(path))
    downloads: set[DataPoint] = set()  # files to download
    if tag is not None:
        aux = set(get_downloads(catalog, "tag", tag))
        if not aux:
            logging.warning(f"No datasets found for --revision={tag}")
        downloads.update(aux)
    if year is not None:
        aux = set(get_downloads(catalog, "year", year))
        if not aux:
            logging.warning(f"No datasets found for --year={year}")
        downloads.update(aux)

    # list file_sizes to download
    total = file_size(sum(x.size for x in downloads))
    logging.info(f"Queue download: {total:>6} Total")
    if not user_consent("Do you wish to proceed?", ask):
        logging.info("OK, bye")
        sys.exit(0)

    # download files
    for x in downloads:
        x.download()
        x.unpack(ask)
        if cleanup:
            x.cleanup()


if __name__ == "__main__":
    opts = parse_arguments(sys.argv[1:] or ["--help"])
    main(
        opts.catalog,
        tag=opts.tag,
        year=opts.year,
        data=opts.data,
        domain=opts.domain,
        ask=opts.ask,
        cleanup=opts.cleanup,
    )
