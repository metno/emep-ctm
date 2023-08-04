#! /usr/bin/env python3

import argparse
import csv
import datetime
import json
import logging
import pathlib
import sys
import netCDF4

import numpy

def country_definitions():
    countryfile = pathlib.Path(__file__).parent.resolve() / 'countries.json'
    if not countryfile.exists():
        logging.error("countries.json not in the directory of this program.\n Please copy or generate using 'country2json.py Country_ml.f90'")
        exit(1)
    with countryfile.open() as fh:
        ccs = json.load(fh)
    cc_by_code = dict((cc['code'], cc) for cc in ccs)
    return cc_by_code

#GNFR N14, ignoring NT=natural
class Emissions_TXT2NC_CV():
    lon_start = -29.95
    lon_size = 1200
    lon_res = 0.1
    lat_start = 30.05
    lat_size = 520
    lat_res = 0.1
    sectors = [
        "A_PublicPower",
        "B_Industry",
        "C_OtherStationaryComb",
        "D_Fugitive",
        "E_Solvents",
        "F_RoadTransport",
        "G_Shipping",
        "H_Aviation",
        "I_Offroad",
        "J_Waste",
        "K_AgriLivestock",
        "L_AgriOther",
        "M_Other"
    ]
    sectormap = {f"N14 {x}": i for i,x in enumerate(sectors)}

    comps = {
        "SOx": {
            "model": "sox",
            "species": "sox",
            "units": "tonnes/year",
            "molecular_weight": 64.,
            "molecular_weight_units": "g mole-1",
        },
        "NOx": {
            "model": "nox",
            "species": "nox",
            "units": "tonnes/year",
            "molecular_weight": 46.,
            "molecular_weight_units": "g mole-1",
        },
        "NH3": {
            "model": "nh3",
            "species": "nh3",
            "units": "tonnes/year",
            "molecular_weight": 17.,
            "molecular_weight_units": "g mole-1",
        },
        "NMVOC": {
            "model": "voc",
            "species": "voc",
            "units": "tonnes/year",
        },
        "CO": {
            "model": "co",
            "species": "co",
            "units": "tonnes/year",
            "molecular_weight": 28.,
            "molecular_weight_units": "g mole-1",
        },
        "PM10": {
            "model": "pm10",
            "species": "pm10",
            "units": "tonnes/year",
        },
        "PMcoarse": {
            "model": "pmco",
            "species": "pmco",
            "units": "tonnes/year",
        },
        "PM2_5": {
            "model": "pm25",
            "species": "pm25",
            "units": "tonnes/year",
        },
        #TBD ask agnes what to do with BC-files
        "BC": {
            "model": "ec",
            "species": "ec",
            "units": "tonnes/year",
        },
    }

    def __init__(self, year) -> None:
        self.year = year
        self.iso2 = country_definitions()
        self.countries = set()
        self.arrays = {}
        return

    def _lon2pos(self, lon: float) -> int:
        return round((lon - self.lon_start) * 10) # 0.1 deg starting at -29.95
    def _lat2pos(self, lat: float) -> int:
        return round((lat - self.lat_start) * 10) # 0.1 deg starting at 30.05

    def _array_name(self, pollutant: str, country: str=None) -> str:
        mapname = self.comps[pollutant]["model"]
        if not country is None:
            mapname += "_" + country
        return mapname

    def _getmap(self, pollutant: str, country: str=None) -> numpy.array:
        mapname = self._array_name(pollutant, country)
        if not mapname in self.arrays:
            if not country is None:
                self.arrays[mapname] = numpy.zeros((len(self.sectors),self.lat_size,self.lon_size), dtype='f4')
            else:
                self.arrays[mapname] = numpy.zeros((self.lat_size,self.lon_size), dtype='f4')
        return self.arrays[mapname]


    def _add_value(self, pollutant, country, sector, lon, lat, emis) -> None:
        ilon = self._lon2pos(lon)
        ilat = self._lat2pos(lat)
        isec = self.sectormap[sector]
        self.countries.add(country)
        self._getmap(pollutant)[ilat, ilon] += emis
        self._getmap(pollutant, country)[isec, ilat, ilon] += emis
        return

    def _create_data_var(self, nc: netCDF4.Dataset, pollutant: str, country: str=None) -> None:
        mapname = self._array_name(pollutant, country)
        if not mapname in self.arrays:
            print(f"no data for {mapname}")
        if not country is None:
            dims = ('time', 'sector', 'lat', 'lon')
            chunks = (1,13,40,40)
        else:
            dims = ('time', 'lat', 'lon')
            chunks = (1,40,40)
        v = nc.createVariable(mapname, 'f4',
                              dims,
                              zlib=True,
                              # compression='zlib', # use in newer netCDF4 versions
                              complevel=4,
                              chunksizes=chunks)
        for att, val in self.comps[pollutant].items():
            if att == 'model':
                continue
            if not (att == 'species' and country is None):
                v.setncattr(att, val)
        if not country is None:
            v.setncattr('country_ISO', country)
            v.setncattr('countrycode', numpy.array([self.iso2[country]['icode']], 'i4'))
        return


    def write_netcdf(self, filename) -> None:
        with netCDF4.Dataset(filename, 'w') as nc:
            nc.Conventions = "CF-1.6"
            args = " ".join(sys.argv)
            if len(args) > 100:
                args = args[:100] + "..."
            nc.history = f"{datetime.datetime.now():%Y-%m-%dT%H:%M:%S} creation with {args}"
            nc.projection = "lon lat"
            nc.periodicity = "yearly"
            nc.SECTORSq_NAME = "GNFR"
            nc.createDimension('lon', self.lon_size)
            nc.createDimension('lat', self.lat_size)
            nc.createDimension('sector', len(self.sectors))
            nc.createDimension('time', 1)

            v = nc.createVariable('lon', 'f4', ('lon'))
            v.standard_name = "longitude"
            v.units = "degrees_east"

            v = nc.createVariable('lat', 'f4', ('lat'))
            v.standard_name = "latitude"
            v.units = "degrees_north"

            v = nc.createVariable('sector', 'i4', ('sector'))
            v.long_name = "GNFR sector index"
            secs = numpy.arange(1, len(self.sectors)+1, dtype='i4')
            v.flag_values = secs
            v.flag_meanings = " ".join(self.sectors)

            v = nc.createVariable('time', 'i4', ('time'))
            v.units = f"days since {self.year}-06-30 00:00:00 +00:00"

            for comp in self.comps:
                compname = self.comps[comp]["model"]
                if compname in self.arrays:
                    for country in self.countries:
                        self._create_data_var(nc, comp, country)
                    self._create_data_var(nc, comp)

            # and the data
            nc['lon'][:] = (numpy.arange(0, self.lon_size, dtype='f4') * self.lon_res) + self.lon_start
            nc['lat'][:] = (numpy.arange(0, self.lat_size, dtype='f4') * self.lat_res) + self.lat_start
            nc['sector'][:] = secs
            nc['time'][:] = numpy.zeros(1)

            for varname, array in self.arrays.items():
                nc[varname][:] = array
        return


    def add_file(self, file) -> None:
        with open(file, newline='') as fh:
            csvreader = csv.DictReader(fh, fieldnames="ISO2;YEAR;SECTOR;POLLUTANT;LONGITUDE;LATITUDE;UNIT;EMISSION".split(';'), delimiter=';')
            for row in csvreader:
                if row["ISO2"].startswith("#"):
                    continue
                if not row["ISO2"] in self.iso2:
                    print(f"unknown country {row['ISO2']} in file: {file}", file=sys.stderr)
                    continue
                if int(row['YEAR']) != self.year:
                    print(f"wrong year {row['YEAR']} in file: {file}", file=sys.stderr)
                    continue
                if row['SECTOR'] == 'N14 National Total':
                    continue # not needed
                if not row['SECTOR'] in self.sectormap:
                    print(f"unknown sector {row['SECTOR']} in file: {file}", file=sys.stderr)
                    continue
                if not row['POLLUTANT'] in self.comps:
                    print(f"unknown pollutant {row['POLLUTANT']} in file: {file}", file=sys.stderr)
                    continue
                if row['UNIT'] != 'Mg':
                    print(f"unknown unit {row['UNIT']} in file: {file}", file=sys.stderr)
                    continue
                longitude = float(row['LONGITUDE'])
                latitude = float(row['LATITUDE'])
                emis = float(row['EMISSION'])
                self._add_value(row["POLLUTANT"], row["ISO2"], row["SECTOR"], longitude, latitude, emis)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='read CEIP emission files and convert them to netcdf CV (country-variable) format')
    parser.add_argument('inputfiles', nargs='+', help='all CEIP txt emission files for one year')
    parser.add_argument('--out', help='output netcdf file', required=True)
    parser.add_argument('--year', type=int, help='emission year to convert', required=True)
    args = parser.parse_args()


    emis2nc = Emissions_TXT2NC_CV(args.year)
#    print(emis2nc.iso2['ES'], file=sys.stderr)
    for f in args.inputfiles:
        print(f"reading {f}")
        emis2nc.add_file(f)
    emis2nc.write_netcdf(args.out)
