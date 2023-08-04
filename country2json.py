#! /usr/bin/env python3

"""
read the models Country_mod.f90 and extract the countries in a json file
"""

import json
import logging
import os
import re
import sys

def clean_string(txt: str) -> str:
    txt = txt.strip()
    if txt.startswith('"'):
        txt = txt.replace('"', '')
    elif txt.startswith("'"):
        txt = txt.replace("'", "")
    txt = txt.strip()
    return txt

def true_false_str(txt: str) -> bool:
    txt = txt.strip().lower()
    if txt.startswith('f'):
        return False
    if txt.startswith('t'):
        return True
    return int(txt) != 0


def country2json(filename):
    """
    countries in filename are defined either as
        cc( code,gains, icode ,is_sea, timezone,timezone_h,timefac,cname)
    or as
        add_Country(ix,ic,code,gains,icode,is_sea,timezone,timezone_h,timefac,cname)
    """
    countries = []
    with open(filename, 'rt') as fh:
        cc_pattern = re.compile(r'[^\w]cc\s*\(([^)]*)\)')
        add_country_pattern = re.compile(r'add_Country\s*\(([^)]*)\)')
        code_dubs = set()
        icode_dubs = set()
        for line in fh:
            if line.strip().startswith('!'):
                continue
            if m := cc_pattern.search(line):
                #print(m[0])
                row = m[1].split(',')
            elif m := add_country_pattern.search(line):
                row = m[1].split(',')[2:]
            else:
                continue

            if len(row) < 8:
                print(f"missing elements in line: {line}, {row}")
                continue
            if (row[1].strip() == 'gains'):
                # code-line, not definition-line
                continue

            try:
                country = {
                    "cname": clean_string(row[7]),
                    "code": clean_string(row[0]),
                    "gains": clean_string(row[1]),
                    "icode": int(row[2]),
                    "is_sea": true_false_str(row[3]),
                    "timezone": int(row[4]),
                    "timezone_h": int(row[5]),
                    "timefac": int(row[6]),
                }
            except Exception as e:
                logging.warning(f"{e} in line: {line}")
                continue
            countries.append(country)

            if country["code"] in code_dubs:
                logging.warning(f"Duplicated code: {country['code']}")
            else:
                code_dubs.add(country["code"])

            if country["icode"] in icode_dubs:
                logging.warning(f"Duplicated icode: {country['icode']} for {country['code']}")
            else:
                icode_dubs.add(country["icode"])

    json.dump(countries, sys.stdout)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"usage: {sys.argv[0]} Country_ml.f90", file=sys.stderr)
        exit(1)
    else:
        country2json(sys.argv[1])

