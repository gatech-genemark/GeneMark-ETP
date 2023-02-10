#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2020, Georgia Institute of Technology, USA
#
# Add "anchored" flag to all CDS in GeneMark predictions. Status "anchored
# 1_1" is added to fully anchored exons, partially anchored exons (only one
# splice site/start/stop codon anchored) are flagged with "anchored 1_0" or
# "anchored 0_1" flag.
# ==============================================================


import argparse
import csv
import re


def extractFeatureGtf(text, feature):
    regex = feature + ' "([^"]+)"'
    return re.search(regex, text).groups()[0]


def leftSignature(row):
    return row[0] + "_" + row[3] + "_" + row[6]


def rightSignature(row):
    return row[0] + "_" + row[4] + "_" + row[6]


def loadHints(gffFile):
    spliceSites = set()
    starts = set()
    stops = set()
    for row in csv.reader(open(gffFile), delimiter='\t'):
        if row[2].lower() == "intron":
            spliceSites.add(row[0] + "_" + str(int(row[3]) - 1) + "_" + row[6])
            spliceSites.add(row[0] + "_" + str(int(row[4]) + 1) + "_" + row[6])
        elif row[2].lower() == "start_codon":
            if row[6] == '+':
                starts.add(leftSignature(row))
            else:
                starts.add(rightSignature(row))
        elif row[2].lower() == "stop_codon":
            if row[6] == '+':
                stops.add(rightSignature(row))
            else:
                stops.add(leftSignature(row))

    return spliceSites, starts, stops


def flagInternal(row, spliceSites):
    left = "0"
    right = "0"

    if leftSignature(row) in spliceSites:
        left = "1"
    if rightSignature(row) in spliceSites:
        right = "1"

    return left, right


def flagInitial(row, spliceSites, starts):
    left = "0"
    right = "0"

    if row[6] == "+":
        if leftSignature(row) in starts:
            left = "1"
        if rightSignature(row) in spliceSites:
            right = "1"
    else:
        if leftSignature(row) in spliceSites:
            left = "1"
        if rightSignature(row) in starts:
            right = "1"

    return left, right


def flagTerminal(row, spliceSites, stops):
    left = "0"
    right = "0"

    if row[6] == "+":
        if rightSignature(row) in stops:
            right = "1"
        if leftSignature(row) in spliceSites:
            left = "1"
    else:
        if rightSignature(row) in spliceSites:
            right = "1"
        if leftSignature(row) in stops:
            left = "1"

    return left, right


def flagSingle(row, starts, stops):
    left = "0"
    right = "0"

    if row[6] == "+":
        if leftSignature(row) in starts:
            left = "1"
        if rightSignature(row) in stops:
            right = "1"
    else:
        if leftSignature(row) in stops:
            left = "1"
        if rightSignature(row) in starts:
            right = "1"

    return left, right


def flag(genemark, spliceSites, starts, stops):
    for row in csv.reader(open(genemark), delimiter='\t'):
        if row[2] != "CDS":
            print("\t".join(row))
            continue

        cdsType = extractFeatureGtf(row[8], "cds_type")

        left = "0"
        right = "0"

        if cdsType.lower() == "internal":
            left, right = flagInternal(row, spliceSites)
        elif cdsType.lower() == "initial":
            left, right = flagInitial(row, spliceSites, starts)
        elif cdsType.lower() == "terminal":
            left, right = flagTerminal(row, spliceSites, stops)
        elif cdsType.lower() == "single":
            left, right = flagTerminal(row, starts, stops)

        if left == "0" and right == "0":
            print("\t".join(row))
        else:
            flag = left + "_" + right
            print("\t".join(row) + ' anchored "' + flag + '";')


def main():
    args = parseCmd()
    spliceSites, starts, stops = loadHints(args.hints)
    flag(args.genemark, spliceSites, starts, stops)


def parseCmd():

    parser = argparse.ArgumentParser(description='Add "anchored" flag to all\
                                     CDS in GeneMark predictions. Status\
                                     "anchored 1_1" is added to fully anchored\
                                     exons, partially anchored exons (only one\
                                     splice site/start/stop codon anchored)\
                                     are flagged with "anchored 1_0" or \
                                     "anchored 0_1" flag.')

    parser.add_argument('genemark', metavar='genemark.gtf', type=str,
                        help='GeneMark prediction file.')

    parser.add_argument('hints', metavar='hints.gff', type=str,
                        help='Hints used to guide GeneMark training.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
