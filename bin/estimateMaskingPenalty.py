#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2022, Georgia Institute of Technology, USA
#
# Estimate optimal repeat masking penalty
# ==============================================================


import argparse
import subprocess
import tempfile
import os
import csv
import sys
import shutil
import re
import logging


def temp(prefix, suffix):
    tmp = tempfile.NamedTemporaryFile("w", delete=False, dir=".",
                                      prefix=prefix, suffix=suffix)
    return tmp


def extractFeatureGtf(text, feature):
    regex = feature + ' "([^"]+)"'
    return re.search(regex, text).groups()[0]


def systemCall(cmd):
    if subprocess.call(["bash", "-c", cmd]) != 0:
        sys.exit('error: Program exited due to an ' +
                 'error in command: ' + cmd)


def checkOutput(cmd):
    return subprocess.check_output(["bash", "-c", cmd]).decode()


def loadGeneBorders(pred):
    borders = {}
    for row in csv.reader(open(pred), delimiter='\t'):
        if len(row) != 9:
            continue
        if row[2] != "CDS":
            continue
        geneID = extractFeatureGtf(row[8], "gene_id")
        c = row[0]
        if c not in borders:
            borders[c] = {}
        if geneID not in borders[c]:
            borders[c][geneID] = [int(row[3]), int(row[4])]
        else:
            borders[c][geneID][0] = min(borders[c][geneID][0], int(row[3]))
            borders[c][geneID][1] = max(borders[c][geneID][1], int(row[4]))
    return borders


def saveSingleChrom(seq, borders, margin, outF):
    for ID, gene in borders.items():
        s = gene[0] - 1 - margin
        e = gene[1] + margin

        if s < 0:
            s = 0
        # Too long won't cause an issue
        outF.write(">" + ID + "\n")
        outF.write(seq[s:e] + "\n")


def saveContigs(fastaFile, borders, margin):
    chrom = ''
    seq = ''
    fasta = open(fastaFile, "r")
    outF = temp("contigs", ".fasta")
    for line in fasta:
        if len(line) == 1:
            continue

        if line[0] == '>':
            if chrom != '':
                if chrom in borders:
                    saveSingleChrom(seq, borders[chrom], margin, outF)
            chrom = line.split()[0][1:]
            seq = ''
        else:
            seq += line.strip()

    if chrom in borders:
        saveSingleChrom(seq, borders[chrom], margin, outF)

    outF.close()
    return os.path.abspath(outF.name)


def rescaleReference(pred, borders, margin):
    outF = temp("reference", ".gtf")
    for row in csv.reader(open(pred), delimiter='\t'):
        if len(row) != 9:
            continue
        geneID = extractFeatureGtf(row[8], "gene_id")
        offset = borders[row[0]][geneID][0] - 1 - margin
        if offset < 0:
            offset = 0
        row[0] = geneID
        row[3] = str(int(row[3]) - offset)
        row[4] = str(int(row[4]) - offset)
        outF.write("\t".join(row) + "\n")
    outF.close()
    return os.path.abspath(outF.name)


class PenaltyEstimator():

    def __init__(self, args):
        self.referencePred = args.pred
        self.preparePredictionContigs(args.seq, args.pred, args.margin)
        self.mod = args.model
        self.args = args
        self.pGeneCache = {}
        self.pExonCache = {}

    def preparePredictionContigs(self, seq, pred, margin):
        borders = loadGeneBorders(pred)
        self.reference = rescaleReference(pred, borders, margin)
        self.contigs = saveContigs(seq, borders, margin)

    def estimate(self):
        self.baselinePenalty = self.argMaxReliable(self.args.penaltyMin,
                                                   self.args.penaltyMax,
                                                   self.args.startingStep)
        logging.info("Selected baseline penalty for the maximum # of " +
                     "correct reliable predictions: " +
                     str(self.baselinePenalty))

        maxReliable = self.pExonCache[str(self.baselinePenalty)]
        # Select the maximum value of a penalty which preserves the specified
        # fraction of baseline reliable exons out of the already computed vals
        minP = self.baselinePenalty
        p = round(self.baselinePenalty + self.args.minStep, 2)
        while p <= self.args.penaltyMax:
            if str(p) in self.pExonCache:
                if self.pExonCache[str(p)] / maxReliable >= self.args.minRelFrac:
                    minP = p
            p = round(p + self.args.minStep, 2)

        # Select the next pre-computed value as an upper bound for the search
        nextP = self.args.penaltyMax
        p = round(minP + self.args.minStep, 2)
        while p <= self.args.penaltyMax:
            if str(p) in self.pExonCache:
                nextP = p
                break
            p = round(p + self.args.minStep, 2)

        selected = self.largestAllowed(minP, nextP, self.args.startingStep)
        logging.info("Masking penalty was set to " + str(selected))
        return selected

    def largestAllowed(self, minP, maxP, step):
        if step < self.args.minStep:
            sys.exit("Error: Unexpected step during largest p. estimation.")
        maxReliable = self.pExonCache[str(self.baselinePenalty)]
        bestPenalty = minP

        p = minP
        while p <= maxP:
            correct = self.predictWithPenalty(p)
            if correct / maxReliable >= self.args.minRelFrac:
                bestPenalty = p
            p = round(p + step, 2)

        if step / 2 < self.args.minStep - 1e-10:
            return bestPenalty

        if bestPenalty + step < maxP:
            maxP = round(bestPenalty + step, 2)

        return self.largestAllowed(bestPenalty, maxP, round(step / 2, 2))

    def argMaxReliable(self, minP, maxP, step):
        logging.info("Finding masking penalty maximizing the number of " +
                     "correctly predicted reliable exons in range from " +
                     str(minP) + " to " + str(maxP) + " with step " +
                     str(step))

        if step < self.args.minStep:
            sys.exit("Error: Unexpected step during max estimation.")

        maxR = 0
        bestPenalty = minP

        p = minP
        while p <= maxP:
            correct = self.predictWithPenalty(p)
            if correct > maxR:
                maxR = correct
                bestPenalty = p
            p = round(p + step, 2)

        if step / 2 < self.args.minStep - 1e-10:
            return bestPenalty

        # Put the current penalty in the middle of the search space
        if bestPenalty - step > minP:
            minP = round(bestPenalty - step, 2)
        if bestPenalty + step < maxP:
            maxP = round(bestPenalty + step, 2)

        return self.argMaxReliable(minP, maxP, round(step / 2, 2))

    def predictWithPenalty(self, penalty):
        if str(penalty) in self.pExonCache:
            return self.pExonCache[str(penalty)]

        logging.info("Running prediction with masking penalty = " +
                     str(penalty))

        dirpath = tempfile.mkdtemp(prefix="gmes", dir='.')
        os.chdir(dirpath)

        systemCall(self.args.GMES_PATH + "gmes_petap.pl --seq " + 
                   self.contigs + " --soft_mask 1000 --max_mask 40000 " +
                   " --predict_with " + self.mod + " --cores " +
                   self.args.threads + " --mask_penalty " + str(penalty))

        compareOut = checkOutput(self.args.binDir +
                                 "/compare_intervals_exact.pl --f1 "
                                 + self.reference +
                                 " --f2 genemark.gtf --gene | head -2 |" +
                                 "tail -1 | cut -f2")

        compareOutExon = checkOutput(self.args.binDir +
                                     "/compare_intervals_exact.pl --f1 "
                                     + self.reference +
                                     " --f2 genemark.gtf | head -2 |" +
                                     "tail -1 | cut -f2")

        os.chdir('..')
        shutil.rmtree(dirpath)
        self.pGeneCache[str(penalty)] = int(compareOut)
        self.pExonCache[str(penalty)] = int(compareOutExon)
        return self.pExonCache[str(penalty)]

    def scan(self):
        p = self.args.penaltyMin
        print("penalty\tgenes\texons")
        while p <= self.args.penaltyMax:
            self.predictWithPenalty(p)
            print("\t".join([str(p),
                             str(self.pGeneCache[str(p)]),
                             str(self.pExonCache[str(p)])]))
            p = round(p + self.args.startingStep, 2)

    def cleanup(self):
        os.remove(self.reference)
        os.remove(self.contigs)


def main():
    args = parseCmd()
    logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s',
                        datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)
    estimator = PenaltyEstimator(args)
    if not args.scan:
        print(estimator.estimate())
        estimator.cleanup()
    else:
        estimator.scan()


def parseCmd():

    parser = argparse.ArgumentParser(description='Estimate optimal repeat \
                                     masking penalty.')

    parser.add_argument('pred', metavar="reliable.gtf", type=str,
                        help='Reliable gene predictions which serve as a \
                        ground truth during penalty estimation.')

    parser.add_argument('seq', metavar="sequence.sofmasked.fasta", type=str,
                        help='Softmasked genomic sequence')

    parser.add_argument('model', metavar="gmhmm.mod", type=str,
                        help='GM HMM prediction model.')

    parser.add_argument('--margin', type=int, default=1000,
                        help='Sequence margin around reliable genes')

    parser.add_argument('--threads', type=int, default=1,
                        help='Number of CPU threads to use.')

    parser.add_argument('--penaltyMin', type=float, default=0,
                        help='Minimum value for the output penalty.')

    parser.add_argument('--penaltyMax', type=float, default=0.2,
                        help='Maximum value for the output penalty.')

    parser.add_argument('--startingStep', type=float, default=0.04,
                        help='Initial search step.')

    parser.add_argument('--minStep', type=float, default=0.01,
                        help='Smallest search step.')

    parser.add_argument('--minRelFrac', type=float, default=0.998,
                        help='Find the maximum value preserving at least \
                        MINRELFRAC of the maximum number of reliable exons \
                        correctly predicted.')

    parser.add_argument('--GMES_PATH', type=str, default='',
                        help='Path to folder with gmes_petap.pl')

    parser.add_argument('--scan', action='store_true', help='Scan all penalty\
        values from min to max with the specified STARTINGSTEP and save the\
        results. This DOES NOT estimate any masking penalty.')

    args = parser.parse_args()

    if args.minStep < 0.01:
        sys.exit("error: the minimum specified step size is too small. " +
                 "Minimum allowed step size is 0.01.\n")

    if args.minStep > args.startingStep:
        sys.exit("error: the initial search step must be >= minStep\n")

    if args.GMES_PATH != '':
        args.GMES_PATH = os.path.abspath(args.GMES_PATH) + '/'
    args.model = os.path.abspath(args.model)
    args.threads = str(args.threads)
    args.binDir = os.path.abspath(os.path.dirname(__file__))

    return args


if __name__ == '__main__':
    main()
