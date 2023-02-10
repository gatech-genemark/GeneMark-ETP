#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2019, Georgia Institute of Technology, USA
#
# Print RNA introns which are in conflict (partially overlap) with GeneMark
# predicted introns. Only GeneMark introns which do not coincide with any RNA
# intron are considered. The RNA alternatives still need to agree with original
# GeneMark prediction in the following ways: 1. The alternative RNA intron lies
# within region defined by outer boundaries of exons adjacent to original
# prediction 2. The alternative RNA intron does not change the reading frame
# oforiginal prediction. 3. The alternative RNA intron does not introduce any
# new stop codons in the prediction.
# ==============================================================


import argparse
import csv
import re
import sys
import tempfile
import os
import subprocess


def extractFeatureGtf(text, feature):
    regex = feature + ' "([^"]+)"'
    if re.search(regex, text):
        return re.search(regex, text).groups()[0]
    else:
        return None


class Intron:
    """Storage for introns
       Informtation about phase and coordinates of bordering
       exons is either loaded from input directly (in case
       the input file was already processed by this script
       before) or it is added to introns based on neighboring
       exon in gtf file.
    """
    def __init__(self, row, prevExon):
        self.row = row

        self.left = int(row[3])
        self.right = int(row[4])
        self.leftBoundary = extractFeatureGtf(self.row[8], "leftBoundary")
        self.rightBoundary = extractFeatureGtf(self.row[8], "rightBoundary")
        self.leftPhase = extractFeatureGtf(self.row[8], "leftPhase")
        self.rightPhase = extractFeatureGtf(self.row[8], "rightPhase")

        if prevExon:
            assert int(prevExon[4]) + 1 == int(row[3])
            self.leftBoundary = int(prevExon[3])
            self.leftPhase = prevExon[7]
            self.__addToColumn(self.leftBoundary, "leftBoundary")
            self.__addToColumn(self.leftPhase, "leftPhase")

    def getSignature(self):
        return self.row[0] + "_" + self.row[3] + "_" + self.row[4] + self.row[6]

    def getScore(self):
        score = int(self.row[5])
        if score == 0:
            score = 1
        return score

    def print(self):
        print("\t".join(self.row))

    def printToFile(self, file):
        file.write("\t".join(self.row) + "\n")

    def __addToColumn(self, variable, name):
        if variable:
            self.row[8] += " " + name + " \"" + str(variable) + "\";"

    def __eq__(self, other):
        return self.getSignature() == other.getSignature()

    def addExonBefore(self, exon):
        pass

    def addExonAfter(self, exon):
        assert int(self.row[4]) + 1 == int(exon[3])
        self.rightBoundary = int(exon[4])
        self.rightPhase = exon[7]
        self.__addToColumn(self.rightBoundary, "rightBoundary")
        self.__addToColumn(self.rightPhase, "rightPhase")

    def checkAlternative(self, alternative, args):
        # I tried to apply these rules in different order and the result is
        # that all these rules are important for High-quality filtering of
        # false positives, none is redundant.

        if int(self.leftBoundary) >= alternative.left or \
           int(self.rightBoundary) <= alternative.right:
            return False

        if not args.checkFrameAndORF:
            if args.conflictingPredictions:
                file = open(args.conflictingPredictions, "a")
                self.printToFile(file)
                file.close()
            return True

        if not self.checkFrame(alternative):
            return False

        strand = self.row[6]
        left = getSequence(args.genome, self.row[0], self.leftBoundary,
                           int(alternative.left) - 1, strand)
        right = getSequence(args.genome, self.row[0], int(alternative.right) + 1,
                            self.rightBoundary, strand)

        newExon = None
        if strand == '+':
            newExon = left + right
        else:
            newExon = right + left

        phase = self.leftPhase
        if strand == '-':
            phase = self.rightPhase

        if checkORF(newExon, int(phase)):
            if args.conflictingPredictions:
                file = open(args.conflictingPredictions, "a")
                self.printToFile(file)
                file.close()
            return True

        return False

    def checkFrame(self, alternative):
        """Check whether the change from this to alternative intron
        preserves original reading frame
        """

        # Check how many nucleotides are inserted (positive)
        # or deleted (negative)
        leftShift = alternative.left - self.left
        # Same check on the right -- logic is reversed because exon extension
        # happens when hint ends before ab initio end
        rightShift = self.right - alternative.right

        insertionLength = leftShift + rightShift
        # If the total extension/deletion of exon is divisible by 3,
        # it does not corrupt the reading frame
        if abs(insertionLength) % 3 == 0:
            return True
        return False


def getSequence(genome, chrom, start, end, strand):
    start = str(start)
    end = str(end)

    complement = ""
    if strand == '-':
        complement = " -i "

    cmd = "samtools faidx " + complement + genome + " " + chrom + ":" + start\
          + "-" + end + " | tail -n +2 | tr -d \"\n\""

    return subprocess.check_output(cmd, shell=True).decode()


def checkORF(seq, phase):
    """Check if the sequence contains a stop codon
       return True if check passes
    """
    length = len(seq)
    i = phase
    while i + 3 < length:
        if checkStop(seq[i:i + 3]):
            return False
        i += 3
    return True


def checkStop(seq):
    seq = seq.upper()
    if seq == "TAA" or seq == "TGA" or seq == "TAG":
        return True
    return False


def loadIntrons(inputFile):
    introns = []
    foundIntron = False
    prevCDS = None
    for row in csv.reader(open(inputFile), delimiter='\t'):
        if row[2].lower() == "intron":
            introns.append(Intron(row, prevCDS))
            foundIntron = True

        if row[2].lower() == "cds":
            if foundIntron:
                introns[-1].addExonAfter(row)
                foundIntron = False
            prevCDS = row

    return introns


def removeMatchingIntrons(candidates, hints):
    """Discard introns matching between genemark and rnaIntrons
    Such introns are already in agreement, thus not candidates
    for fixing.
    """
    mismatchingCandidates = []
    mismatchingHints = []

    matchingSignatures = set()
    candidateSignatures = set()
    for candidate in candidates:
        candidateSignatures.add(candidate.getSignature())

    for hint in hints:
        if not hint.getSignature() in candidateSignatures:
            mismatchingHints.append(hint)
        else:
            matchingSignatures.add(hint.getSignature())

    for candidate in candidates:
        if not candidate.getSignature() in matchingSignatures:
            mismatchingCandidates.append(candidate)

    return mismatchingCandidates, mismatchingHints


def filterHints(hints, score):
    """Discard introns with low coverage scores
    """
    filteredHints = []

    for candidate in hints:
        if candidate.getScore() >= score:
            filteredHints.append(candidate)

    return filteredHints


def printIntrons(introns, out):
    for intron in introns:
        intron.printToFile(out)


def processIntersectionResult(intersectOut, args):
    for row in csv.reader(open(intersectOut), delimiter='\t'):
        candidate = Intron(row[0:9], None)
        alternative = Intron(row[9:18], None)

        if candidate.checkAlternative(alternative, args):
            alternative.print()


def printGoodAlternatives(candidates, hints, args):
    """Identify candidates for fixing -- ab initio introns overlapped
    by RNA seq introns.

    For each ab initio candidate, return a list of overlapping introns
    """

    canOut = tempfile.NamedTemporaryFile(dir=".", delete=False, mode="w")
    printIntrons(candidates, canOut)
    canOut.close()

    hintsOut = tempfile.NamedTemporaryFile(dir=".", delete=False, mode="w")
    printIntrons(hints, hintsOut)
    hintsOut.close()

    intersectOut = tempfile.NamedTemporaryFile(dir=".", delete=False).name
    cmd = "bedtools intersect -a " + canOut.name + " -b " + \
        hintsOut.name + " -wa -wb -s > " + intersectOut
    subprocess.call(cmd, shell=True)

    os.remove(canOut.name)
    os.remove(hintsOut.name)

    processIntersectionResult(intersectOut, args)
    os.remove(intersectOut)


def main():

    args = parseCmd()
    candidates = loadIntrons(args.genemark)
    hints = loadIntrons(args.rnaIntrons)

    hints = filterHints(hints, args.minIntronScore)

    candidates, hints = removeMatchingIntrons(candidates, hints)

    if args.otherIntrons:
        otherIntrons = loadIntrons(args.otherIntrons)
        candidates, _ = removeMatchingIntrons(candidates, otherIntrons)

    printGoodAlternatives(candidates, hints, args)


def parseCmd():

    description = "Print RNA introns which are in conflict \
        (partially overlap) with GeneMark predicted introns. \
        Only GeneMark introns which do not coincide with any RNA \
        intron are considered. The RNA alternatives still need to \
        agree with original GeneMark prediction in the following \
        ways: 1. The alternative RNA intron lies within region \
        defined by outer boundaries of exons adjacent to original \
        prediction 2. The alternative RNA intron does not change \
        the reading frame oforiginal prediction. 3. The alternative \
        RNA intron does not introduce any new stop codons in the prediction.\
        Filters 2. and 3. are optional."

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('genemark', metavar='genemark.gtf', type=str,
                        help='GeneMark gtf predictions')

    parser.add_argument('rnaIntrons', metavar='rnaIntrons.gff', type=str,
                        help='Introns predicted from RNA mapping')

    parser.add_argument('--otherIntrons', metavar='otherIntrons.gff', type=str,
                        help='Other, non-RNA, protein hints (for example \
        protein hints. Ab initio predictions supported by these other hints \
        will not be modified.')

    parser.add_argument('--genome', metavar='genome.fasta', type=str,
                        help='Genomic sequence')

    parser.add_argument('--checkFrameAndORF',  default=False, action='store_true',
                        help='Check that alternative introns do not change a \
                        reading frame and do not introduce any new stop codons.')

    parser.add_argument('--conflictingPredictions', type=str, help='Print \
                        conflicting GeneMark predictions to this file.')

    parser.add_argument('--minIntronScore', type=int, help='Ignore RNA \
                        introns with score lower than this.', default=1)

    args = parser.parse_args()

    if args.checkFrameAndORF and not args.genome:
        sys.exit("Error, --genome option must be specified when " \
                 "--checkFrameAndORF is used.")

    return args


if __name__ == '__main__':
    main()
