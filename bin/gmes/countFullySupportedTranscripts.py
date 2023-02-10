#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2020, Georgia Institute of Technology, USA
#
# Count how many transcripts in the input set are fully supported by evidence
# hints. All introns in a transcript must be supported by external evidence. In
# case of single-exon transcripts, both start and stop codon must be supported
# by external evidence. On top of these criteria, transcripts must be complete.
#
# ==============================================================


import argparse
import predictionAnalysis as analysis


def main():
    args = parseCmd()

    prediction = analysis.PredictionAnalysis(args.prediction, args.hints)
    print(prediction.getFullySupportedTranscriptCount())


def parseCmd():

    parser = argparse.ArgumentParser(description='Count how many transcripts \
        in the input set are fully supported by evidence hints. All introns \
        in a transcript must be supported by external evidence. In case of \
        single-exon transcripts, both start and stop codon must be supported \
        by external evidence. On top of these criteria, transcripts must be \
        complete.')

    parser.add_argument('prediction', metavar='prediction.gtf', type=str,
                        help='Prediction file.')

    parser.add_argument('hints', metavar='hints.gff', type=str,
                        help='File with external hints.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
