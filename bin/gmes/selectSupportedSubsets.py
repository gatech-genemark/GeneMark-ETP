#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2020, Georgia Institute of Technology, USA
#
# Split supplied gene set into different subsets based on the amount of support
# by external evidence. The script relies on the "anchored" attribute in the 
# gtf input. GeneMark predictions already come with this attribute when run in
# ET or EP+ mode. To manually add this attribute to gtf, use the
# flag_anchored_elements.py script.
# ==============================================================


import argparse
import statistics.predictionAnalysis as analysis


def main():
    args = parseCmd()

    prediction = analysis.PredictionAnalysis(args.genemark)
    prediction.saveFullySupportedGenes(args.fullSupport)
    prediction.saveAnySupportedGenes(args.anySupport)
    prediction.saveNoSupportedGenes(args.noSupport)


def parseCmd():

    parser = argparse.ArgumentParser(description='Split supplied gene set \
        into different subsets based on the amount of support by external \
        evidence. The script relies on the "anchored" attribute in the gtf \
        input. GeneMark predictions already come with this attribute when run \
        in ET or EP+ mode. To manually add this attribute to gtf, use the \
        flag_anchored_elements.py script.')

    parser.add_argument('genemark', metavar='genemark.gtf', type=str,
                        help='GeneMark prediction file.')
    parser.add_argument('--fullSupport', type=str,
                        help='Output genes fully supported by external \
        evidence to this file. All introns in a gene must be supported by \
        external evidence. In case of single-exon genes, both start and stop \
        codon must be supported by external evidence.')
    parser.add_argument('--anySupport', type=str,
                        help='Output genes with any external support to this \
        file. At least one intron, start or stop codon of a predicted gene \
        must be supported.')
    parser.add_argument('--noSupport', type=str,
                        help='Output genes with no external support to this \
        file')

    return parser.parse_args()


if __name__ == '__main__':
    main()
