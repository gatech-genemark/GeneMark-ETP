#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2020, Georgia Institute of Technology, USA
#
# Generate GeneMark annotation report
#
# The code for outlier aware histogram comes from
# https://github.com/bdoughty/outlier-aware-histogram
# ==============================================================


import argparse
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import statistics.predictionAnalysis as analysis
from matplotlib.backends.backend_pdf import PdfPages
import statistics.motifs as motifs


def printBasicStatistics(prediction, report):
    fig = plt.figure()
    fig.text(0.5, 0.9, 'GeneMark Prediction Report', ha='center', va='center',
             size=20)

    fig.text(0.1, 0.8, 'General Statistics', ha='left', va='top', size=15)

    geneCount = prediction.getGeneCount()
    anySupport = prediction.getAnySupportedGeneCount()
    fullSupport = prediction.getFullySupportedGeneCount()
    noSupport = geneCount - anySupport
    complete = prediction.getCompleteCount()
    partial = prediction.getIncompleteCount()

    text = 'Gene count: ' + str(prediction.getGeneCount()) + "\n"
    text += '    Single-exon genes: ' + str(prediction.getSingleGeneCount()) + "\n"
    text += '    Multi-exon genes: ' + str(prediction.getMultiGeneCount()) + "\n"
    text += "\n"
    text += 'Introns per gene: ' + str(round(prediction.getIntronsPerGene(), 2)) + "\n"
    text += 'Introns per multi-exon gene: ' + str(round(prediction.getIntronsPerMultiGene(), 2)) + "\n"
    text += "\n"
    text += "Genes fully supported by external evidence: " + str(fullSupport) + \
            " (" + str(round(100 * fullSupport / geneCount, 2)) + "%)\n"
    text += "Genes partially supported by external evidence: " + str(anySupport) + \
            " (" + str(round(100 * anySupport / geneCount, 2)) + "%)\n"
    text += "Genes unsupported by any external evidence: " + str(noSupport) + \
            " (" + str(round(100 * noSupport / geneCount, 2)) + "%)\n"
    text += "\n"
    text += "Complete genes: " + str(complete) + \
            " (" + str(round(100 * complete / geneCount, 2)) + "%)\n"
    text += "Partial genes: " + str(partial) + \
            " (" + str(round(100 * partial / geneCount, 2)) + "%)\n"

    fig.text(0.1, 0.72, text, ha='left', va='top', size=9, linespacing=1.5)

    report.savefig(fig)


def getBounds(data, zScore):
    std = np.std(data)
    median = np.median(data)
    return (median + zScore * std)


def histogram(data, report, title, xlabel, zScore,
              minimum=-1, maximum=-1, bins=-1):
    data = np.asarray(data)
    upper = getBounds(data, zScore)

    if upper > data.max():
        upper = data.max()
        upper_outliers = False
    else:
        upper_outliers = True

    fig = plt.figure()

    color = "tab:blue"

    if bins == -1:
        bins = 'auto'
    elif bins == "single":
        bins = range(1, int(upper))
        plt.xticks(bins)
        color = "forestgreen"
    else:
        sys.error("Error: Unexpected bins argument: " + bins)

    if minimum != -1:
        n, bins, patches = plt.hist(data, range=(minimum, maximum), bins=bins,
                                    align='left')
    else:
        n, bins, patches = plt.hist(data, range=(data.min(), upper), bins=bins,
                                    align='left', color=color)

    if upper_outliers:
        n_upper_outliers = (data > upper).sum()
        patches[-1].set_height(patches[-1].get_height() + n_upper_outliers)
        patches[-1].set_facecolor('m')
        patches[-1].set_label('Range of upper outliers: (' + str(int(upper)) + ', ' + str(data.max()) + ')')

    if upper_outliers:
        plt.legend()

    ax = plt.gca()
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Count")
    report.savefig(fig)
    return data.min(), upper


def printIntronHistogram(prediction, report):
    data = prediction.getIntronLengths()
    histogram(data, report, "Histogram of intron lengths", "Intron length", 1)

def printIntergenicHistogram(prediction, report):
    data = prediction.getIntergenicLengths()
    histogram(data, report, "Histogram of intergenic lengths", "Intergenic length", 3)


def printExonHistogram(prediction, report):
    minimum, maximum = histogram(prediction.getExonLengths(), report,
                                 "Histogram of exon lengths", "Exon length", 3)
    histogram(prediction.getExonLengthsByType("Initial"), report,
              "Histogram of intial exon lengths", "Initial exon length", 3,
              minimum, maximum)
    histogram(prediction.getExonLengthsByType("Internal"), report,
              "Histogram of internal exon lengths", "Internal exon length", 3,
              minimum, maximum)
    histogram(prediction.getExonLengthsByType("Terminal"), report,
              "Histogram of terminal exon lengths", "Terminal exon length", 3,
              minimum, maximum)
    histogram(prediction.getExonLengthsByType("Single"), report,
              "Histogram of single exon lengths", "Single exon length", 3,
              minimum, maximum)


def printGeneHistogram(prediction, report):
    data = prediction.getGeneLengths()
    histogram(data, report, "Histogram of gene lengths", "Gene length", 3)


def printExonsPerGene(prediction, report):
    data = prediction.getExonsPerGene()
    histogram(data, report, "Exons per gene", "Exon number", 5,
              bins='single')


def printMotif(model, report, header, title):
    motif = motifs.logoFromMat(model, header, title)
    if motif:
        report.savefig(motif)


def printMotifs(model, report):
    printMotif(model, report, "$INI_MAT", "Start motif")
    printMotif(model, report, "$DONOR_0_MAT", "GT donor motif (phase 0)")
    printMotif(model, report, "$DONOR_1_MAT", "GT donor motif (phase 1)")
    printMotif(model, report, "$DONOR_2_MAT", "GT donor motif (phase 2)")
    printMotif(model, report, "$DONOR_GC_0_MAT", "GC donor motif (phase 0)")
    printMotif(model, report, "$DONOR_GC_1_MAT", "GC donor motif (phase 1)")
    printMotif(model, report, "$DONOR_GC_2_MAT", "GC donor motif (phase 2)")
    printMotif(model, report, "$BRANCH_MAT", "Branch point motif")
    dist = motifs.BPlengthDist(model)
    if dist:
        report.savefig(dist)
    printMotif(model, report, "$ACC_BP_0_MAT", "Short acceptor motif (phase 0)")
    printMotif(model, report, "$ACC_BP_1_MAT", "Short acceptor motif (phase 1)")
    printMotif(model, report, "$ACC_BP_2_MAT", "Short acceptor motif (phase 2)")
    printMotif(model, report, "$ACCEPTOR_0_MAT", "Acceptor motif (phase 0)")
    printMotif(model, report, "$ACCEPTOR_1_MAT", "Acceptor motif (phase 1)")
    printMotif(model, report, "$ACCEPTOR_2_MAT", "Acceptor motif (phase 2)")
    printMotif(model, report, "$TERM_TAG_MAT", "Amber motif")
    printMotif(model, report, "$TERM_TGA_MAT", "Umber motif")
    printMotif(model, report, "$TERM_TAA_MAT", "Ochre motif")



def main():
    args = parseCmd()
    report = PdfPages(args.output)

    prediction = analysis.PredictionAnalysis(args.genemark)
    printBasicStatistics(prediction, report)
    printGeneHistogram(prediction, report)
    printExonHistogram(prediction, report)
    printExonsPerGene(prediction, report)
    printIntronHistogram(prediction, report)
    printIntergenicHistogram(prediction, report)

    printMotifs(args.model, report)

    report.close()


def parseCmd():

    parser = argparse.ArgumentParser(description='Generate GeneMark \
        annotation report.')

    parser.add_argument('genemark', metavar='genemark.gtf', type=str,
                        help='GeneMark prediction file.')
    parser.add_argument('model', metavar='gmhmm.mod', type=str,
                        help='GeneMark model file (usually located in the \
                        GeneMark "output" folder).')
    parser.add_argument('output', metavar='output.pdf', type=str,
                        help='Name of the output pdf report.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
