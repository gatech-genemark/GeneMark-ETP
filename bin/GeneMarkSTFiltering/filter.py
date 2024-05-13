#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2022, Georgia Institute of Technology, USA
#
# Filter GMST predictions
# ==============================================================


import argparse
import csv
import re
import os
import tempfile
import subprocess
import sys
import math
import shutil
from collections import OrderedDict
import logging
import removeConflictingPredictions

tempFiles = []
binDir = ""
diamondDir = ""
prothintDir = ""
outputFolder = ""


def extractFeatureGtf(text, feature):
    regex = feature + ' "([^"]+)"'
    return re.search(regex, text).groups()[0]


def getSignature(row):
    return f'{row[0]}_{row[3]}_{row[4]}'


def temp(prefix, suffix):
    tmp = tempfile.NamedTemporaryFile("w", delete=False,
                                      dir=outputFolder + "/tmp",
                                      prefix=prefix, suffix=suffix)
    tempFiles.append(tmp.name)
    return tmp


def cleanup():
    for file in tempFiles:
        os.remove(file)
    shutil.rmtree(outputFolder + "/tmp")


def systemCall(cmd):
    logging.info("Running the following system call: " + cmd)
    if subprocess.call(["bash", "-c", cmd]) != 0:
        logging.error('Program exited due to an error in command: ' + cmd)
        logging.error('Check stderr for more details.')
        sys.exit(1)


class Prediction():
    def __init__(self, ID, geneID, completeStatus, predClass, logOdds):
        self.ID = ID
        self.geneID = geneID
        self.stop = False
        self.complete = True
        self.predClass = predClass
        self.length = 0
        self.logOdds = float(logOdds)
        self.bitScore = -sys.maxsize
        self.AAIdentity = None
        if completeStatus == "partial":
            self.complete = False

    def addStop(self):
        self.stop = True

    def updateLength(self, row):
        self.length += int(row[4]) - int(row[3]) + 1

    def isPartialAt5End(self):
        return self.stop and self.predClass == "upLORF"


class Predictions(object):
    def __init__(self, gtfFile=None):
        if gtfFile:
            self.loadPredictions(gtfFile)
        self.diamondOut = None

    def loadPredictions(self, gtfFile):
        self.gtfFile = gtfFile
        self.predictions = {}
        for row in csv.reader(open(gtfFile), delimiter='\t'):
            if len(row) != 9:
                continue

            ID = extractFeatureGtf(row[8], "transcript_id")
            geneID = extractFeatureGtf(row[8], "gene_id")
            completeStatus = extractFeatureGtf(row[8], "status")
            predClass = extractFeatureGtf(row[8], "class")

            if ID not in self.predictions:
                self.predictions[ID] = Prediction(ID, geneID, completeStatus,
                                                  predClass, row[5])

            if row[2] == "stop_codon":
                self.predictions[ID].addStop()

            if row[2] == "CDS":
                self.predictions[ID].updateLength(row)

        return self.predictions

    def filterComplete(self, diamondDB, args):
        if not self.diamondOut:
            self.diamondOut = runDiamond(args.genome, self.gtfFile, diamondDB,
                                         args.topProteins, args.threads)

        notFound = set(self.predictions.keys())
        selected = set()
        for ID, query in self.diamondOut.items():
            # Save bitscore for later use
            self.predictions[ID].bitScore, \
                self.predictions[ID].AAIdentity = query.getBestScores()

            if query.checkFullProteinMatch(args):
                selected.add(ID)
                notFound.remove(ID)

        return self.selectSubsetByIds(selected), \
            self.selectSubsetByIds(notFound)

    def filterIncomplete(self, diamondDB, args, shorterPredictions):
        if not self.diamondOut:
            self.diamondOut = runDiamond(args.genome, self.gtfFile, diamondDB,
                                         args.topProteins, args.threads)

        selected = set()
        selectedShortened = set()
        for ID, query in self.diamondOut.items():
            # Because we are re-using the diamond result
            if ID not in self.predictions:
                continue

            # Save bitscore for later use. Save the score for shortened
            # selection as well in case it is selected
            bitScore, AAIdentity = query.getBestScores()
            self.predictions[ID].bitScore = bitScore
            self.predictions[ID].AAIdentity = AAIdentity
            shorterPredictions.predictions[ID].bitScore = bitScore
            shorterPredictions.predictions[ID].AAIdentity = AAIdentity

            # Whether to use the full partial prediction or shorten it to be
            # safe in case it is actually complete.
            if query.checkFullPartialMatch(args):
                selected.add(ID)
            else:
                if query.checkEndMatch(args):
                    selectedShortened.add(ID)

        return self.selectSubsetByIds(selected), \
            shorterPredictions.selectSubsetByIds(selectedShortened)

    def filterByLength(self, args):
        selected = set()
        for ID, pred in self.predictions.items():
            if pred.length >= args.minUnsupportedLength:
                selected.add(ID)
        return self.selectSubsetByIds(selected)

    def filterByLogOdds(self, args):
        selected = set()
        for ID, pred in self.predictions.items():
            if pred.logOdds >= args.minUnsupportedLogOdds:
                selected.add(ID)
        return self.selectSubsetByIds(selected)

    def selectSubsetByClass(self, classes):
        outputPreds = {}

        for ID, pred in self.predictions.items():
            if pred.predClass in classes:
                outputPreds[ID] = pred

        output = Predictions()
        output.setPredictions(outputPreds, self.gtfSubset(outputPreds))
        return output

    def selectSubsetByIds(self, ids):
        outputPreds = {}

        for ID, pred in self.predictions.items():
            if ID in ids:
                outputPreds[ID] = pred

        output = Predictions()
        output.setPredictions(outputPreds, self.gtfSubset(outputPreds))
        return output

    def selectSubsetWithStops(self):
        outputPreds = {}

        for ID, pred in self.predictions.items():
            if pred.stop:
                outputPreds[ID] = pred

        output = Predictions()
        output.setPredictions(outputPreds, self.gtfSubset(outputPreds))
        return output

    def gtfSubset(self, preds):
        outputGmstFile = temp("gtfFile", ".gtf")
        for row in csv.reader(open(self.gtfFile), delimiter='\t'):
            if len(row) != 9:
                continue
            ID = extractFeatureGtf(row[8], "transcript_id")
            if ID in preds:
                outputGmstFile.write("\t".join(row) + "\n")
        outputGmstFile.close()
        return outputGmstFile.name

    def setPredictions(self, preds, gtfFile):
        self.gtfFile = gtfFile
        self.predictions = preds

    def print(self, outName):
        shutil.copyfile(self.gtfFile, outName)

    def getIDs(self):
        return set(self.predictions.keys())

    def storeDiamondOut(self, diamondOut):
        self.diamondOut = diamondOut

    def filterIsoforms(self, threshold):
        bestBitScores = getBestBitScoresPerGene(self.predictions)

        selected = set()
        for ID, pred in self.predictions.items():
            if pred.bitScore >= bestBitScores[pred.geneID] * threshold:
                selected.add(ID)
        return self.selectSubsetByIds(selected)

    def getBestIsoforms(self):
        # Select isoforms with higher bitscore. If the bitscore of all isoforms
        # is not set, use length
        representatives = {}
        for pred in self.predictions.values():
            if pred.geneID not in representatives:
                representatives[pred.geneID] = pred
                continue

            if pred.bitScore > representatives[pred.geneID].bitScore:
                representatives[pred.geneID] = pred

            if representatives[pred.geneID].bitScore == -sys.maxsize:
                # None of the predictions so far had a bitscore set, use
                # length as a criterion instead
                if pred.length > representatives[pred.geneID].length:
                    representatives[pred.geneID] = pred

        selected = set()
        for pred in representatives.values():
            selected.add(pred.ID)
        return self.selectSubsetByIds(selected)

    def getLongestIsoforms(self):
        representatives = {}
        for pred in self.predictions.values():
            if pred.geneID not in representatives:
                representatives[pred.geneID] = pred
            elif pred.length > representatives[pred.geneID].length:
                representatives[pred.geneID] = pred

        selected = set()
        for pred in representatives.values():
            selected.add(pred.ID)
        return self.selectSubsetByIds(selected)

    def removeIfLongerIsoform(self, longestIsoforms):
        # Remove if there is another isoform of the same gene
        # (in longestIsoforms)
        selected = set()
        gene2length = {}
        for pred in longestIsoforms.predictions.values():
            gene2length[pred.geneID] = pred.length
        for pred in self.predictions.values():
            if pred.length >= gene2length[pred.geneID]:
                selected.add(pred.ID)
        return self.selectSubsetByIds(selected)


class Hit():
    def __init__(self, row):
        self.query = row[0]
        self.target = row[1]
        self.AAIdentity = float(row[2])
        self.length = int(row[3])
        self.qstart = int(row[7])
        self.qend = int(row[8])
        self.tstart = int(row[9])
        self.tend = int(row[10])
        self.bitScore = float(row[12])
        self.qlen = int(row[13])
        self.tlen = int(row[14])
        self.qLeft = self.qlen - self.qend + 1
        self.tLeft = self.tlen - self.tend + 1
        self.normalizedBitScore = self.bitScore / self.length


class Query():
    def __init__(self, row):
        self.ID = row[0]
        self.length = int(row[13])
        self.hits = OrderedDict()

    def addHit(self, row):
        target = row[1]
        self.hits[target] = Hit(row)

    def getBestScores(self):
        bitScore = -10000
        AAIdentity = -10000
        for hit in self.hits.values():
            bitScore = max(hit.bitScore, bitScore)
            AAIdentity = max(hit.AAIdentity, AAIdentity)
        return bitScore, AAIdentity

    def checkFullProteinMatch(self, args):
        for _, hit in self.hits.items():
            if hit.qstart <= args.maxqstart and \
                    abs(hit.qstart - hit.tstart) <= args.startdiff and \
                    hit.qLeft <= args.maxqend and \
                    abs(hit.qLeft - hit.tLeft) <= args.enddiff:
                return True
        return False

    def checkFullPartialMatch(self, args):
        # Only look at the first hit because this needs to be strict
        hit = self.hits[list(self.hits)[0]]
        if hit.qstart == 1 and hit.qLeft <= args.maxqend and \
           abs(hit.qLeft - hit.tLeft) <= args.enddiff:
            return True
        return False

    def checkEndMatch(self, args):
        for hit in self.hits.values():
            if hit.qLeft <= args.maxqend and \
               abs(hit.qLeft - hit.tLeft) <= args.enddiff:
                return True
        return False

    def checkIfPartial(self, shorterQuery, args):
        i = 1
        for hitID, hit in self.hits.items():
            if hitID not in shorterQuery.hits:
                if i == 1:
                    # If the best match does not even hit the shorter query
                    return True
                continue
            shorterHit = shorterQuery.hits[hitID]

            AAIdentityLogDiff = math.log(hit.AAIdentity /
                                         shorterHit.AAIdentity)
            weightedAADiff = args.partialAAIdentityWeight * AAIdentityLogDiff
            tStartDiff = shorterHit.tstart - hit.tstart

            # The partialness score characterizes how much more of the target
            # protein is covered in the alignment against the longer qurey.
            # This difference is compared against how many N-terminus AAs
            # are unaligned in the alignment of the longer query.
            # AA identity difference between the longer and shorter alignments
            # can further offset the score.
            partialnessScore = tStartDiff + 1 - hit.qstart + weightedAADiff

            if (hit.qstart == 1) or (
               partialnessScore > 0 and
               shorterHit.tstart != 1):
                return True
            i += 1
        return False


def selectSubsetRawGmst(gmst, subset):
    output = temp("gmstSubset", ".gtf")
    for row in csv.reader(open(gmst), delimiter='\t'):
        if len(row) != 9:
            continue
        if row[0] in subset:
            output.write("\t".join(row) + "\n")
    output.close()
    return output.name


def gmst2Hints(gmst, tseq, ggtf, complete=False, lorf=False, extend=False,
               subset=set()):

    if len(subset) != 0:
        gmst = selectSubsetRawGmst(gmst, subset)

    gtfOut = temp("gmst", ".gtf").name
    completeFlag = ''
    if complete:
        completeFlag = "--complete"
    lorfFlag = ''
    if lorf:
        lorfFlag = "--long"

    if lorf and complete:
        completeFlag = ""
        lorfFlag = "--LORFStart"

    extendFlag = ''
    if extend:
        extendFlag = "--upmax"

    systemCall(f'{binDir}/gms2hints.pl --tseq {tseq} --ggtf {ggtf} \
               --tgff {gmst} --out {gtfOut} \
               {completeFlag} {lorfFlag} {extendFlag}')

    return Predictions(gtfOut)


def makeDiamondDB(proteins, threads):
    db = temp("diamondDB", ".dmnd").name

    threadsArg = ""
    if threads:
        threadsArg = f"--threads {threads}"
    logging.info(f'Making diamond database from {proteins}')
    systemCall(f'{diamondDir}diamond makedb --in {proteins} -d {db} \
               {threadsArg}')
    return db


def runDiamond(genome, gmst, diamondDB, maxTargets, threads):
    query = temp("queryProteins", ".fasta").name
    output = temp("diamond", ".out").name
    systemCall(f'{binDir}/proteins_from_gtf.pl --seq {genome} --annot {gmst} \
               --out {query}')

    threadsArg = ""
    if threads:
        threadsArg = f"--threads {threads}"
    logging.info(f'Running DIAMOND with the query: {query}')
    systemCall(f'{diamondDir}diamond blastp --query {query} --db {diamondDB} \
               --outfmt 6 qseqid sseqid pident length positive mismatch \
               gapopen qstart qend sstart send evalue bitscore qlen slen \
               --max-target-seqs {maxTargets} --more-sensitive \
               {threadsArg} > {output}')

    return loadDiamondResult(output)


def loadDiamondResult(diamondOut):
    queries = {}
    for row in csv.reader(open(diamondOut), delimiter='\t'):
        if row[0] not in queries:
            queries[row[0]] = Query(row)
        queries[row[0]].addHit(row)
    return queries


def selectPartialCandidates(rawQueries, shorterQueries, args):
    partial = set()
    complete = set()
    for ID, shorterQuery in shorterQueries.items():
        if rawQueries[ID].checkIfPartial(shorterQuery, args):
            partial.add(ID)
        else:
            complete.add(ID)
    return partial, complete


def selectProteinSubset(alignedQueries, queryList, proteins):
    outputFile = temp("proteins", ".fasta")
    targets = set()
    for ID, query in alignedQueries.items():
        if ID in queryList:
            for hit in query.hits:
                targets.add(hit)

    selected = False
    p = open(proteins, "r")
    for line in p:
        if line[0] == ">":
            selected = False
            ID = line.split()[0][1:]
            if ID in targets:
                selected = True
        if selected:
            outputFile.write(line)

    outputFile.close()
    return outputFile.name


def classifyPartial(longerPredictions, shorterPredictions,
                    allProteinsDB, args):
    longerDiamond = runDiamond(args.genome, longerPredictions.gtfFile,
                               allProteinsDB, args.topProteins, args.threads)
    proteinsForShorter = selectProteinSubset(longerDiamond,
                                             shorterPredictions.getIDs(),
                                             args.proteins)
    shorterDiamond = runDiamond(args.genome, shorterPredictions.gtfFile,
                                makeDiamondDB(proteinsForShorter,
                                              args.threads), 0, args.threads)
    partialSet = set()
    completeSet = set()
    for ID, shorterQuery in shorterDiamond.items():
        if ID in longerDiamond and \
           longerDiamond[ID].checkIfPartial(shorterQuery, args):
            partialSet.add(ID)
        else:
            completeSet.add(ID)

    partial = longerPredictions.selectSubsetByIds(partialSet)
    complete = shorterPredictions.selectSubsetByIds(completeSet)
    partial.storeDiamondOut(longerDiamond)
    return partial, complete


def combinePredictions(predictionsList):
    combinedGtf = temp("combined", ".gtf")
    combinedPredictions = {}
    for predictions in predictionsList:
        combinedPredictions.update(predictions.predictions)
        with open(predictions.gtfFile, "r") as f:
            for line in f:
                combinedGtf.write(line)
    combinedGtf.close()
    result = Predictions()
    result.setPredictions(combinedPredictions, combinedGtf.name)
    return result


def getBestBitScoresPerGene(predictions):
    bestBitScores = {}
    for pred in predictions.values():
        if pred.geneID not in bestBitScores:
            bestBitScores[pred.geneID] = pred.bitScore
        else:
            bestBitScores[pred.geneID] = max(pred.bitScore,
                                             bestBitScores[pred.geneID])
    return bestBitScores


def selectRepresentatives(complete, incomplete):
    # Determine how to represent genes which have both complete and
    # incomplete transcripts. The selection is based on protein alignment score
    completeScores = getBestBitScoresPerGene(complete.predictions)
    incompleteScores = getBestBitScoresPerGene(incomplete.predictions)

    selectedComplete = set()
    selectedIncomplete = set()
    for ID, pred in complete.predictions.items():
        if pred.geneID not in incompleteScores:
            selectedComplete.add(ID)
        elif completeScores[pred.geneID] >= incompleteScores[pred.geneID]:
            selectedComplete.add(ID)

    for ID, pred in incomplete.predictions.items():
        if pred.geneID not in completeScores:
            selectedIncomplete.add(ID)
        elif completeScores[pred.geneID] < incompleteScores[pred.geneID]:
            selectedIncomplete.add(ID)

    return complete.selectSubsetByIds(selectedComplete), \
        incomplete.selectSubsetByIds(selectedIncomplete)


def selectUniqueGenes(existing, candidates):
    # Select unique genes in the "candidates" set which are not already
    # present in the "existing" set
    existingGenes = set()
    selectedCandidates = set()
    for pred in existing.predictions.values():
        existingGenes.add(pred.geneID)

    for ID, pred in candidates.predictions.items():
        if pred.geneID not in existingGenes:
            selectedCandidates.add(ID)

    return candidates.selectSubsetByIds(selectedCandidates)


def runProtHint(outDirName, predictions, args):
    # Run ProtHint on a set of predictions
    if not os.path.exists(outDirName):
        os.makedirs(outDirName)
    seeds = outDirName + "/seeds.gtf"
    predictions.print(seeds)

    threadsArg = ""
    if args.threads:
        threadsArg = f"--threads {args.threads}"
    systemCall(f'{prothintDir}prothint.py {args.genome} {args.proteins} \
               --geneSeeds {seeds} --workdir {outDirName} --longGene 30000000 \
               --longProtein 15000000 {threadsArg}')


def setup(args):
    global binDir
    global outputFolder
    global diamondDir
    global prothintDir
    outputFolder = args.outputFolder
    diamondDir = args.DIAMOND_PATH
    prothintDir = args.PROTHINT_PATH
    binDir = os.path.abspath(os.path.dirname(__file__))
    if not os.path.exists(args.outputFolder):
        os.makedirs(args.outputFolder)
    if not os.path.isdir(outputFolder + "/tmp"):
        os.mkdir(outputFolder + "/tmp")
    logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s',
                        datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO,
                        handlers=[
                            logging.FileHandler(outputFolder + "/log"),
                            logging.StreamHandler()
                        ])


def main():
    args = parseCmd()
    setup(args)
    logging.info("Starting the GMST filtering and classification.")

    rawPredictions = gmst2Hints(args.gmst, args.tseq, args.ggtf)
    # Save this for later filtering.
    longestIsoformsAll = rawPredictions.getLongestIsoforms()
    rawPredictions = rawPredictions.selectSubsetWithStops()
    allProteinsDB = makeDiamondDB(args.proteins, args.threads)

    logging.info("Preparing complete alternatives for incomplete predictions.")
    upLorfPredictions = rawPredictions.selectSubsetByClass(["upLORF"])
    shortenedPredictions = gmst2Hints(args.gmst, args.tseq, args.ggtf,
                                      subset=upLorfPredictions.getIDs(),
                                      complete=True)

    logging.info("Classifying incomplete predictions")
    upLorfP, upLorfC = classifyPartial(upLorfPredictions,
                                       shortenedPredictions,
                                       allProteinsDB, args)
    logging.info("Filtering incomplete predictions")
    upLorfPFull, upLorfPShort = upLorfP.filterIncomplete(allProteinsDB,
                                                         args,
                                                         shortenedPredictions)
    logging.info("Filtering shortened complete predictions")
    upLorfCFiltered, _ = upLorfC.filterComplete(allProteinsDB, args)

    logging.info("Filtering complete GMST predictions")
    comleteClasses = ["LORF_UPSTOP"]
    if not args.extend_LORF_NOUPSTOP:
        comleteClasses.append("LORF_NOUPSTOP")
    completePredictions = rawPredictions.selectSubsetByClass(comleteClasses)
    completeF, completeOut = completePredictions.filterComplete(allProteinsDB,
                                                                args)

    logging.info("Extending predictions with non 5'-most starts (when " +
                 "supported by proteins) to the 5'-most start.")
    sORFPredictions = rawPredictions.selectSubsetByClass(["sORF_NOUPSTOP",
                                                          "sORF_UPSTOP"])
    extendedSORFP = gmst2Hints(args.gmst, args.tseq, args.ggtf,
                               subset=sORFPredictions.getIDs(), lorf=True)
    extFiltered, notFound = extendedSORFP.filterComplete(allProteinsDB, args)
    sORFCandidates = sORFPredictions.selectSubsetByIds(notFound.getIDs())
    sORFFiltered, _ = sORFCandidates.filterComplete(allProteinsDB, args)

    logging.info("Filtering incomplete predictions with no start candidate")
    noORFPreds = rawPredictions.selectSubsetByClass(["noORF"])
    # Use the ones that are fully supported from qstart
    noORFPPFull, noORFPPShort = noORFPreds.filterIncomplete(allProteinsDB,
                                                            args,
                                                            noORFPreds)

    complete = combinePredictions([extFiltered, sORFFiltered, completeF,
                                  upLorfCFiltered])
    incomplete = combinePredictions([upLorfPFull, upLorfPShort, noORFPPFull])

    if args.extend_LORF_NOUPSTOP:
        logging.info("Classifying complete predictions with not upstream " +
                     "stop found in UTRs")
        completeNoUpS = rawPredictions.selectSubsetByClass(["LORF_NOUPSTOP"])
        extendedPredictions = gmst2Hints(args.gmst, args.tseq, args.ggtf,
                                         subset=completeNoUpS.getIDs(),
                                         extend=True)
        noUpStopP, noUpStopC = classifyPartial(extendedPredictions,
                                               completeNoUpS,
                                               allProteinsDB, args)
        noUpStFull, noUpStPShort = noUpStopP.filterIncomplete(allProteinsDB,
                                                              args,
                                                              completeNoUpS)
        noUpStopCFiltered, _ = noUpStopC.filterComplete(allProteinsDB, args)

        complete = combinePredictions([complete, noUpStopCFiltered])
        incomplete = combinePredictions([incomplete, noUpStFull,
                                         noUpStPShort])

    logging.info("Selecting representative isoforms.")
    incomplete1 = incomplete.getBestIsoforms()
    completeReps, incompleteReps = selectRepresentatives(complete,
                                                         incomplete1)

    incompleteReps.print(outputFolder + "/" + "incomplete.gtf")

    logging.info("Filtering complete alternative isoforms.")
    completeFiltered = completeReps.filterIsoforms(args.altIsoformsThreshold)

    logging.info("Adding reliable genes not fully supported by proteins.")
    unsupported = completeOut.selectSubsetByClass(["LORF_UPSTOP"])
    # Only select genes not represent in complete/incomplete groups
    unsupported = selectUniqueGenes(completeFiltered, unsupported)
    unsupported = selectUniqueGenes(incompleteReps, unsupported)
    unsupported = unsupported.filterByLength(args)
    # Isoform selection by bitscore was tested (for transcripts that do have
    # some protein alignment), but it led to worse results
    unsupported = unsupported.getLongestIsoforms()
    unsupported.print(outputFolder + "/unsupported_old.gtf")
    # Remove if there is another isoform of the same gene (not with the
    # LORF UPSTOP flag) which is longer. This often indicates an assem. error
    unsupported = unsupported.removeIfLongerIsoform(longestIsoformsAll)
    unsupported.print(outputFolder + "/unsupported.gtf")
    #unsupported = unsupported.filterByLogOdds(args)

    protHintOutFolder = outputFolder + "/prothint_unsupported"
    runProtHint(protHintOutFolder, unsupported, args)
    logging.info("Removing ProtHint-conflicting prediction in the" +
                 " unsupported gene set.")
    removeConflictingPredictions.remove(protHintOutFolder + "/seeds.gtf",
                                        protHintOutFolder + "/Spaln/spaln.gff",
                                        100,
                                        protHintOutFolder + "/no_confl.gtf")
    unsupNoConflicts = Predictions(protHintOutFolder + "/no_confl.gtf")
    unsupNoConflicts.print(outputFolder + "/unsupported_noConfl.gtf")
    unsupNoConflicts = unsupNoConflicts.filterByLogOdds(args)
    unsupNoConflicts.print(outputFolder + "/unsupported_noConfl_logodds.gtf")

    completeExt = combinePredictions([completeFiltered, unsupNoConflicts])
    completeExt.print(outputFolder + "/" + "complete.gtf")

    if args.printIntermediate:
        o = outputFolder + "/intermediate/"
        if not os.path.exists(o):
            os.makedirs(o)

        complete.print(o + "completeAll.gtf")
        incomplete.print(o + "incompleteAll.gtf")
        incomplete1.print(o + "incompleteAll_one_iso.gtf")
        completeReps.print(o + "completeReps.gtf")

        unsupported.print(o + "/unsupported.gtf")
        unsupNoConflicts.print(o + "/unsupported_no_conflicts.gtf")
        completeFiltered.print(o + "/" + "completeFullSupport.gtf")

        completePredictions.print(o + "LORF_Unfiltered.gtf")
        completeF.print(o + "LORF_Filtered.gtf")

        upLorfPredictions.print(o + "upLORF_All.gtf")
        shortenedPredictions.print(o + "upLORF_Shortened_All.gtf")
        upLorfP.print(o + "upLORF_Partial_Unfiltered.gtf")
        upLorfPFull.print(o + "upLORF_Partial_Full_Filtered.gtf")
        upLorfPShort.print(o + "upLORF_Partial_Shortened_Filtered.gtf")
        upLorfC.print(o + "upLORF_Complete_Unfiltered.gtf")
        upLorfCFiltered.print(o + "upLORF_Complete_Filtered.gtf")

        sORFPredictions.print(o + "sORF_All.gtf")
        extendedSORFP.print(o + "LORF_Start_All.gtf")
        extFiltered.print(o + "LORF_Start_Filtered.gtf")
        sORFCandidates.print(o + "sORF_Candidates.gtf")
        sORFFiltered.print(o + "sORF_Filtered.gtf")

        noORFPreds.print(o + "noORF.gtf")
        noORFPPFull.print(o + "noORF_Filtered_Full.gtf")
        noORFPPShort.print(o + "noORF_Filtered_Short.gtf")

        if args.extend_LORF_NOUPSTOP:
            completeNoUpS.print(o + "NoUpStop_All.gtf")
            extendedPredictions.print(o + "NoUpStop_Extended_All.gtf")
            noUpStopP.print(o + "NoUpStop_Partial_Unfiltered.gtf")
            noUpStFull.print(o + "noUpStop_Partial_Full_Filtered.gtf")
            noUpStPShort.print(o + "noUpStop_Partial_Shortened_Filtered.gtf")
            noUpStopC.print(o + "noUpStop_Complete_Unfiltered.gtf")
            noUpStopCFiltered.print(o + "noUpStop_Complete_Filtered.gtf")

    logging.info("Deleting temporary files.")
    cleanup()
    if not args.printIntermediate:
        shutil.rmtree(protHintOutFolder)
    logging.info("GMST filtering and classification finished.")


def parseCmd():

    parser = argparse.ArgumentParser(description='Filter GMST predictions')

    parser.add_argument('gmst', metavar="gmst.gtf", type=str, help='GMST \
        predictions on transcript sequences')
    parser.add_argument('tseq', metavar='transcripts.fasta', type=str,
                        help='Transcript nucleotide sequences.')
    parser.add_argument('ggtf', metavar='transcripts.gff', type=str,
                        help='Transcript coordinates.')
    parser.add_argument('proteins', metavar='proteins.fasta', type=str,
                        help='Protein database')
    parser.add_argument('genome', metavar='genome.fasta', type=str,
                        help='Protein database')

    parser.add_argument('--outputFolder', default="out", type=str,
                        help='Output folder name')

    parser.add_argument('--maxqstart', type=int, default=999999,
                        help='Description')

    parser.add_argument('--startdiff', type=int, default=5,
                        help='Description')

    parser.add_argument('--maxqend', type=int, default=999999,
                        help='Description')

    parser.add_argument('--enddiff', type=int, default=20,
                        help='Description')

    parser.add_argument('--topProteins', type=int, default=25,
                        help='How many target proteins to use during\
                        filtering')

    parser.add_argument('--partialAAIdentityWeight', type=float, default=1000)

    parser.add_argument('--minUnsupportedLength', type=int, default=300,
                        help='Length threshold for HC genes not supported by\
                        proteins. Default = 300')

    parser.add_argument('--minUnsupportedLogOdds', type=float, default=50,
                        help='Log-odds score threshold for HC genes not \
                        supported by proteins. Default = 50')

    parser.add_argument('--altIsoformsThreshold', type=float,
                        help='All reported complete HC isoforms of a gene\
                        must have an alignment bitScore highter than \
                        ALTISOFORMSTHRESHOLD * best bitscore of any \
                        complete isoform in the gene.',
                        default=0.8)

    parser.add_argument('--printIntermediate', action='store_true',
                        help='Print intermediate filtering results. Useful \
                        for debugging and collectin statistics.')

    parser.add_argument('--extend_LORF_NOUPSTOP', action='store_true',
                        help='Try to extend complete predictions with no\
                        upstream stop in the UTR region. If the extended\
                        prediction has better support than the default one,\
                        classify it as incomplete. This option is turned\
                        off by default because GMST is usually correct in\
                        prediction complete genes (even when the upstream\
                        stop codon is not present).')

    parser.add_argument('--DIAMOND_PATH', type=str, default='',
                        help='Path to folder with DIAMOND binary. If not\
                        specified, the program will look for DIAMOND in the\
                        $PATH.')

    parser.add_argument('--PROTHINT_PATH', type=str, default='',
                        help='Path to folder with ProtHint binary. If not\
                        specified, the program will look for ProtHint in the\
                        $PATH.')

    parser.add_argument('--threads', type=int,
                        help='Number of threads used by DIAMOND and ProtHint\
                        By default, all available threads are used.')

    args = parser.parse_args()

    if args.DIAMOND_PATH != '':
        args.DIAMOND_PATH = os.path.abspath(args.DIAMOND_PATH) + '/'

    if args.PROTHINT_PATH != '':
        args.PROTHINT_PATH = os.path.abspath(args.PROTHINT_PATH) + '/'

    return args


if __name__ == '__main__':
    main()
