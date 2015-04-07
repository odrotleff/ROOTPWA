#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.core

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="generate phase space Monte Carlo events"
	                                )

	parser.add_argument("outputFileName", type=str, metavar="fileName", help="path to output file")
	parser.add_argument("-c", type=str, metavar="configFileName", dest="configFileName", default="./rootpwa.config", help="path to config file (default: './rootpwa.config')")
	parser.add_argument("-b", type=int, metavar="#", dest="binID", default=0, help="bin ID of fit (default: 0)")
	parser.add_argument("-n", type=int, metavar="#", dest="nEvents", default=0, help="maximum number of events to process (default: all)")
	parser.add_argument("-s", type=int, metavar="#", dest="seed", default=0, help="random seed (default: 0)")
# 	parser.add_argument("-r", type=int, metavar="#", dest="nEventsRenorm", default=0, help="number of events to renormalize to (default: no renormalization)")
	parser.add_argument("-w", type=str, metavar="path", dest="weightsFileName", default="", help="path to MC weight file for de-weighting (default: none)")
	parser.add_argument("-v", "--verbose", help="verbose; print debug output (default: false)", action="store_true")
	args = parser.parse_args()

	pyRootPwa.utils.stdoutisatty = sys.stdout.isatty()
	pyRootPwa.utils.stderrisatty = sys.stderr.isatty()

	printErr  = pyRootPwa.utils.printErr
	printWarn = pyRootPwa.utils.printWarn
	printSucc = pyRootPwa.utils.printSucc
	printInfo = pyRootPwa.utils.printInfo
	printDebug = pyRootPwa.utils.printDebug

	config = pyRootPwa.rootPwaConfig()
	if not config.initialize(args.configFileName):
		pyRootPwa.utils.printErr("loading config file '" + args.configFileName + "' failed. Aborting...")
		sys.exit(1)
	pyRootPwa.core.particleDataTable.readFile(config.pdgFileName)
	fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
	if not fileManager:
		pyRootPwa.utils.printErr("loading the file manager failed. Aborting...")
		sys.exit(1)

	ampFileList = fileManager.getAmpFilePaths(args.binID, pyRootPwa.core.eventMetadata.REAL)
	if not ampFileList:
		printErr("could not retrieve valid amplitude file list. Aborting...")
		sys.exit(1)

	psIntegralPath  = fileManager.getIntegralFilePath(args.binID, pyRootPwa.core.eventMetadata.GENERATED)
	accIntegralPath = fileManager.getIntegralFilePath(args.binID, pyRootPwa.core.eventMetadata.ACCEPTED)
	fitResult = pyRootPwa.pwaFit(ampFileList, psIntegralPath, accIntegralPath, args.nEvents, args.weightsFileName)

	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFileName, "RECREATE")
	if not outputFile:
		printErr("cannot open output file '" + args.outputFileName + "'. Aborting...")
		sys.exit(1)
	nmbBytes = fitResult.Write(args.integralName)
	outputFile.Close()
	if nmbBytes == 0:
		printErr("problems writing integral to TKey '" + args.integralName + "' "
		         + "in file '" + args.outputFileName + "'")
		sys.exit(1)
	else:
		printSucc("wrote integral to TKey '" + args.integralName + "' "
		          + "in file '" + args.outputFileName + "'")
