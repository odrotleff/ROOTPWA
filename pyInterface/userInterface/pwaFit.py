#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.core


if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="pwa fit executable"
	                                )

	parser.add_argument("outputFileName", type=str, metavar="fileName", help="path to output file")
	parser.add_argument("-c", type=str, metavar="configFileName", dest="configFileName", default="./rootpwa.config", help="path to config file (default: './rootpwa.config')")
	parser.add_argument("-b", type=int, metavar="#", dest="binID", default=0, help="bin ID of fit (default: 0)")
	parser.add_argument("-s", type=int, metavar="#", dest="seed", default=0, help="random seed (default: 0)")
	parser.add_argument("-w", type=str, metavar="path", dest="waveListFileName", default="", help="path to wavelist file (default: none)")
	parser.add_argument("-S", type=str, metavar="path", dest="startValFileName", default="", help="path to start value fit result file (default: none)")
	parser.add_argument("-r", type=int, metavar="#", dest="rank", default=1, help="rank of spin density matrix (default: 1)")
	parser.add_argument("-A", type=int, metavar="#", dest="accEventsOverride", default=0, help="number of input events to normalize acceptance to (default: use number of events from acceptance integral file)")
	parser.add_argument("-H", "--checkHessian", help="check analytical Hessian eigenvalues (default: false)", action="store_true")
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
		printErr("loading config file '" + args.configFileName + "' failed. Aborting...")
		sys.exit(1)
	pyRootPwa.core.particleDataTable.readFile(config.pdgFileName)
	fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
	if not fileManager:
		printErr("loading the file manager failed. Aborting...")
		sys.exit(1)

	ampFileList = fileManager.getAmpFilePaths(args.binID, pyRootPwa.core.eventMetadata.REAL)
	if not ampFileList:
		printErr("could not retrieve valid amplitude file list. Aborting...")
		sys.exit(1)
	binningMap = fileManager.getBinFromID(args.binID)
	massBinCenter = (binningMap['mass'][1] + binningMap['mass'][0]) / 2.

	psIntegralPath  = fileManager.getIntegralFilePath(args.binID, pyRootPwa.core.eventMetadata.GENERATED)
	accIntegralPath = fileManager.getIntegralFilePath(args.binID, pyRootPwa.core.eventMetadata.ACCEPTED)

	likelihood = pyRootPwa.core.pwaLikelihood()
	if (not args.verbose):
		likelihood.setQuiet(True)
	if (not likelihood.init(
	                        args.rank,
	                        ampFileList,
	                        massBinCenter,
	                        args.waveListFileName,
	                        psIntegralPath,
	                        accIntegralPath,
	                        args.accEventsOverride)):
		printErr("could not initialize likelihood. Aborting...")
		sys.exit(1)
	fitResult = pyRootPwa.pwaFit(
	                             likelihood,
	                             seed = args.seed,
	                             massBinLower = binningMap['mass'][0],
	                             massBinUpper = binningMap['mass'][1],
	                             startValFileName = args.startValFileName,
	                             checkHessian = args.checkHessian,
	                             verbose = args.verbose
	                             )
	printInfo("writing result to '" + args.outputFileName + "'")
	valTreeName   = "pwa"
	valBranchName = "fitResult_v2"
	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFileName, "UPDATE")
	if ((not outputFile) or outputFile.IsZombie()):
		printErr("cannot open output file '" + args.outputFileName + "'. Aborting...")
		sys.exit(1)
	tree = outputFile.Get(valTreeName)
	if (not tree):
		printInfo("file '" + args.outputFileName + "' is empty. "
		        + "creating new tree '" + valTreeName + "' for PWA result.")
		tree = pyRootPwa.ROOT.TTree(valTreeName, valTreeName)
		if not fitResult.branch(tree, valBranchName):
			printErr("failed to create new branch '" + valBranchName + "' in file '" + args.outputFileName + "'.")
			sys.exit(1)
	else:
		fitResult.setBranchAddress(tree, valBranchName)
	tree.Fill()
	nmbBytes = tree.Write()
	outputFile.Close()
	if nmbBytes == 0:
		printErr("problems writing integral to TKey 'fitResult' "
		       + "in file '" + args.outputFileName + "'")
		sys.exit(1)
	else:
		printSucc("wrote integral to TKey 'fitResult' "
		        + "in file '" + args.outputFileName + "'")
