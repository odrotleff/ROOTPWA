#!/usr/bin/env python

import sys
import argparse
import pyRootPwa
import pyRootPwa.core

def getFitResultFromFile(fitResultFileName):
	fitResultFile = ROOT.TFile.Open(fitResultFileName, 'READ')
	if not fitResultFile:
		print("Could not open generated fit result file '" + fitResultFileName + "'.")
		return None
	fitResultTree = fitResultFile.Get('pwa')
	result = pyRootPwa.core.fitResult()
	result.setBranchAddress(fitResultTree, 'fitResult_v2')
	if fitResultTree.GetEntries() != 1:
		print('More than one fit result in TTree, somebody should probably implement this properly...')
		fitResultFile.Close()
		return None
	fitResultTree.GetEntry(0)
	if not result.converged():
		print("Fit not converged for fit result file '" + fitResultFileName + "'.")
		fitResultFile.Close()
		return None
	return result

def readWaveList(waveListFileName, keyFiles):
	pyRootPwa.utils.printInfo("reading amplitude names and thresholds from wave list file "
	          + "'" + waveListFileName + "'.")
	waveDescThres = []
	if waveListFileName == '':
		for waveName in keyFiles:
			waveDesc = pyRootPwa.core.waveDescription()
			waveDesc.parseKeyFile(keyFiles[waveName])
			waveDescThres.append( (waveName, waveDesc, 0.) )
	else:
		with open(waveListFileName, 'r') as waveListFile:
			lineNmb = 0
			for line in waveListFile:
				if (line[0] == '#'):  # comments start with #
					continue
				line = line.replace('\n', '')
				lineArray = line.split(" ")
				if(len(lineArray) >= 1 and len(lineArray) <= 2):
					waveName = lineArray[0]
					if(len(lineArray) == 1):
						threshold = "0"
					else:
						threshold = lineArray[1]
					waveDesc = pyRootPwa.core.waveDescription()
					waveDesc.parseKeyFile(keyFiles[waveName])
					waveDescThres.append( (waveName, waveDesc, float(threshold)) )
				else:
					pyRootPwa.utils.printWarn("cannot parse line '" + line + "' in wave list file "
					          + "'" + waveListFileName + "'.")
	#  			if (_debug):
	#  				printDebug("reading line " + lineNmb + 1 + ": " + waveName + ", "
	#  				           + "threshold = " + threshold + " MeV/c^2")
				lineNmb += 1
			pyRootPwa.utils.printInfo("read " + str(lineNmb) + " lines from wave list file " + "'" + waveListFileName + "'")
	return waveDescThres

def initLikelihood(fileManager, massBinCenter, binID, waveListFileName, accEventsOverride, rank, cauchyPriors, verbose, genIntFilename="", accIntFilename="", addBinningMap={}, evtFileName=""):
	waveDescThres = readWaveList(waveListFileName, fileManager.getKeyFiles())
	likelihood = pyRootPwa.core.pwaLikelihood()
	likelihood.useNormalizedAmps(True)
	if cauchyPriors:
		likelihood.setPriorType(pyRootPwa.core.HALF_CAUCHY)
	if not verbose:
		likelihood.setQuiet()
	if (not likelihood.init(waveDescThres,
	                        rank,
	                        massBinCenter)):
		printErr("could not initialize likelihood. Aborting...")
		return False

	normIntegralFileName  = fileManager.getIntegralFilePath(binID, pyRootPwa.core.eventMetadata.GENERATED)
	if not genIntFilename == "":
		normIntegralFileName = genIntFilename
	accIntegralFileName = fileManager.getIntegralFilePath(binID, pyRootPwa.core.eventMetadata.ACCEPTED)
	if not accIntFilename == "":
		accIntegralFileName = accIntFilename
	normIntFile = ROOT.TFile.Open(normIntegralFileName, "READ")
	if len(normIntFile.GetListOfKeys()) != 1:
		pyRootPwa.utils.printWarn("'" + normIntegralFileName + "' does not contain exactly one TKey.")
		sys.exit(1)
	normIntMatrix = normIntFile.Get(pyRootPwa.core.ampIntegralMatrix.integralObjectName)
	if (not likelihood.addNormIntegral(normIntMatrix)):
		pyRootPwa.utils.printErr("could not add normalization integral. Aborting...")
		sys.exit(1)
	normIntFile.Close()
	accIntFile = ROOT.TFile.Open(accIntegralFileName, "READ")
	if len(accIntFile.GetListOfKeys()) != 1:
		pyRootPwa.utils.printWarn("'" + normIntegralFileName + "' does not contain exactly one TKey.")
		sys.exit(1)
	accIntMatrix = accIntFile.Get(pyRootPwa.core.ampIntegralMatrix.integralObjectName)
	if (not likelihood.addAccIntegral(accIntMatrix, accEventsOverride)):
		pyRootPwa.utils.printErr("could not add acceptance integral. Aborting...")
		sys.exit(1)
	accIntFile.Close()

	for wave in waveDescThres:
		waveName = wave[0]
		ampFileName = ampFileList[waveName]
		ampFile = ROOT.TFile.Open(ampFileName, "READ")
		if not ampFile:
			pyRootPwa.utils.printErr("could not open amplitude file '" + ampFileName + "'.")
			return False
		ampMeta = pyRootPwa.core.amplitudeMetadata.readAmplitudeFile(ampFile, waveName)
		if not ampMeta:
			pyRootPwa.utils.printErr("could not get metadata for waveName '" + waveName + "'.")
			return False
		evtFile = ROOT.TFile.Open(evtFileName, "READ")
		if not evtFile:
			pyRootPwa.utils.printErr("could not open amplitude file '" + evtFileName + "'.")
			return False
		evtMeta = pyRootPwa.core.eventMetadata.readEventFile(evtFile)
		if not evtMeta:
			pyRootPwa.utils.printErr("could not get metadata for event file '" + evtFileName + "'.")
			return False
		if (not likelihood.addAmplitude(ampMeta, addBinningMap, evtMeta)):
			pyRootPwa.utils.printErr("could not add amplitude '" + waveName + "'. Aborting...")
			return False
	if (not likelihood.finishInit()):
		pyRootPwa.utils.printErr("could not finish initialization of likelihood. Aborting...")
		sys.exit(1)
	return likelihood

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="pwa fit executable"
	                                )

	parser.add_argument("outputFileName", type=str, metavar="fileName", help="path to output file")
	parser.add_argument("-c", type=str, metavar="configFileName", dest="configFileName", default="./rootpwa.config", help="path to config file (default: './rootpwa.config')")
	parser.add_argument("-f", type=str, dest="fitResultPath", default="", help="path to fit result")
	parser.add_argument("-b", type=int, metavar="#", dest="binID", default=0, help="bin ID of fit (default: 0)")
	parser.add_argument("-B", dest="addBin", action='append', help="additional binning in the form 'binningVariable;lowerBound;upperBound' (e.g. 'mass;1000;1100')."+
	                                                             "You can use the argument multiple times for multiple binning variables")
	parser.add_argument("-s", type=int, metavar="#", dest="seed", default=0, help="random seed (default: 0)")
	parser.add_argument("-w", type=str, metavar="path", dest="waveListFileName", default="", help="path to wavelist file (default: none)")
	parser.add_argument("-A", type=int, metavar="#", dest="accEventsOverride", default=0, help="number of input events to normalize acceptance to (default: use number of events from acceptance integral file)")
 	parser.add_argument("-C", "--cauchyPriors", help="use half-Cauchy priors (default: false)", action="store_true")
	parser.add_argument("-g", type=str, metavar="integralPath", dest="genIntFilename", default="", help="phase space integral file override")
	parser.add_argument("-a", type=str, metavar="integralPath", dest="accIntFilename", default="", help="acceptance integral file override")
	parser.add_argument("-v", "--verbose", help="verbose; print debug output (default: false)", action="store_true")
	args = parser.parse_args()

	pyRootPwa.utils.stdoutisatty = sys.stdout.isatty()
	pyRootPwa.utils.stderrisatty = sys.stderr.isatty()

	printErr  = pyRootPwa.utils.printErr
	printWarn = pyRootPwa.utils.printWarn
	printSucc = pyRootPwa.utils.printSucc
	printInfo = pyRootPwa.utils.printInfo
	printDebug = pyRootPwa.utils.printDebug
	ROOT = pyRootPwa.ROOT

	config = pyRootPwa.rootPwaConfig()
	if not config.initialize(args.configFileName):
		printErr("loading config file '" + args.configFileName + "' failed. Aborting...")
		sys.exit(1)
	pyRootPwa.core.particleDataTable.readFile(config.pdgFileName)
	fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
	if not fileManager:
		printErr("loading the file manager failed. Aborting...")
		sys.exit(1)

	ampFileList = fileManager.getAmplitudeFilePaths(args.binID, pyRootPwa.core.eventMetadata.REAL)
	if not ampFileList:
		printErr("could not retrieve valid amplitude file list. Aborting...")
		sys.exit(1)
	if args.fitResultPath == "":
		printErr("no fit result given. Aborting...")
		sys.exit(1)
	result = getFitResultFromFile(args.fitResultPath)
	binningMap = fileManager.getBinFromID(args.binID)
	addBinningMap = pyRootPwa.utils.binningMapFromArgList(args.addBin)
	if not addBinningMap:
		printWarn("received no valid additional binning map argument")
	eventFile = fileManager.getDataFile(args.binID, pyRootPwa.core.eventMetadata.REAL)
	eventFileName = eventFile.dataFileName
	likelihood = initLikelihood(fileManager, result.massBinCenter(), args.binID, args.waveListFileName, args.accEventsOverride, result.rank(), args.cauchyPriors, args.verbose, args.genIntFilename, args.accIntFilename, addBinningMap, eventFileName)

	pars = []

	for i in range(likelihood.nmbPars()):
		parName = likelihood.parName(i);
		pars.append(result.fitParameter(parName))
	covMatrix = likelihood.CovarianceMatrix(pars)
	if args.verbose:
		covMatrix.Print()

	newResult = pyRootPwa.core.fitResult()
	newResult.fill(result.nmbEvents(),
	               result.normNmbEvents(),
	               result.massBinCenter(),
	               result.logLikelihood(),
	               result.rank(),
	               result.prodAmps(),
	               result.prodAmpNames(),
	               covMatrix,
	               result.fitParCovIndices(),
	               result.normIntegralMatrix(),
	               result.acceptedNormIntegralMatrix(),
	               result.phaseSpaceIntegralVector(),
	               result.converged(),
	               True)
	valTreeName   = "pwa"
	valBranchName = "fitResult_v2"
	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFileName, "NEW")
	if (not outputFile) or outputFile.IsZombie():
		pyRootPwa.utils.printErr("cannot open output file '" + args.outputFileName + "'. Aborting...")
		sys.exit(1)
	printInfo("file '" + args.outputFileName + "' is empty. "
	        + "creating new tree '" + valTreeName + "' for PWA result.")
	tree = pyRootPwa.ROOT.TTree(valTreeName, valTreeName)
	if not newResult.branch(tree, valBranchName):
		printErr("failed to create new branch '" + valBranchName + "' in file '" + args.outputFileName + "'.")
		sys.exit(1)
	tree.Fill()
	nmbBytes = tree.Write()
	outputFile.Close()
