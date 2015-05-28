import os

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

def readWaveList(waveListFileName, keyFiles):
	pyRootPwa.utils.printInfo("reading amplitude names and thresholds from wave list file "
	          + "'" + waveListFileName + "'.")
	with open(waveListFileName, 'r') as waveListFile:
# 	if (not waveListFile) {
# 		printErr << "cannot open file '" << waveListFileName << "'. aborting." << endl;
# 		throw;
# 	}
		waveNames = []
		waveDescriptions = []
		waveThresholds = []
		lineNmb = 0
		for line in waveListFile:
			if (line[0] == '#'):  # comments start with #
				continue
			line = line.replace('\n', '')
			lineArray = line.split(" ")
			if(len(lineArray) >= 1 and len(lineArray) <= 2):
				waveName = lineArray[0]
				waveNames.append(waveName)
				if(len(lineArray) == 1):
					threshold = 0
				else:
					threshold = lineArray[1]
				waveDesc = pyRootPwa.core.waveDescription()
				waveDesc.parseKeyFile(keyFiles[waveName])
				waveDescriptions.append(waveDesc)
				waveThresholds.append(float(threshold))
			else:
				pyRootPwa.utils.printWarn("cannot parse line '" + line + "' in wave list file "
				          + "'" + waveListFileName + "'.")
#  			if (_debug):
#  				printDebug("reading line " + lineNmb + 1 + ": " + waveName + ", "
#  				           + "threshold = " + threshold + " MeV/c^2")
			lineNmb += 1
	pyRootPwa.utils.printInfo("read " + str(lineNmb) + " lines from wave list file " + "'" + waveListFileName + "'")
	return (waveNames, waveDescriptions, waveThresholds)


def pwaFit(ampFileList, normIntegralFileName, accIntegralFileName, binningMap, waveListFileName, keyFiles, seed=0, maxNmbEvents=0, startValFileName="", accEventsOverride=0, checkHessian=False, rank=1, verbose=False):
	treeDict = {}
	(waveNames, waveDescriptions, waveThresholds) = readWaveList(waveListFileName, keyFiles)
	massBinCenter = (binningMap['mass'][1] + binningMap['mass'][0]) / 2. # YOU CAN DO BETTER

	likelihood = pyRootPwa.core.pwaLikelihood()
	likelihood.useNormalizedAmps(True)
	if (not verbose):
		likelihood.setQuiet()
	if (not likelihood.init(
	                        rank,
	                        massBinCenter,
	                        waveDescriptions,
	                        waveThresholds)):
		printErr("could not initialize likelihood. Aborting...")
		sys.exit(1)

	normIntFile = ROOT.TFile.Open(normIntegralFileName, "READ")
	if len(normIntFile.GetListOfKeys()) != 1:
		pyRootPwa.utils.printWarn("'" + normIntegralFileName + "' does not contain exactly one TKey.")
		return False
	normIntMatrix = normIntFile.Get(pyRootPwa.core.ampIntegralMatrix.integralObjectName)
	likelihood.addNormIntegral(normIntMatrix)
	normIntFile.Close()
	accIntFile = ROOT.TFile.Open(accIntegralFileName, "READ")
	if len(accIntFile.GetListOfKeys()) != 1:
		pyRootPwa.utils.printWarn("'" + normIntegralFileName + "' does not contain exactly one TKey.")
		return False
	accIntMatrix = accIntFile.Get(pyRootPwa.core.ampIntegralMatrix.integralObjectName)
	likelihood.addAccIntegral(accIntMatrix)
	accIntFile.Close()

	for waveName in waveNames:
		ampFileName = ampFileList[waveName]
		ampFile = ROOT.TFile.Open(ampFileName, "READ")
		if not ampFile:
			pyRootPwa.utils.printErr("could not open amplitude file '" + ampFileName + "'.")
		meta = ampFile.Get(waveName + ".meta")
		if not meta:
			pyRootPwa.utils.printErr("could not get metadata for waveName '" + waveName + "'.")
			return False
		tree = ampFile.Get(waveName + ".amp")
		if not tree:
			pyRootPwa.utils.printErr("could not get amplitude tree for waveName '" + waveName + "'.")
			return False
		likelihood.addAmplitude(tree, meta)
	likelihood.rescaleIntegrals()
	lowerBound = binningMap[binningMap.keys()[0]][0]
	upperBound = binningMap[binningMap.keys()[0]][1]
	fitResult = pyRootPwa.core.pwaFit(
	                                  likelihood = likelihood,
	                                  seed = seed,
	                                  massBinMin = lowerBound,
	                                  massBinMax = upperBound,
	                                  checkHessian = checkHessian,
	                                  verbose = verbose
	                                  )
	return fitResult
	return False
