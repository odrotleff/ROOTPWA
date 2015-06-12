#!/usr/bin/env python

import sys
import argparse
import numpy as np
import pyRootPwa
import pyRootPwa.core
from nuts import nuts6

def readWaveList(waveListFileName, keyFiles):
	pyRootPwa.utils.printInfo("reading amplitude names and thresholds from wave list file "
	          + "'" + waveListFileName + "'.")
	with open(waveListFileName, 'r') as waveListFile:
#	if (not waveListFile) {
#		printErr << "cannot open file '" << waveListFileName << "'. Aborting..." << endl;
#		throw;
#	}
		waveDescThres = []
		lineNmb = 0
		for line in waveListFile:
			if (line[0] == '#'):  # comments start with #
				continue
			line = line.replace('\n', '')
			lineArray = line.split(" ")
			if(len(lineArray) >= 1 and len(lineArray) <= 2):
				waveName = lineArray[0]
				if(len(lineArray) == 1):
					threshold = 0
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

def getStartValues(likelihood, seed):
	nmbPars = likelihood.nmbPars()
	random = pyRootPwa.ROOT.TRandom3(seed)
	startValues = np.zeros(nmbPars)
	sqrtNmbEvts = np.sqrt(likelihood.nmbEvents())
	for i in xrange(nmbPars):
		parName = likelihood.parName(i)
		# check if parameter needs to be fixed
		if (not likelihood.parFixed(i)):
			startVal = random.Uniform(defaultStartValue, sqrtNmbEvts);
			if(random.Rndm() > 0.5):
				startVal *= -1.;
			pyRootPwa.utils.printInfo("    setting parameter [" + str(i) + "] " + parName + " = %10f" % startVal)
			startValues[i] = startVal
		else:
			pyRootPwa.utils.printErr("Fixed variables not supported!")
			sys.exit(1)
	return startValues

def initLikelihood(fileManager, massBinCenter, binID, waveListFileName, accEventsOverride, rank, cauchyPriors, verbose):
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
	accIntegralFileName = fileManager.getIntegralFilePath(binID, pyRootPwa.core.eventMetadata.ACCEPTED)
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
		meta = pyRootPwa.core.amplitudeMetadata.readAmplitudeFile(ampFile, waveName)
		if not meta:
			pyRootPwa.utils.printErr("could not get metadata for waveName '" + waveName + "'.")
			return False
		if (not likelihood.addAmplitude(meta)):
			pyRootPwa.utils.printErr("could not add amplitude '" + waveName + "'. Aborting...")
			return False
	if (not likelihood.finishInit()):
		pyRootPwa.utils.printErr("could not finish initialization of likelihood. Aborting...")
		sys.exit(1)
	return likelihood


def FdFNoCauchy(par):
	likeli, grad = likelihood.FdF(par.tolist()) # convert numpy array to list
	grad = [-elem for elem in grad]
	return -likeli, np.asarray(grad) # convert list to numpy array

def FdFWithCauchy(par):
	likeli, grad = likelihoodWithCauchy.FdF(par.tolist()) # convert numpy array to list
	grad = [-elem for elem in grad]
	return -likeli, np.asarray(grad) # convert list to numpy array

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="pwa fit executable"
	                                )

	parser.add_argument("outputFileName", type=str, metavar="fileName", help="path to output file")
	parser.add_argument("-c", type=str, metavar="configFileName", dest="configFileName", default="./rootpwa.config", help="path to config file (default: './rootpwa.config')")
	parser.add_argument("-b", type=int, metavar="#", dest="binID", default=0, help="bin ID of fit (default: 0)")
	parser.add_argument("-s", type=int, metavar="#", dest="seed", default=0, help="random seed (default: 0)")
	parser.add_argument("-w", type=str, metavar="path", dest="waveListFileName", default="", help="path to wavelist file (default: none)")
	parser.add_argument("-N", type=int, metavar="#", dest="nmbSamples", default=10000, help="number of samples (default: 10,000)")
	parser.add_argument("-r", type=int, metavar="#", dest="rank", default=1, help="rank of spin density matrix (default: 1)")
	parser.add_argument("-A", type=int, metavar="#", dest="accEventsOverride", default=0, help="number of input events to normalize acceptance to (default: use number of events from acceptance integral file)")
 	parser.add_argument("-C", "--cauchyPriors", help="use half-Cauchy priors (default: false)", action="store_true")
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

	defaultStartValue = 0.01
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
	binningMap = fileManager.getBinFromID(args.binID)
	massBinCenter = (binningMap['mass'][1] + binningMap['mass'][0]) / 2.
	likelihood = initLikelihood(fileManager, massBinCenter, args.binID, args.waveListFileName, args.accEventsOverride, args.rank, args.cauchyPriors, args.verbose)

	print "setting start values"
	startValuesNoCauchy = getStartValues(likelihood, args.seed)

	M = args.nmbSamples
	Madapt = 5000
	delta = 0.2
	print('Running HMC without cauchy priors with dual averaging and trajectory length %0.2f...' % delta)
	samples, lnprob, epsilon = nuts6(FdFNoCauchy, M, Madapt, startValuesNoCauchy, delta)
	print('Done. Final epsilon = %f.' % epsilon)
# 	print('Running HMC with cauchy priors with dual averaging and trajectory length %0.2f...' % delta)
# 	samplesWithCauchy, lnprobWithCauchy, epsilonWithCauchy = nuts6(FdFWithCauchy, M, Madapt, startValuesWithCauchy, delta)
# 	print('Done. Final epsilonWithCauchy = %f.' % epsilonWithCauchy)

	if len(samples) != len(lnprob):
		print("dieded")
		raise Exception("asdf")
	tups = [ (lnprob[i], samples[i]) for i in xrange(len(lnprob))]
	tups = sorted(tups, key=lambda x: x[0], reverse=True)
 	tups = tups[:int(len(tups)*0.68)]
 	print("cut = " + str(tups[-1][0]))
# 	print('Mean')
# 	print (np.mean(samples, axis=0))
# 	print('Stddev')
# 	print (np.std(samples, axis=0))

	nmbPars = likelihood.nmbPars()
	valTreeName = "samples"
	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFileName, "NEW")
	if ((not outputFile) or outputFile.IsZombie()):
		pyRootPwa.utils.printErr("cannot open output file '" + args.outputFileName + "'. Aborting...")
		sys.exit(1)
	tree = pyRootPwa.ROOT.TTree(valTreeName, valTreeName)
	arr = []
	for par in xrange(nmbPars):	
		arr.append(np.zeros(1, dtype='d'))
	for par in xrange(nmbPars):
		tree.Branch("b_" + str(par),arr[par], "b_" + str(par) + "/D")
 	probval = np.zeros(1,dtype='d')
 	tree.Branch("lnprob",probval,"prob/D")
	for j in xrange(len(samples)):
		sample = samples[j]
		for i in xrange(nmbPars):
			arr[i][0] = sample[i]
 		probval[0] = lnprob[j]
		tree.Fill()
	tree.Write()

	(prodAmps, parIndices, prodAmpNames) = likelihood.buildProdAmpArrays([0. for _ in range(nmbPars)])
	fitParCovMatrix                      = pyRootPwa.ROOT.TMatrixD(0, nmbPars, 0, nmbPars)
# 	fitParCovIndices                     = [(likelihood.parName(p),  for p in range(nmbPars)]
	print prodAmps
	print prodAmpNames
	print "likelihood:" + str(tups[0][0])
# 	print fitParCovIndices
	fitResult = pyRootPwa.core.fitResult()
	fitResult.fill(likelihood.nmbEvents(),
	               1,
	               massBinCenter,
	               tups[0][0],
	               args.rank,
	               prodAmps,
	               prodAmpNames,
	               fitParCovMatrix,
	               parIndices,
	               pyRootPwa.core.complexMatrix(0,0),
	               pyRootPwa.core.complexMatrix(0,0),
	               [],
	               True,
	               False)
	valTreeName   = "pwa"
	valBranchName = "fitResult_v2"
	printInfo("file '" + args.outputFileName + "' is empty. "
	        + "creating new tree '" + valTreeName + "' for PWA result.")
	tree = pyRootPwa.ROOT.TTree(valTreeName, valTreeName)
	if not fitResult.branch(tree, valBranchName):
		printErr("failed to create new branch '" + valBranchName + "' in file '" + args.outputFileName + "'.")
		sys.exit(1)
	tree.Fill()
	nmbBytes = tree.Write()
	outputFile.Close()
