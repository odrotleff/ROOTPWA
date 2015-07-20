#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.core

# tPrime Binning
# 0.100000-0.112853
# 0.112853-0.127471
# 0.127471-0.144385
# 0.144385-0.164401
# 0.164401-0.188816
# 0.188816-0.219907
# 0.219907-0.262177
# 0.262177-0.326380
# 0.326380-0.448588
# 0.448588-0.724294
# 0.724294-1.000000

def readWaveList(waveListFileName):
	pyRootPwa.utils.printInfo("reading amplitude names and thresholds from wave list file "
	                        + "'" + waveListFileName + "'.")
	try:
		with open(waveListFileName, 'r') as waveListFile:
			waveNamesFromWavelist = []
			waveThresholds = []
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
					waveNamesFromWavelist.append(waveName)
					waveThresholds.append(float(threshold))
				else:
					pyRootPwa.utils.printWarn("cannot parse line '" + line + "' in wave list file "
					                        + "'" + waveListFileName + "'.")
	#  			if (_debug):
	#  				printDebug("reading line " + lineNmb + 1 + ": " + waveName + ", "
	#  				           + "threshold = " + threshold + " MeV/c^2")
				lineNmb += 1
	except IOError:
		pyRootPwa.utils.printError("cannot open file '" + waveListFileName + "'. Aborting...")
		sys.exit(1)
	pyRootPwa.utils.printInfo("read " + str(lineNmb) + " lines from wave list file " + "'" + waveListFileName + "'")
	return (waveNamesFromWavelist, waveThresholds)

if __name__ == "__main__":

	pyRootPwa.utils.stdoutisatty = sys.stdout.isatty()
	pyRootPwa.utils.stderrisatty = sys.stderr.isatty()

	# parse command line arguments
	parser = argparse.ArgumentParser(
	                                 description="calculates decay amplitudes "
	                                             "for given wave for events in "
	                                             "input data files and "
	                                             "writes amplitudes to file"
	                                )

	parser.add_argument("outputFileNameBase", type=str, metavar="fileName", help="base to output file")
	parser.add_argument("-c", type=str, metavar="configFileName", default="rootpwa.config", dest="configFileName", help="path to config file (default: ./rootpwa.config)")
	parser.add_argument("-n", type=int, metavar="#", default=-1, dest="maxNmbEvents",  help="maximum number of events to read (default: all)")
	parser.add_argument("-b", type=int, metavar="massBin", default=-1, dest="massBin", help="mass bin to be calculated (default: all)")
	parser.add_argument("-B", type=str, metavar="tPrimeBinning", default="", dest="tPrimeBinning", help="tPrime binning (default: none)")
	parser.add_argument("-e", type=str, metavar="eventsType", default="all", dest="eventsType", help="events type to be calculated ('generated' or 'accepted', default: all)")
	parser.add_argument("-f", "--no-progress-bar", action="store_true", dest="noProgressBar", help="disable progress bars (decreases computing time)")
	parser.add_argument("-k", "--keyfiles", type=str, metavar="keyfiles", dest="keyfiles", nargs="*", help="keyfiles to calculate amplitude for (overrides settings from the config file)")
	parser.add_argument("-w", type=str, metavar="wavelistFileName", default="", dest="wavelistFileName", help="path to wavelist file (default: none)")
	parser.add_argument("-v", action="store_true", dest="debug", help="verbose; print debug output (default: false)")
	args = parser.parse_args()

	pyRootPwa.ROOT.gROOT.ProcessLine("#include <complex>")
	config = pyRootPwa.rootPwaConfig()
	if not config.initialize(args.configFileName):
		pyRootPwa.utils.printErr("loading config file '" + args.configFileName + "' failed. Aborting...")
		sys.exit(1)
	pyRootPwa.core.particleDataTable.readFile(config.pdgFileName)
	fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
	if not fileManager:
		pyRootPwa.utils.printErr("loading the file manager failed. Aborting...")
		sys.exit(1)

	waveList = []
	if (not args.wavelistFileName==""):
		(waveList, waveThresholds) = readWaveList(args.wavelistFileName)
	if (len(waveList) == 0):
		waveList = fileManager.getWaveNameList()

	binIDList = fileManager.getBinIDList()
	if (not args.massBin == -1):
		binIDList = [args.massBin]

	eventsTypes = []
	if (args.eventsType == "generated"):
		eventsTypes = [ pyRootPwa.core.eventMetadata.GENERATED ]
	elif (args.eventsType == "accepted"):
		eventsTypes = [ pyRootPwa.core.eventMetadata.ACCEPTED ]
	elif (args.eventsType == "all"):
		eventsTypes = [ pyRootPwa.core.eventMetadata.GENERATED,
		                pyRootPwa.core.eventMetadata.ACCEPTED ]
	else:
		pyRootPwa.utils.printErr("Invalid events type given ('" + args.eventsType + "'). Aborting...")
		sys.exit(1)

	tPrimeBinning = [0.1,0.112853,0.127471,0.144385,0.164401,0.188816,0.219907,0.262177,0.32638,0.448588,0.724294,1.]
	if not args.tPrimeBinning == "":
		tPrimeBinning = [float(bin) for bin in args.tPrimeBinning.split(";")]

	for binID in binIDList:
		for eventsType in eventsTypes:
			dataFile = fileManager.getDataFile(binID, eventsType)
			if not dataFile:
				continue
			print dataFile.dataFileName
			if not pyRootPwa.core.getTbinnedIntegralsFromKeyFiles(args.outputFileNameBase,
			                                                      fileManager.getKeyFilePaths(),
			                                                      [dataFile.dataFileName],
			                                                      tPrimeBinning,
			                                                      args.maxNmbEvents):
				pyRootPwa.utils.printWarn("could not calculate integrals.")
