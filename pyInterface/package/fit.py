import os

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

def pwaFit(ampFileList, normIntegralFileName, accIntegralFileName, maxNmbEvents=0, weightFileName=""):
	treeDict = {}
	waveNames = []
	ampFiles = []

	normIntFile = ROOT.TFile.Open(normIntegralFileName, "READ")
	if len(normIntFile.GetListOfKeys()) != 1:
		pyRootPwa.utils.printWarn("'" + normIntegralFileName + "' does not contain exactly one TKey.")
		return False
	normIntMatrix = normIntFile.Get(normIntFile.GetListOfKeys()[0].GetName())
	accIntFile = ROOT.TFile.Open(normIntegralFileName, "READ")
	if len(accIntFile.GetListOfKeys()) != 1:
		pyRootPwa.utils.printWarn("'" + normIntegralFileName + "' does not contain exactly one TKey.")
		return False
	accIntMatrix = accIntFile.Get(normIntFile.GetListOfKeys()[0].GetName())

	for ampFileName in ampFileList:
		ampFile = ROOT.TFile.Open(ampFileName, "READ")
		ampFiles.append(ampFile)
		if not ampFile:
			pyRootPwa.utils.printErr("could not open amplitude file '" + ampFileName + "'.")
		foundAmpKey = False
		for key in ampFile.GetListOfKeys():
			if not key:
				pyRootPwa.utils.printWarn("NULL pointer to TKey in file '" + ampFileName + "'.")
				continue
			keyName = key.GetName()
			keyWithoutExt, keyExt = os.path.splitext(keyName)
			if keyExt == ".amp":
				foundAmpKey = True
				tree = ampFile.Get(keyName)
				meta = ampFile.Get(keyWithoutExt + ".meta")
				if not meta:
					pyRootPwa.utils.printErr("could not get metadata for waveName '" + keyWithoutExt + "'.")
					del ampFiles
					return False
				waveNames.append(meta.objectBaseName())
				treeDict[meta.objectBaseName()] = tree
		if not foundAmpKey:
			pyRootPwa.utils.printWarn("no TKey in file '" + ampFileName + "'.")
	fitResult = pyRootPwa.core.pwaFit(treeDict, normIntMatrix, accIntMatrix)
	del ampFiles
	return fitResult
