#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

def getDiffListAbs(r1, i1, r2, i2):
	diff = []
	for i in range(len(r1)):
		diff.append(abs(r1[i] - r2[i]) + abs(i1[i] - i2[i]))
	return diff

def getDiffListRel(r1, i1, r2, i2):
	diff = []
	for i in range(len(r1)):
		diff.append((abs(r1[i] - r2[i]) + abs(i1[i] - i2[i])) / (abs(r1[i] + r2[i]) + abs(i1[i] + i2[i])))
	return diff

def saveHistogram(path, name, title, lst):
	minHist = minimum(lst)
	maxHist = maximum(lst)
	if maxHist == 0: maxHist = 10**(-15)
	hist = ROOT.TH1D(name, title, 100, minHist, maxHist)
	for val in lst:
		hist.Fill(val)
	histFile = ROOT.TFile.Open(path + "/" + name + ".root", "RECREATE")
	hist.Write()
	histFile.Close()
	return hist

def minimum(lst):
	currMin = lst[0]
	for val in lst[1:]:
		if val < currMin: currMin = val
	return currMin

def maximum(lst):
	currMax = lst[0]
	for val in lst[1:]:
		if val > currMax: currMax = val
	return currMax

def extractRealImagListsFromAmpFile(fileName):
	ampFile = ROOT.TFile(fileName, "READ")
	tree = None
	for currKey in ampFile.GetListOfKeys():
		if currKey.GetName()[-3:] == "amp": tree = ampFile.Get(currKey.GetName())

	real = []
	imag = []
	for currEvent in tree:
		leaf  = currEvent.__getattr__(pyRootPwa.core.amplitudeMetadata.amplitudeLeafName)
		real.append(leaf.amp().real())
		imag.append(leaf.amp().imag())
	return (real, imag)

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="compare amplitude files"
	                                )

	parser.add_argument("file1", type=str, metavar="<filename>", help="path to first amplitude file")
	parser.add_argument("file2", type=str, metavar="<filename>", help="path to second amplitude file")
	args = parser.parse_args()

	pyRootPwa.utils.stdoutisatty = sys.stdout.isatty()
	pyRootPwa.utils.stderrisatty = sys.stderr.isatty()

	real1, imag1 = extractRealImagListsFromAmpFile(args.file1)
	real2, imag2 = extractRealImagListsFromAmpFile(args.file2)
	if not len(real1) == len(real2):
		pyRootPwa.utils.printErr("both trees have to have the same number of entries. Aborting...")
		sys.exit(1)

	diffAbs = getDiffListAbs(real1, imag1, real2, imag2)
	diffRel = getDiffListRel(real1, imag1, real2, imag2)
	diffAbsSum = 0
	diffRelAvg = 0
	for currentDiff in diffAbs:
		diffAbsSum += currentDiff
	for currentDiff in diffRel:
		diffRelAvg += currentDiff
	diffRelAvg /= len(diffRel)

	pyRootPwa.utils.printSucc("comparing the files '" + args.file1 + "' and '" + args.file2 + "' done.")
	pyRootPwa.utils.printSucc("sum of absolute delta: " + str(diffAbsSum))
	pyRootPwa.utils.printSucc("average of relative delta: " + str(diffRelAvg))

	saveHistogram("Test", "absolute delta", diffAbs)
