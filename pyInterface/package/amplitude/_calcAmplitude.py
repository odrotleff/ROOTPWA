
import array

import pyRootPwa
import pyRootPwa.utils

def calcAmplitudes(inFile, keyfile, outfile):

	prodKinParticles = inFile.Get(pyRootPwa.config.prodKinPartNamesObjName)
	decayKinParticles = inFile.Get(pyRootPwa.config.decayKinPartNamesObjName)

	inTree = inFile.Get(pyRootPwa.config.inTreeName)

	prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
	decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")

	inTree.SetBranchAddress(pyRootPwa.config.prodKinMomentaLeafName, prodKinMomenta)
	inTree.SetBranchAddress(pyRootPwa.config.decayKinMomentaLeafName, decayKinMomenta)

	pythonAdmin = pyRootPwa.pythonAdministrator()
	if not pythonAdmin.constructAmplitude(keyfile):
		pyRootPwa.utils.printWarn('Could not construct amplitude for keyfile "' + keyfile + '".')
		return False
	if not pythonAdmin.initKinematicsData(prodKinParticles, decayKinParticles):
		pyRootPwa.utils.printErr('Could not initialize kinematics Data "' + keyfile + '".')
		return False

	progressbar = pyRootPwa.utils.progressBar(0, inTree.GetEntries())
	progressbar.start()
	for treeIndex in range(inTree.GetEntries()):
		inTree.GetEntry(treeIndex)
		if not pythonAdmin.readKinematicsData(prodKinMomenta, decayKinMomenta):
			progressbar.cancel()
			pyRootPwa.utils.printErr('Could not read kinematics data.')
			return False
		amp = pythonAdmin()
		if pyRootPwa.config.outputFileFormat == "ascii":
			outfile.write("(" + str(amp.real) + "," + str(amp.imag) + ")\n")
		elif pyRootPwa.config.outputFileFormat == "binary":
			arrayAmp = array.array('d', [amp.real, amp.imag])
			arrayAmp.tofile(outfile)
		elif pyRootPwa.config.outputFileFormat == "root":
			pyRootPwa.utils.printErr('Root output file format not implemented yet!')
			sys.exit(1)
		else:
			raise Exception('Something is wrong, this should have been checked in the initialization of the configuration!')
		progressbar.update(treeIndex)

	return True
