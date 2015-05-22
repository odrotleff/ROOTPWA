#include"calcMultipleAmplitudes.h"
#include <boost/numeric/conversion/cast.hpp>
#include"calcAmplitude.h"


bool rpwa::hli::getIntegralsFromKeyFiles(			const std::string			integralName,
								const std::string			outFileName,
								const std::vector<std::string> 		&keyFiles, 
								const std::vector<std::string> 		&eventFiles, 
								const int 				maxNmbEvents){

	if (keyFiles.size()==0){
		printWarn<<"no keyFiles given, abort"<<std::endl;
		return false;
	};

	std::vector<ampIntegralMatrix> integralMatrix = std::vector<ampIntegralMatrix>(1,ampIntegralMatrix());
	std::vector<rpwa::isobarAmplitudePtr> amplitudes = getAmplitudesFromKeyFiles(keyFiles);
	std::vector<std::string> names = waveNamesFromKeyFiles(keyFiles, true);
	std::vector<double> tBinning(0);
	integralMatrix[0].initialize(names);
	if (eventFiles.size()==0){
		printWarn<<"no eventFiles given, abort"<<std::endl;
		return false;
	};
	for (size_t nFile = 0;nFile<eventFiles.size();++nFile){
		TFile* inFile = new TFile(eventFiles[nFile].c_str(),"READ");
		const eventMetadata* actData = eventMetadata::readEventFile(inFile);
		int nEvents = integralMatrix[0].nmbEvents();
		if (nEvents >= maxNmbEvents and maxNmbEvents!=-1){
			break;
		};
		if (not calcTbinnedIntegralsFromEventTree(actData, amplitudes, integralMatrix, std::vector<double>(0), maxNmbEvents-nEvents, 0)){
			printWarn<< "error in the integration, abort"<<std::endl;
			return false;
		};
	};
	TFile* outFile = new TFile(outFileName.c_str(), "RECREATE");
	integralMatrix[0].Write(integralName.c_str());
	outFile->Close();
	return true;
};



std::vector<std::complex<double> > rpwa::hli::evaluateAmplitudes(	std::vector<rpwa::isobarAmplitudePtr> 	&amplitudes, 
									TClonesArray	 			prodKinematics, 
									TClonesArray 				decayKinematics){
	std::vector<std::complex<double> > returnValue(amplitudes.size());
	for (size_t amp=0;amp<amplitudes.size();++amp){
		if (not amplitudes[amp]->decayTopology()->readKinematicsData(prodKinematics,decayKinematics)){
			printErr<<"could not set kinematics for";
			amplitudes[amp]->decayTopology()->print(printErr);
			return std::vector<std::complex<double> >();
		};
		returnValue[amp] = amplitudes[amp]->amplitude();
	};
	return returnValue;
};

bool rpwa::hli::initAmplitudesKinematics(	std::vector<rpwa::isobarAmplitudePtr> 	&amplitudes, 
						std::vector<std::string> 		prodNames, 
						std::vector<std::string> 		decayNames){
	bool returnValue = true;
	for(size_t amp=0;amp<amplitudes.size();++amp){
		if(not amplitudes[amp]->decayTopology()->initKinematicsData(prodNames,decayNames)){
			printErr<<"could not initialize kinematics for ";
			amplitudes[amp]->decayTopology()->print(printErr);
			returnValue = false;
		};
	};
	return returnValue;
};

std::vector<rpwa::isobarAmplitudePtr> rpwa::hli::getAmplitudesFromKeyFiles(const std::vector<std::string> &keyFiles){
	size_t nWaves = keyFiles.size();
	std::vector<rpwa::isobarAmplitudePtr>returnValue(nWaves,rpwa::isobarAmplitudePtr());
	waveDescription description = waveDescription();
	rpwa::isobarDecayTopologyPtr topology;
	for (size_t wave=0;wave<nWaves;++wave){
		description.parseKeyFile(keyFiles[wave]);
		description.constructAmplitude(returnValue[wave]);
		returnValue[wave]->init();
		description.constructDecayTopology(topology);
	};
	return returnValue;
};

std::vector<std::string> rpwa::hli::waveNamesFromKeyFiles(	const std::vector<std::string> 	&keyFiles, 
								bool 				newConvention){
	size_t nWaves = keyFiles.size();
	std::vector<std::string> names(nWaves);
	waveDescription description = waveDescription();
	rpwa::isobarDecayTopologyPtr topology;
	for (size_t wave=0;wave<nWaves;++wave){
		description.parseKeyFile(keyFiles[wave]);	
		description.constructDecayTopology(topology);
		names[wave] = description.waveNameFromTopology(*topology,newConvention);
	};
	return names;
};

template<typename T>
void printVector(std::vector<T> vec){
	
	std::cout<<"[";
	for (size_t i=0;i<vec.size();++i){
		std::cout<<vec[i];
	};
	std::cout<<"]"<<std::endl;
};

bool rpwa::hli::calcTbinnedIntegralsFromEventTree(	const rpwa::eventMetadata* 			eventMeta, 
							std::vector<rpwa::isobarAmplitudePtr> 		&amplitudes, 
							std::vector<rpwa::ampIntegralMatrix>		&matrix,
							std::vector<double>				tBinning,
							const long int 					maxNmbEvents,
							const long int					startEvent,
							const std::string& 				treePerfStatOutFileName,
							const long int 					treeCacheSize){
	

	bool tBinned = true;
	if (tBinning.size() ==0){
		tBinned = false;
	}else if(tBinning.size() != matrix.size()+1){
		printErr<<"number of tbins does not correpond to the number of matrices"<<std::endl;
		return false;
	};
	size_t nWaves = amplitudes.size();
	if (nWaves ==0){
		printWarn<<"no amplitudes given, abort."<<std::endl;
		return false;
	};
	for (size_t amp=0;amp<nWaves;++amp){
		if (not amplitudes[amp]){
			printWarn<<"null pointer to isobar decay amplitude #"<<amp<<". cannot process tree."<<std::endl;
			return false;
		};
	};
	TTree* tree = (*eventMeta).eventTree();
	if(not tree) {
		printErr << "event tree not found." << std::endl;
		return false;
	};

	// create branch pointers and leaf variables
	TBranch*      prodKinMomentaBr  = 0;
	TBranch*      decayKinMomentaBr = 0;
	TClonesArray* prodKinMomenta    = 0;
	TClonesArray* decayKinMomenta   = 0;

	// connect leaf variables to tree branches
	tree->SetBranchAddress(eventMetadata::productionKinematicsMomentaBranchName.c_str(),  &prodKinMomenta,  &prodKinMomentaBr );
	tree->SetBranchAddress(eventMetadata::decayKinematicsMomentaBranchName.c_str(), &decayKinMomenta, &decayKinMomentaBr);

	tree->SetCacheSize(treeCacheSize);
	tree->AddBranchToCache(eventMetadata::productionKinematicsMomentaBranchName.c_str(),  true);
	tree->AddBranchToCache(eventMetadata::decayKinematicsMomentaBranchName.c_str(), true);
	tree->StopCacheLearningPhase();

	TTreePerfStats* treePerfStats = 0;
	if(treePerfStatOutFileName != "") {
		treePerfStats = new TTreePerfStats("ioPerf", tree);
	};

	// loop over events
	if(not initAmplitudesKinematics(amplitudes,(*eventMeta).productionKinematicsParticleNames(), (*eventMeta).decayKinematicsParticleNames())) {
		printWarn << "problems initializing input data. cannot read input data." << std::endl;
		return false;
	};
	const long int    nmbEventsTree     = tree->GetEntries();
	const long int    nmbEvents         = ((maxNmbEvents > 0) ? std::min(maxNmbEvents, nmbEventsTree) : nmbEventsTree);
	for(long int eventIndex = startEvent; eventIndex < nmbEvents; ++eventIndex){
		if(tree->LoadTree(eventIndex) < 0) {
			break;
		};

		prodKinMomentaBr->GetEntry(eventIndex);
		decayKinMomentaBr->GetEntry(eventIndex);

		size_t tBin=0;

		if (tBinned){

		};
			

		if(not prodKinMomenta or not decayKinMomenta) {
			printWarn << "at least one of the input data arrays is a null pointer: "
			          << "        production kinematics: " << "momenta = " << prodKinMomenta  << std::endl
			          << "        decay kinematics:      " << "momenta = " << decayKinMomenta << std::endl
			          << "skipping event." << std::endl;
			return false;
		};
		std::vector<std::complex<double> > ampl = evaluateAmplitudes(amplitudes,*prodKinMomenta,*decayKinMomenta);
		matrix[tBin].addEventAmplitudes(ampl); // Here chose the right t' bin
		if (eventIndex%1000==0){
			printInfo<<eventIndex<<std::endl;
		};	
	};
	if(treePerfStats){
		treePerfStats->SaveAs(treePerfStatOutFileName.c_str());
		delete treePerfStats;
	};

	return true;
};
