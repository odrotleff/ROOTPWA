#include"calcMultipleAmplitudes.h"
#include <boost/numeric/conversion/cast.hpp>
#include"calcTprime.h"


std::vector<double> rpwa::hli::getTprimesFromFiles(					const std::string			fileName,
											const std::string			keyFile,
											const long int				maxNmbEvents){

	waveDescription description = waveDescription();
	rpwa::isobarDecayTopologyPtr topology;
	description.parseKeyFile(keyFile);
	description.constructDecayTopology(topology);
	TFile* inFile = new TFile(fileName.c_str(),"READ");
	const eventMetadata* actData = eventMetadata::readEventFile(inFile);
	return getTprimesFromTree(actData, topology,maxNmbEvents);
};


std::vector<double> rpwa::hli::getTprimesFromTree(					const rpwa::eventMetadata* 		eventMeta,
											rpwa::isobarDecayTopologyPtr		topology,
											const long int				maxNmbEvents,
											const std::string& 			treePerfStatOutFileName,
											const long int 				treeCacheSize){

	if (not topology){
		printWarn<<"no topology given, abort."<<std::endl;
		return std::vector<double>();
	};
	if (not eventMeta){
		printWarn<<"problem with data file, abort"<<std::endl;	

	};
	TTree* tree = (*eventMeta).eventTree();
	if(not tree) {
		printErr << "event tree not found." << std::endl;
		return std::vector<double>();
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

	if(not topology->initKinematicsData((*eventMeta).productionKinematicsParticleNames(), (*eventMeta).decayKinematicsParticleNames())){
		printWarn << "problems initializing input data. cannot read input data." << std::endl;
		return std::vector<double>();
	};
	const long int    nmbEventsTree     = tree->GetEntries();
	const long int    nmbEvents         = ((maxNmbEvents > 0) ? std::min(maxNmbEvents, nmbEventsTree) : nmbEventsTree);
	std::vector<double> returnVector(nmbEvents);
	for(long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex){
		if(tree->LoadTree(eventIndex) < 0) {
			break;
		};

		prodKinMomentaBr->GetEntry(eventIndex);
		decayKinMomentaBr->GetEntry(eventIndex);

		if(not prodKinMomenta or not decayKinMomenta) {
			printWarn << "at least one of the input data arrays is a null pointer: "
			          << "        production kinematics: " << "momenta = " << prodKinMomenta  << std::endl
			          << "        decay kinematics:      " << "momenta = " << decayKinMomenta << std::endl
			          << "skipping event." << std::endl;
			return std::vector<double>();
		};
		if (not topology->readKinematicsData(*prodKinMomenta,*decayKinMomenta)){
			printWarn<< "could not set kinematics for event "<<eventIndex<<std::endl;
			return returnVector;
		};
		topology->calcIsobarLzVec();

		returnVector[eventIndex] = topology->getTprime();
		if (eventIndex%1000==0){
			printInfo<<eventIndex<<std::endl;
		};	
	};
	if(treePerfStats){
		treePerfStats->SaveAs(treePerfStatOutFileName.c_str());
		delete treePerfStats;
	};

	return returnVector;
};
