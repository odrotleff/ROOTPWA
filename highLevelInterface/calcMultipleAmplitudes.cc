#include"calcMultipleAmplitudes.h"
#include <boost/numeric/conversion/cast.hpp>
#include"calcAmplitude.h"
#include<sstream>


bool rpwa::hli::getIntegralsFromKeyFiles(                       const std::string                       integralName,
                                                                const std::string                       outFileName,
                                                                const std::vector<std::string>          &keyFiles, 
                                                                const std::vector<std::string>          &eventFiles, 
                                                                const int                               maxNmbEvents){

	if (keyFiles.size()==0){
		printWarn<<"no keyFiles given, abort"<<std::endl;
		return false;
	};

	std::vector<ampIntegralMatrix> integralMatrix = std::vector<ampIntegralMatrix>(1,ampIntegralMatrix());
	std::vector<rpwa::isobarAmplitudePtr> amplitudes = getAmplitudesFromKeyFiles(keyFiles);
	std::vector<std::string> names = waveNamesFromKeyFiles(keyFiles);
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
		if (not calcBinnedIntegralsFromEventTree(actData, amplitudes, integralMatrix, std::vector<std::map<std::string,std::pair<double,double> > >(0), maxNmbEvents-nEvents, 0)){
			printWarn<< "error in the integration, abort"<<std::endl;
			return false;
		};
	};
	TFile* outFile = new TFile(outFileName.c_str(), "RECREATE");
	integralMatrix[0].Write(integralName.c_str());
	outFile->Close();
	return true;
};

bool rpwa::hli::getTbinnedIntegralsFromKeyFiles(                const std::string                       integralNameBase,
                                                                const std::string                       outFileName,
                                                                const std::vector<std::string>          &keyFiles,
                                                                const std::vector<std::string>          &eventFiles,
                                                                const std::vector<double>               &tBinning,
                                                                const int                               maxNmbEvents){

	if (tBinning.size() <2){
		printWarn<<"no t' bins given"<<std::endl;
		return false;
	};
	std::vector<std::map<std::string,std::pair<double,double> > > binning;
	for (size_t t=1;t<tBinning.size();++t){
		std::map<std::string,std::pair<double,double> > bin;
		bin["tPrime"] = std::pair<double,double>(tBinning[t-1],tBinning[t]);
		binning.push_back(bin);
	};

	if (keyFiles.size()==0){
		printWarn<<"no keyFiles given, abort"<<std::endl;
		return false;
	};

	size_t nTbin = binning.size();

	std::vector<ampIntegralMatrix> integralMatrix = std::vector<ampIntegralMatrix>(nTbin,ampIntegralMatrix());
	std::vector<rpwa::isobarAmplitudePtr> amplitudes = getAmplitudesFromKeyFiles(keyFiles);
	std::vector<std::string> names = waveNamesFromKeyFiles(keyFiles);
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
		if (not calcBinnedIntegralsFromEventTree(actData, amplitudes, integralMatrix, binning , maxNmbEvents-nEvents, 0)){
			printWarn<< "error in the integration, abort"<<std::endl;
			return false;
		};
	};
	TFile* outFile = new TFile(outFileName.c_str(), "RECREATE");
	for (size_t b=0;b<nTbin;++b){
		std::stringstream integralName;
		integralName<<integralNameBase<<"_"<<tBinning[b]<<"_"<<tBinning[b+1];
		integralMatrix[b].Write(integralName.str().c_str());
	};
	outFile->Close();
	integralMatrix[0].writeAscii("./integrals_new_method.txt");
	return true;
};   



std::vector<std::complex<double> > rpwa::hli::evaluateAmplitudes(       std::vector<rpwa::isobarAmplitudePtr>   &amplitudes, 
                                                                        TClonesArray                            prodKinematics, 
                                                                        TClonesArray                            decayKinematics){
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

bool rpwa::hli::initAmplitudesKinematics(       std::vector<rpwa::isobarAmplitudePtr>   &amplitudes, 
                                                std::vector<std::string>                prodNames, 
                                                std::vector<std::string>                decayNames){
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

std::vector<std::string> rpwa::hli::waveNamesFromKeyFiles(      const std::vector<std::string>  &keyFiles){
	size_t nWaves = keyFiles.size();
	std::vector<std::string> names(nWaves);
	waveDescription description = waveDescription();
	rpwa::isobarDecayTopologyPtr topology;
	for (size_t wave=0;wave<nWaves;++wave){
		description.parseKeyFile(keyFiles[wave]);
		description.constructDecayTopology(topology);
		names[wave] = description.waveNameFromTopology(*topology);
	};
	return names;
};

bool rpwa::hli::checkBinnings(std::vector<std::map<std::string,std::pair<double, double> > > binnings){

	size_t nBinnings = binnings.size();
	if (nBinnings == 0){
		printWarn<<"no binnings given"<<std::endl;
		return true; // Reurn true nevetheless, since this does not corrupt the calculation
	};
	std::vector<std::string> variables;
	typedef std::map<std::string, std::pair<double,double> >::iterator it_type;
	for (it_type iterator = binnings[0].begin(); iterator != binnings[0].end(); iterator++){
		variables.push_back(iterator->first);
	};
	size_t size = variables.size();
	for (size_t b=1;b<nBinnings;++b){
		if (not binnings[b].size() == size){
			printErr<<"different number of variables given"<<std::endl;
			return false;
		};
		for (size_t v=0;v<size;++v){
			if (binnings[b].count(variables[v]) == 0){
				printErr<<"variable "<<variables[v]<<" not in binning #"<<b<<std::endl;
				return false;
			};
		};
	};
	return true;
};

bool rpwa::hli::calcBinnedIntegralsFromEventTree(       const rpwa::eventMetadata*                      eventMeta, 
                                                        std::vector<rpwa::isobarAmplitudePtr>           &amplitudes, 
                                                        std::vector<rpwa::ampIntegralMatrix>            &matrix,
                                                        std::vector<std::map<std::string,std::pair<double,double> > >binning,
                                                        const long int                                  maxNmbEvents,
                                                        const long int                                  startEvent,
                                                        const std::string&                              treePerfStatOutFileName,
                                                        const long int                                  treeCacheSize){

	bool binned = true;
	if (binning.size() == 0){
		binned = false;
	}else if(binning.size() != matrix.size()){
		printErr<<"number of tbins does not correpond to the number of matrices"<<std::endl;
		return false;
	};
	size_t nWaves = amplitudes.size();
	size_t nMatrix= matrix.size();
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

	std::vector<std::string> variables(0);
	if (binned){
		if (not checkBinnings(binning)){
			printErr<<"problem with the binnings"<<std::endl;
			return false;
		};
		typedef std::map<std::string, std::pair<double,double> >::iterator it_type;
		for (it_type iterator = binning[0].begin(); iterator != binning[0].end(); iterator++){
			variables.push_back(iterator->first);
		};
	};
	size_t nVar = variables.size();

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

	std::vector<double> binVariables(nVar,0.);
	std::vector<TBranch*> variablesBranches(nVar,0);
	for (size_t v=0;v<nVar;++v){
		variablesBranches[v] = tree->GetBranch(variables[v].c_str());
		if (not variablesBranches[v]){
			printErr<<"branch "<<variables[v]<<" does not exist. abort integration"<<std::endl;
			return false;
		};
		variablesBranches[v]->SetAddress(&binVariables[v]);
	};

	// connect leaf variables to tree branches
	tree->SetBranchAddress(eventMetadata::productionKinematicsMomentaBranchName.c_str(),  &prodKinMomenta,  &prodKinMomentaBr );
	tree->SetBranchAddress(eventMetadata::decayKinematicsMomentaBranchName.c_str(), &decayKinMomenta, &decayKinMomentaBr);

	tree->SetCacheSize(treeCacheSize);
	tree->AddBranchToCache(eventMetadata::productionKinematicsMomentaBranchName.c_str(),  true);
	tree->AddBranchToCache(eventMetadata::decayKinematicsMomentaBranchName.c_str(), true);
	for(size_t v=0;v<nVar;++v){
		tree->AddBranchToCache(variables[v].c_str());
	};

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
		if(not prodKinMomenta or not decayKinMomenta) {
			printWarn << "at least one of the input data arrays is a null pointer: "
			          << "        production kinematics: " << "momenta = " << prodKinMomenta  << std::endl
			          << "        decay kinematics:      " << "momenta = " << decayKinMomenta << std::endl
			          << "abort integration." << std::endl;
			return false;
		};

		for (size_t v=0;v<nVar;++v){
			variablesBranches[v]->GetEntry(eventIndex);
		};		
		bool inAtAll = false;
		std::vector<bool> inBin(nMatrix,false);
		for (size_t b=0;b<nMatrix;++b){
			bool isIn = true;
			for (size_t v=0;v<nVar;++v){
				if (binning[b][variables[v]].first > binVariables[v] or binning[b][variables[v]].second < binVariables[v]){
					isIn = false;
					break;
				};
			};
			if (isIn){
				inBin[b] = true;
				if (not inAtAll){
					inAtAll = true;
				};
			};
		};
		if (binned and not inAtAll){
			continue; // Do not calculate amplitudes if not needed
		};

		std::vector<std::complex<double> > ampl = evaluateAmplitudes(amplitudes,*prodKinMomenta,*decayKinMomenta);
		if (binned){
			for (size_t m=0;m<nMatrix;++m){
				if (inBin[m]){
					matrix[m].addEventAmplitudes(ampl);
				};
			};
		}else{
			matrix[0].addEventAmplitudes(ampl); 
		};
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
