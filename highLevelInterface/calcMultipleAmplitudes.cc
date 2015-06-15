#include"calcMultipleAmplitudes.h"
#include <boost/numeric/conversion/cast.hpp>
#include"calcAmplitude.h"

void rpwa::hli::firstTryAgain(){
	particleDataTable::readFile("/nfs/hicran/project/compass/analysis/fkrinner/ROOTPWA/particleData/particleDataTable.txt");
	TFile* inFile = new TFile("/nfs/mds/user/no72doz/3pi--+/data/1160_1170.generated.root","READ");
	const eventMetadata* data = eventMetadata::readEventFile(inFile);
	std::cout<<">>><<<"<<data<<std::endl;
	std::vector<std::string> keyFiles;
	keyFiles.push_back("/nfs/hicran/project/compass/analysis/fkrinner/ROOTPWA/userAnalysisWorkspace/3pi.--+/keyfiles/keyfiles/1-3+00+f0980_30_pi-.key");
	keyFiles.push_back("/nfs/hicran/project/compass/analysis/fkrinner/ROOTPWA/userAnalysisWorkspace/3pi.--+/keyfiles/keyfiles/1-2+02-f21270_12_pi-.key");
	keyFiles.push_back("/nfs/hicran/project/compass/analysis/fkrinner/ROOTPWA/userAnalysisWorkspace/3pi.--+/keyfiles/keyfiles/1-2-01+rho31690_33_pi-.key");
	keyFiles.push_back("/nfs/hicran/project/compass/analysis/fkrinner/ROOTPWA/userAnalysisWorkspace/3pi.--+/keyfiles/keyfiles/1-2-01+f01500_20_pi-.key");
	std::vector<rpwa::isobarAmplitudePtr> amplitudes = getAmplitudesFromKeyFiles(keyFiles);
	std::vector<std::string> names = waveNamesFromKeyFiles(keyFiles, true);
	std::vector<double> tBinning(0);
	std::vector<ampIntegralMatrix> mattt = std::vector<ampIntegralMatrix>(1,ampIntegralMatrix());
	size_t count =0;
	ofstream writeOut("out.txt");;
	for (size_t nEvt=1;nEvt<1414;++nEvt){
		for (size_t i=0;i<mattt.size();++i){
			mattt[i].initialize(names);
		};
		calcTbinnedIntegralsFromEventTree(data,amplitudes,mattt,tBinning,count+nEvt,count);
		count+=nEvt;
		writeOut<<nEvt<<" "<<std::real(mattt[0].element(3,0))<<" "<<std::imag(mattt[0].element(3,0))<<" "<<std::real(mattt[0].element(3,1))<<" "<<std::imag(mattt[0].element(3,1))<<" "<<std::real(mattt[0].element(3,2))<<" "<<std::imag(mattt[0].element(3,2))<<" "<<std::real(mattt[0].element(3,3))<<" "<<std::imag(mattt[0].element(3,3))<<std::endl;
	};
	writeOut.close();
//	mattt[0].writeAscii("./integralsSelf.txt");

	std::cout<<"Voilà!"<<std::endl;
};

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

std::vector<double> rpwa::hli::getTprimesFromFiles(					const std::string			fileName,
											const std::string			keyFile,
											const long int				maxNmbEvents){

	std::vector<std::string> vectorKey(1,keyFile);
	std::vector<isobarAmplitudePtr> amplitudes = getAmplitudesFromKeyFiles(vectorKey);
	TFile* inFile = new TFile(fileName.c_str(),"READ");
	const eventMetadata* actData = eventMetadata::readEventFile(inFile);

	return getTprimesFromTree(actData, amplitudes[0],maxNmbEvents);
};


std::vector<double> rpwa::hli::getTprimesFromTree(					const rpwa::eventMetadata* 		eventMeta,
											isobarAmplitudePtr			amplitude,
											const long int				maxNmbEvents,
											const std::string& 			treePerfStatOutFileName,
											const long int 				treeCacheSize){

	if (not amplitude){
		printWarn<<"no amplitude given, abort."<<std::endl;
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

	if(not amplitude->decayTopology()->initKinematicsData((*eventMeta).productionKinematicsParticleNames(), (*eventMeta).decayKinematicsParticleNames())){
		printWarn << "problems initializing input data. cannot read input data." << std::endl;
		return std::vector<double>();
	};
	const long int    nmbEventsTree     = tree->GetEntries();
	const long int    nmbEvents         = ((maxNmbEvents > 0) ? std::min(maxNmbEvents, nmbEventsTree) : nmbEventsTree);
	std::vector<double> ret(nmbEvents);
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
		if (not amplitude->decayTopology()->readKinematicsData(*prodKinMomenta,*decayKinMomenta)){
			printWarn<< "could not set kinematics for event "<<eventIndex<<std::endl;
			return ret;
		};
		ret[eventIndex] = amplitude->decayTopology()->getTprime();
		if (eventIndex%1000==0){
			printInfo<<eventIndex<<std::endl;
		};	
	};
	if(treePerfStats){
		treePerfStats->SaveAs(treePerfStatOutFileName.c_str());
		delete treePerfStats;
	};

	return ret;





};

