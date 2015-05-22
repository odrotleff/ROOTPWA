#include<iostream>
#include<string>
#include<vector>
#include<complex>

#include"TFile.h"
#include"TTree.h"
#include"TTreePerfStats.h"
#include"TClonesArray.h"

#include"waveDescription.h"
#include"isobarAmplitude.h"
#include"particleDataTable.h"
#include"eventMetadata.h"
#include"ampIntegralMatrix.h"
namespace rpwa {
	namespace hli {
		bool getIntegralsFromKeyFiles(						const std::string			integralName,
											const std::string			outFileName,
											const std::vector<std::string> 		&keyFiles, 
											const std::vector<std::string> 		&eventFiles, 
											const int 				maxNmbEvents = -1);

		std::vector<rpwa::isobarAmplitudePtr> getAmplitudesFromKeyFiles(const std::vector<std::string> &keyFiles);
		std::vector<std::string> waveNamesFromKeyFiles(const std::vector<std::string> &keyFiles,bool newConvention = false);

		std::vector<std::complex<double> > evaluateAmplitudes(std::vector<rpwa::isobarAmplitudePtr> &amplitudes, TClonesArray prodKinematics, TClonesArray decayKinematics);

		bool initAmplitudesKinematics(std::vector<rpwa::isobarAmplitudePtr> &amplitudes, std::vector<std::string> prodNames, std::vector<std::string> decayNames);

		bool calcTbinnedIntegralsFromEventTree(	const eventMetadata* 				eventMeta, 
							std::vector<isobarAmplitudePtr> 		&amplitudes, 
							std::vector<ampIntegralMatrix>			&matrix,
							std::vector<double>				tBinning				= std::vector<double>(0),
							const long int 					maxNmbEvents 				= -1, 
							const long int					startEvent				=  0,
							const std::string& 				treePerfStatOutFileName			= "", 
							const long int 					treeCacheSize 				= 25000000);
	};
};
