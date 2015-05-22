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
		std::vector<double> getTprimesFromTree(					const eventMetadata* 			eventMeta,
											isobarDecayTopologyPtr			topology,
											const long int				maxNmbEvents = -1,
											const std::string& 			treePerfStatOutFileName = "",
											const long int 				treeCacheSize = 25000000);


		std::vector<double> getTprimesFromFiles(				const std::string			fileName,
											const std::string			keyFile,
											const long int				maxNmbEvents);
	};
};
