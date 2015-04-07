#include <TTree.h>

#include <ampIntegralMatrix.h>
#include <fitResult.h>

namespace rpwa {

	namespace hli {

		rpwa::fitResultPtr pwaFit(std::map<std::string, TTree*>& ampTrees,
		                          const rpwa::ampIntegralMatrix& normMatrix,
		                          const rpwa::ampIntegralMatrix& accMatrix,
		                          const double massBinMin,
		                          const double massBinMax,
		                          const std::string waveListFileName,
		                          const std::string startValFileName,
		                          const bool useNormalizedAmps,
		                          const unsigned int rank);
	}

}
