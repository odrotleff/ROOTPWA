
namespace rpwa {

	class pwaFit {

	public:
		int fit(
		        std::map<std::string, TTree*>& ampTrees,
		        rpwa::ampIntegralMatrix& normMatrix,
		        rpwa::ampIntegralMatrix& accMatrix,
		        double massBinMin = 0,
		        double massBinMax = 0,
		        std::string waveListFileName = "",
		        std::string startValFileName = "",
		        bool useNormalizedAmps = false,
		        unsigned int rank = 1);

		static const std::string valTreeName;
		static const std::string valBranchName;
	};

}
