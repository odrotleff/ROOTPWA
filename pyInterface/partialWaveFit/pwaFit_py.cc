
#include "pwaFit_py.h"

#include <TTree.h>

#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;

namespace {

	int pwaFit_fit(rpwa::pwaFit& self,
	               const bp::dict& ampTreesDict,
	               const rpwa::ampIntegralMatrix& normMatrix,
	               const rpwa::ampIntegralMatrix& accMatrix,
	               double massBinMin = 0,
	               double massBinMax = 0,
	               std::string waveListFileName = "",
	               std::string startValFileName = "",
	               bool useNormalizedAmps = false,
	               unsigned int rank = 1)
	{
		std::map<std::string, TTree*>& ampTreesMap;

		bp::list keys = ampTreesDict.keys();
		for(int i = 0; i < len(keys); ++i) {
			bp::object curTree = ampTreesDict[keys[i]];
			if(curTree) {
				ampTreesMap.insert(bp::extract<std::string>(keys[i]), rpwa::py::convertFromPy<TTree*>(ampTreesDict[keys[i]]));
			}
		}
		return self.fit();
	}
}

void rpwa::py::exportFitResult() {

	bp::class_<rpwa::pwaFit>("pwaFit")
		.def("evidenceComponents", &pwaFit_fit);
}
