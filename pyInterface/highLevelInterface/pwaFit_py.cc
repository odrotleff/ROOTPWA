
#include "pwaFit_py.h"

#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;

namespace {

	rpwa::fitResultPtr pwaFit(PyObject* ampTreesDict,
	                          PyObject* normMatrix,
	                          PyObject* accMatrix,
	                          double massBinMin = 0.,
	                          double massBinMax = 0.,
	                          std::string waveListFileName = "",
	                          std::string startValFileName = "",
	                          bool useNormalizedAmps = false,
	                          unsigned int rank = 1)
	{
		std::map<std::string, TTree*> ampTreesMap;
		bp::dict pyDictAmpTreesDict = bp::extract<bp::dict>(ampTreesDict);
		bp::list keys = pyDictAmpTreesDict.keys();
		for(int i = 0; i < bp::len(keys); ++i) {
			bp::object curTree = pyDictAmpTreesDict[keys[i]];
			if(curTree) {
				ampTreesMap.insert(std::pair<std::string, TTree*>(bp::extract<std::string>(keys[i]), rpwa::py::convertFromPy<TTree*>(curTree.ptr())));
			}
		}
		rpwa::ampIntegralMatrix* normMatrixPtr = rpwa::py::convertFromPy(normMatrix);
		rpwa::ampIntegralMatrix* accMatrixPtr = rpwa::py::convertFromPy(accMatrix);
		return rpwa::hli::pwaFit(ampTreesMap, *normMatrixPtr, *accMatrixPtr, massBinMin, massBinMax, waveListFileName, startValFileName, useNormalizedAmps, rank);
//		return rpwa::fitResultPtr(new rpwa::fitResult());
	}

}

void rpwa::py::exportPwaFit()
{

	bp::def(
		"pwaFit"
		, &::pwaFit
		, (bp::arg("ampTreesDict"),
		   bp::arg("normMatrix"),
		   bp::arg("accMatrix"),
		   bp::arg("massBinMin") = 0,
		   bp::arg("massBinMax") = 0,
		   bp::arg("waveListFileName") = "",
		   bp::arg("startValFileName") = "",
		   bp::arg("useNormalizedAmps") = false,
		   bp::arg("rank") = 1)
	);

}
