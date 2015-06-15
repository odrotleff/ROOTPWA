#include"tPrime_py.h"
#include"physUtils.hpp"

#include<TLorentzVector.h>
#include "rootConverters_py.h"
#include "stlContainers_py.h"


namespace bp = boost::python;

double tPrime_py(bp::list beam, bp::list target, bp::list out){
	std::vector<double> convert;

	rpwa::py::convertBPObjectToVector<double>(beam, convert);
	TLorentzVector lvBeam = TLorentzVector(&convert[0]);

	rpwa::py::convertBPObjectToVector<double>(target, convert);
	TLorentzVector lvTarget = TLorentzVector(&convert[0]);

	rpwa::py::convertBPObjectToVector<double>(out, convert);
	TLorentzVector lvOut = TLorentzVector(&convert[0]);

	return rpwa::tPrime(lvBeam,lvTarget,lvOut);
};



void rpwa::py::exportTprime() {
	bp::def("tPrime", &tPrime_py);
};
