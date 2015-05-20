#include"tPrime_py.h"
#include"physUtils.hpp"

#include<TLorentzVector.h>
#include "rootConverters_py.h"
#include "stlContainers_py.h"


namespace bp = boost::python;

double tPrime_py(PyObject* beam, PyObject* target, PyObject* out){
	TLorentzVector* lvBeam   = rpwa::py::convertFromPy<TLorentzVector*>(beam  );
	TLorentzVector* lvTarget = rpwa::py::convertFromPy<TLorentzVector*>(target);
	TLorentzVector* lvOut    = rpwa::py::convertFromPy<TLorentzVector*>(out   );

	return rpwa::tPrime(*lvBeam,*lvTarget,*lvOut);
};


void rpwa::py::exportTprime() {
	bp::def("tPrime", &tPrime_py);
};
