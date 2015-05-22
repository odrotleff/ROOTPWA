#include "calcMultipleAmplitudes_py.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"

#include<iostream>
namespace bp = boost::python;
bool getIntegralsFromKeyFiles_py(	const std::string	integralName,
					const std::string	outFileName,
					const bp::object	keyFilesPy,
					const bp::object	eventFilesPy,
					int			maxNmbEnvents){

		std::vector<std::string> vectorKeyFiles;
		if(not rpwa::py::convertBPObjectToVector<std::string>(keyFilesPy, vectorKeyFiles)) {
			PyErr_SetString(PyExc_TypeError, "could not extract waveNames");
			bp::throw_error_already_set();
		};
		std::vector<std::string> vectorEventFiles;
		if(not rpwa::py::convertBPObjectToVector<std::string>(eventFilesPy, vectorEventFiles)) {
			PyErr_SetString(PyExc_TypeError, "could not extract waveNames");
			bp::throw_error_already_set();
		};
		return rpwa::hli::getIntegralsFromKeyFiles(integralName, outFileName, vectorKeyFiles ,vectorEventFiles ,maxNmbEnvents);
};


void rpwa::py::exportCalcMultipleAmplitudes(){
	bp::def(
		"getIntegralsFromKeyFiles"
		, & getIntegralsFromKeyFiles_py
	);
};
