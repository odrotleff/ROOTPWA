#include "calcMultipleAmplitudes_py.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"

#include<iostream>
namespace bp = boost::python;
void write_voila(){
	rpwa::hli::firstTryAgain();
};

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

bp::list getTprimesFromFiles_py(	const std::string dataFile, const std::string keyFile, const long int maxnNmbEvents){

	return bp::list(rpwa::hli::getTprimesFromFiles(dataFile, keyFile,maxnNmbEvents));
};
					

void rpwa::py::exportCalcMultipleAmplitudes(){
	bp::def(
		"calcMultipleAmplitudes"
		, &write_voila
	);
	

	bp::def(
		"getIntegralsFromKeyFiles"
		, & getIntegralsFromKeyFiles_py
	);

	bp::def(
		"getTprimesFromFiles"
		, &getTprimesFromFiles_py
	);
};
