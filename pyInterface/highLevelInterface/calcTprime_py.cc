#include "calcTprime_py.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"

#include<iostream>
namespace bp = boost::python;

bp::list getTprimesFromFiles_py(	const std::string dataFile, const std::string keyFile, const long int maxnNmbEvents){

	return bp::list(rpwa::hli::getTprimesFromFiles(dataFile, keyFile,maxnNmbEvents));
};
					

void rpwa::py::exportCalcTprime(){
	bp::def(
		"getTprimesFromFiles"
		, &getTprimesFromFiles_py
	);
};
