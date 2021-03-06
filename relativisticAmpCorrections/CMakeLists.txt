#///////////////////////////////////////////////////////////////////////////
#//
#//    Copyright 2010
#//
#//    This file is part of rootpwa
#//
#//    rootpwa is free software: you can redistribute it and/or modify
#//    it under the terms of the GNU General Public License as published by
#//    the Free Software Foundation, either version 3 of the License, or
#//    (at your option) any later version.
#//
#//    rootpwa is distributed in the hope that it will be useful,
#//    but WITHOUT ANY WARRANTY; without even the implied warranty of
#//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#//    GNU General Public License for more details.
#//
#//    You should have received a copy of the GNU General Public License
#//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
#//
#///////////////////////////////////////////////////////////////////////////
#//-------------------------------------------------------------------------
#//
#// Description:
#//      build file for amplitude library
#//
#//
#// Author List:
#//      Boris Grube          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


message_setup_this_dir()


# set include directories
set(INCLUDE_DIRECTORIES
	${CMAKE_CURRENT_SOURCE_DIR}
	${RPWA_DECAYAMPLITUDE_INCLUDE_DIR}
	${RPWA_UTILITIES_INCLUDE_DIR}
	${RPWA_PARTICLEDATA_INCLUDE_DIR}
	SYSTEM
	${Boost_INCLUDE_DIRS}
	${ROOT_INCLUDE_DIR}
	${Libconfig_INCLUDE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into library
set(SOURCES
	TFracNum.cc
	ClebschGordanBox.cc
	TJwfTensor.cc
	TSpinWaveFunction.cc
	TLSAmpl.cc
	TFhh.cc
	)


# library
set(THIS_LIB "RelAmpl")
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	)


# executables
make_executable(CalcAmpl  CalcAmpl.cc  ${THIS_LIB})
make_executable(CalcAmpl2 CalcAmpl2.cc ${THIS_LIB} ${Libconfig_LIBS} ${RPWA_DECAYAMPLITUDE_LIB} ${RPWA_PARTICLEDATA_LIB})
