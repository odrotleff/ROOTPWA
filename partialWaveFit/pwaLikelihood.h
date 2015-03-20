///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------
//
// Description:
//      Likelihood function Object to use with ROOT minimizers
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#ifndef PWALIKELIHOOD_H
#define PWALIKELIHOOD_H


#include <vector>
#include <string>
#include <iostream>

#define BOOST_DISABLE_ASSERTS
#include "boost/multi_array.hpp"
#include "boost/tuple/tuple.hpp"

#include "Math/IFunction.h"
#include "TFile.h"
#include "TH1.h"

#include "sumAccumulators.hpp"
#include "ampIntegralMatrix.h"


class TString;


namespace rpwa {


	class complexMatrix;


	template<typename complexT>  // type of internal variables used for intermediate results
	class pwaLikelihood : public ROOT::Math::IGradientFunctionMultiDim {

	public:

		typedef typename complexT::value_type value_type;

		// define array types
		typedef boost::multi_array<std::string,                    2> waveNameArrayType;    // array for wave names
		typedef boost::multi_array<double,                         2> waveThrArrayType;     // array for wave thresholds
		typedef boost::multi_array<unsigned int,                   2> waveToIntMapType;     // array for mapping of waves to integral indices
		typedef boost::multi_array<unsigned int,                   2> waveToListMapType;    // array for mapping of waves to position in wave list
		typedef boost::multi_array<boost::tuples::tuple<int, int>, 3> ampToParMapType;      // array for mapping of amplitudes to parameters
		typedef boost::multi_array<complexT,                       3> ampsArrayType;        // array for production and decay amplitudes
		typedef boost::multi_array<complexT,                       4> normMatrixArrayType;  // array for normalization matrices
		typedef boost::multi_array<value_type,                     2> phaseSpaceIntType;    // array for phase space integrals


	public:

		// enum for function call counters
		enum functionCallEnum {
			FDF                  = 0,
			GRADIENT             = 1,
			DOEVAL               = 2,
			DODERIVATIVE         = 3,
			NMB_FUNCTIONCALLENUM = 4
		};

		struct functionCallInfo {
			typedef boost::accumulators::accumulator_set
			  <double, boost::accumulators::stats
				<boost::accumulators::tag::sum(boost::accumulators::compensated)> > timeAccType;
			unsigned int nmbCalls;   // number of times function was called
			timeAccType  funcTime;   // time needed to calculate function value(s) (w/o normalization)
			timeAccType  normTime;   // time needed to normalize function value(s)
			timeAccType  totalTime;  // total execution time of function
		};

		pwaLikelihood();
		~pwaLikelihood();

		// overload public IGradientFunctionMultiDim member functions:
		/// clones the function using the default copy constructor
		virtual pwaLikelihood* Clone() const { return new pwaLikelihood(*this); }
		/// returns total number of function parameters (= dimension of the function)
		virtual unsigned int NDim() const { return nmbPars(); }
		/// optimized method to evaluate function value and derivative at a point defined by par at the same time
		virtual void FdF(const double* par,
						 double&       funcVal,
						 double*       gradient) const;
		/// calculates gradient (vector of partial derivatives) of function at point defined by par
		virtual void Gradient(const double* par,
							  double*       gradient) const;

		// overload private IGradientFunctionMultiDim member functions
		virtual double DoEval      (const double* par) const;
		virtual double DoDerivative(const double* par,
									unsigned int  derivativeIndex) const;

		unsigned int             nmbEvents   ()                                    const { return _nmbEvents;               }  ///< returns number of events that enter in the likelihood
		unsigned int             rank        ()                                    const { return _rank;                    }  ///< returns rank of spin density matrix
		unsigned int             nmbWaves    (const int          reflectivity = 0) const;                                      ///< returns total number of waves (reflectivity == 0) or number or number of waves with positive/negative reflectivity; flat wave is not counted!
		unsigned int             nmbPars     ()                                    const { return _nmbPars;                 }  ///< returns total number of parameters
		// std::string              waveName    (const unsigned int waveIndex)        const { return _waveNames[waveIndex];    }  ///< returns name of wave at waveIndex
		std::vector<std::string> waveNames   ()                                    const;  ///< returns vector with all wave names ordered like in input wave list
		std::string              parName     (const unsigned int parIndex)         const { return _parNames[parIndex];      }  ///< returns name of likelihood parameter at parIndex
		double                   parThreshold(const unsigned int parIndex)         const { return _parThresholds[parIndex]; }  ///< returns threshold in GeV/c^2 above which likelihood parameter at parIndex becomes free

		double dLcache(const unsigned int i) const { return _derivCache[i]; }
		unsigned int ncalls(const functionCallEnum func = FDF) const
		{ return _funcCallInfo[func].nmbCalls; }
		double Ltime(const functionCallEnum func = FDF) const
		{ return boost::accumulators::sum(_funcCallInfo[func].funcTime); }
		double Ntime(const functionCallEnum func = FDF) const
		{ return boost::accumulators::sum(_funcCallInfo[func].normTime); }
		//const integral& normInt() const { return _normInt; }

		// modifiers
		void        enableCuda       (const bool enableCuda = true);
		bool        cudaEnabled      () const;
		void        useNormalizedAmps(const bool useNorm = true) { _useNormalizedAmps = useNorm; }
		static void setQuiet         (const bool flag    = true) { _debug             = !flag;   }

		// operations
		void init(const unsigned int                   rank,
		          const std::string&                   waveListFileName,
		          const ampIntegralMatrix&             normMatrix,
		          const ampIntegralMatrix&             accMatrix,
		          std::map<std::string, TTree*>& ampTrees,
		          const unsigned int                   numbAccEvents);  ///< prepares all internal data structures

		void getIntegralMatrices(rpwa::complexMatrix&       normMatrix,
								 rpwa::complexMatrix&       accMatrix,
								 std::vector<double>&       phaseSpaceIntegral) const;

		// note: amplitudes which do not exist in higher ranks are NOT built!
		void buildProdAmpArrays(const double*                       inPar,
								std::vector<std::complex<double> >& prodAmps,
								std::vector<std::pair<int,int> >&   parIndices,
								std::vector<std::string>&           prodAmpNames,
								const bool                          withFlat = false) const;

		std::ostream& print(std::ostream& out = std::cout) const;
		std::ostream& printFuncInfo(std::ostream& out = std::cout) const;
		friend std::ostream& operator << (std::ostream&         out,
		                                  const pwaLikelihood& func) { return func.print(out); }

		std::vector<unsigned int> orderedParIndices() const;  // helper function for backwards-compatibility


	private:

		// helper functions
		void readWaveList       (const std::string& waveListFileName);  ///< reads wave names and thresholds from wave list file
		void buildParDataStruct (const unsigned int rank);              ///< builds parameter data structures
		void readIntegrals      (const ampIntegralMatrix& normIntFileName,
								 const ampIntegralMatrix& accIntFileName);  ///< reads normalization and acceptance integrals from file

		void readDecayAmplitudes(std::map<std::string, TTree*>& trees);  ///< reads decay amplitudes from files in specified directory


		void clear();
		static int getReflectivity(const TString& waveName);

		void reorderIntegralMatrix(const rpwa::ampIntegralMatrix& integral,
								   normMatrixArrayType&           reorderedMatrix) const;

	public:

		void copyFromParArray(const double*  inPar,              // input parameter array
							  ampsArrayType& outVal,             // output values organized as 3D array of complex numbers with [rank][reflectivity][wave index]
							  value_type&    outFlatVal) const;  // output value corresponding to flat wave
		void copyToParArray(const ampsArrayType& inVal,          // values corresponding to production amplitudes [rank][reflectivity][wave index]
							const value_type     inFlatVal,      // value corresponding to flat wave
							double*              outPar) const;  // output parameter array

	private:

		void resetFuncCallInfo() const;

		unsigned int _nmbEvents;        // number of events
		unsigned int _rank;             // rank of spin density matrix
		unsigned int _nmbWaves;         // number of waves
		unsigned int _nmbWavesRefl[2];  // number of negative (= 0) and positive (= 1) reflectivity waves
		unsigned int _nmbWavesReflMax;  // maximum of number of negative and positive reflectivity waves
		unsigned int _nmbPars;          // number of function parameters

	#ifdef USE_CUDA
		bool        _cudaEnabled;        // if true CUDA kernels are used for some calculations
	#endif
		bool        _useNormalizedAmps;  // if true normalized amplitudes are used
		static bool _debug;              // if true debug messages are printed

		unsigned int _numbAccEvents; // number of input events used for acceptance integrals (accepted + rejected!)
		double       _totAcc;        // total acceptance in this bin

		waveNameArrayType        _waveNames;            // wave names [reflectivity][wave index]
		waveThrArrayType         _waveThresholds;       // mass thresholds of waves
		waveToListMapType        _waveToWaveIndex;      // maps wave to its index in wave list
		std::vector<std::string> _parNames;             // function parameter names
		std::vector<double>      _parThresholds;        // mass thresholds of parameters
		ampToParMapType          _prodAmpToFuncParMap;  // maps each production amplitude to the indices
														// of its real and imginary part in the parameter
														// array; negative indices mean that the parameter
														// is not existing due to rank restrictions

		ampsArrayType _decayAmps;  // precalculated decay amplitudes [event index][reflectivity][wave index]

		mutable std::vector<double> _parCache;    // parameter cache for derivative calc.
		mutable std::vector<double> _derivCache;  // cache for derivatives

		// normalization integrals
		normMatrixArrayType _normMatrix;          // normalization matrix w/o acceptance [reflectivity 1][wave index 1][reflectivity 2][wave index 2]
		normMatrixArrayType _accMatrix;           // normalization matrix with acceptance [reflectivity 1][wave index 1][reflectivity 2][wave index 2]
		phaseSpaceIntType   _phaseSpaceIntegral;  // phase space integrals

		mutable functionCallInfo _funcCallInfo[NMB_FUNCTIONCALLENUM];  // collects function call statistics

	};


}

#endif  // TPWALIKELIHOOD_H
