#include "pwaFit.h"

#include <complex>

#include <boost/progress.hpp>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TStopwatch.h"

#include "pwaLikelihood.h"

#include <reportingUtils.hpp>

const std::string valTreeName   = "pwa";
const std::string valBranchName = "fitResult_v2";

rpwa::fitResultPtr
rpwa::hli::pwaFit(std::map<std::string, TTree*>& ampTrees,
	   const rpwa::ampIntegralMatrix& normMatrix,
	   const rpwa::ampIntegralMatrix& accMatrix,
	   const double massBinMin,
	   const double massBinMax,
	   const std::string waveListFileName,
	   const std::string startValFileName,
	   const bool useNormalizedAmps,
	   const unsigned int rank)
{
//	rpwa::printCompilerInfo();
//	rpwa::printLibraryInfo ();
//	rpwa::printGitHash     ();
//	std::cout << std::endl;

	// force loading predefined std::complex dictionary
	// see http://root.cern.ch/phpBB3/viewtopic.php?f=5&t=9618&p=50164
//	gROOT->ProcessLine("#include <complex>");

	// ---------------------------------------------------------------------------

	double             defaultStartValue     = 0.01;
	bool               useFixedStartValues   = false;
	double             startValStep          = 0.0005;
	const unsigned int maxNmbOfIterations    = 20000;
	const unsigned int maxNmbOfFunctionCalls = 40000;
	const bool         runHesse              = true;
	const bool         runMinos              = false;
	int                startValSeed          = 1234567;

	unsigned int numbAccEvents            = 0;                      // number of events used for acceptance integrals
	std::string  minimizerType[2]         = {"Minuit2", "Migrad"};  // minimizer, minimization algorithm
	int          minimizerStrategy        = 1;                      // minimizer strategy
	double       minimizerTolerance       = 1e-10;                  // minimizer tolerance
#if ROOT_VERSION_CODE >= ROOT_VERSION(5, 34, 19)
//	bool         saveMinimizerMemory      = true;
#endif
//	bool         cudaEnabled              = false;                  // if true CUDA kernels are activated
	bool         quiet                    = false;

	// report parameters
	/*printInfo << "running pwaFit with the following parameters:" << std::endl;
	std::cout << "    mass bin [" <<  massBinMin << ", " <<  massBinMax << "] MeV/c^2" << std::endl
	          << "    path to wave list file ......................... '" << waveListFileName << "'" << std::endl
	          << "    path to amplitude directory .................... '" << ampDirName       << "'" << std::endl
	          << "    use .root amplitude files ...................... "  << yesNo(useRootAmps)      << std::endl
	          << "    path to output file ............................ '" << outFileName      << "'" << std::endl
	          << "    path to file with start values ................. '" << startValFileName << "'" << std::endl
	          << "    seed for random start values ................... "  << startValSeed            << std::endl;
	if (useFixedStartValues)
		std::cout << "    using fixed instead of random start values ..... " << defaultStartValue << std::endl
	              << "    use normalization .............................. "  << rpwa::yesNo(useNormalizedAmps) << std::endl
	              << "        path to file with normalization integral ... '" << normIntFileName  << "'" << std::endl
	              << "        path to file with acceptance integral ...... '" << accIntFileName   << "'" << std::endl
	              << "        number of acceptance norm. events .......... "  << numbAccEvents    << std::endl
	              << "    rank of spin density matrix .................... "  << rank                    << std::endl
	              << "    minimizer ...................................... "  << minimizerType[0] << ", " << minimizerType[1] << std::endl
	              << "    minimizer strategy ............................. "  << minimizerStrategy  << std::endl
	              << "    minimizer tolerance ............................ "  << minimizerTolerance << std::endl
	              << "    CUDA acceleration .............................. "  << enDisabled(cudaEnabled) << std::endl
	              << "    quiet .......................................... "  << rpwa::yesNo(quiet) << std::endl;*/

	// ---------------------------------------------------------------------------
	// setup likelihood function
	printInfo << "creating and setting up likelihood function" << std::endl;
	rpwa::pwaLikelihood<std::complex<double> > L;
	if (quiet)
		L.setQuiet();
	L.useNormalizedAmps(useNormalizedAmps);
#ifdef USE_CUDA
//	L.enableCuda(cudaEnabled);
#endif
	L.init(rank, waveListFileName, normMatrix, accMatrix, ampTrees, numbAccEvents);
	if (not quiet)
		std::cout << L << std::endl;
	const unsigned int nmbPar  = L.NDim();
	const unsigned int nmbEvts = L.nmbEvents();

	// ---------------------------------------------------------------------------
	// setup minimizer
	printInfo << "creating and setting up minimizer '" << minimizerType[0] << "' "
	          << "using algorithm '" << minimizerType[1] << "'" << std::endl;
	ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer(minimizerType[0], minimizerType[1]);
	if (not minimizer) {
		printErr << "could not create minimizer. exiting." << std::endl;
		throw;
	}

	// special for Minuit2
#if ROOT_VERSION_CODE >= ROOT_VERSION(5, 34, 19)
//	if(saveMinimizerMemory and dynamic_cast<ROOT::Minuit2::Minuit2Minimizer*>(minimizer)) {
//			((ROOT::Minuit2::Minuit2Minimizer*)minimizer)->SetStorageLevel(0);
//			printInfo << "Minuit2 storage level set to 0." << std::endl;
//	}
#endif

	minimizer->SetFunction        (L);
	minimizer->SetStrategy        (minimizerStrategy);
	minimizer->SetTolerance       (minimizerTolerance);

	// setting the ErrorDef to 1 since the ROOT interface does not
	// Propagate the value. Will do the error rescaling by hand below.
	minimizer->SetErrorDef(1);
	minimizer->SetPrintLevel      ((quiet) ? 0 : 3);
	minimizer->SetMaxIterations   (maxNmbOfIterations);
	minimizer->SetMaxFunctionCalls(maxNmbOfFunctionCalls);

	// ---------------------------------------------------------------------------
	// read in fitResult with start values
	printInfo << "reading start values from '" << startValFileName << "'" << std::endl;
	const double massBinCenter  = (massBinMin + massBinMax) / 2;
	rpwa::fitResult*   startFitResult = NULL;
	bool         startValValid  = false;
	TFile*       startValFile   = NULL;
	if (startValFileName.length() <= 2)
		printWarn << "start value file name '" << startValFileName << "' is invalid. "
		          << "using default start values." << std::endl;
	else {
		// open root file
		startValFile = TFile::Open(startValFileName.c_str(), "READ");
		if (not startValFile or startValFile->IsZombie())
			printWarn << "cannot open start value file '" << startValFileName << "'. "
			          << "using default start values." << std::endl;
		else {
			// get tree with start values
			TTree* tree;
			startValFile->GetObject(valTreeName.c_str(), tree);
			if (not tree)
				printWarn << "cannot find start value tree '"<< valTreeName << "' in file "
				          << "'" << startValFileName << "'" << std::endl;
			else {
				startFitResult = new rpwa::fitResult();
				tree->SetBranchAddress(valBranchName.c_str(), &startFitResult);
				// find tree entry which is closest to mass bin center
				unsigned int bestIndex = 0;
				double       bestMass  = 0;
				for (unsigned int i = 0; i < tree->GetEntriesFast(); ++i) {
					tree->GetEntry(i);
					if (fabs(massBinCenter - startFitResult->massBinCenter()) <= fabs(massBinCenter - bestMass)) {
						bestIndex = i;
						bestMass  = startFitResult->massBinCenter();
					}
				}
				tree->GetEntry(bestIndex);
				startValValid = true;
			}
		}
	}

	// ---------------------------------------------------------------------------
	// set start parameter values
	printInfo << "setting start values for " << nmbPar << " parameters" << std::endl
	          << "    parameter naming scheme is: V[rank index]_[IGJPCMR][isobar spec]" << std::endl;
	unsigned int maxParNameLength = 0;       // maximum length of parameter names
	std::vector<bool> parIsFixed(nmbPar, false);  // memorizes state of variables; ROOT::Math::Minimizer has no corresponding accessor
	{
		// use local instance of random number generator so that other
		// code has no chance of tampering with gRandom and thus cannot
		// affect the reproducability of the start values
		TRandom3 random(startValSeed);
		bool     success = true;
		for (unsigned int i = 0; i < nmbPar; ++i) {
			const std::string parName = L.parName(i);
			if (parName.length() > maxParNameLength)
				maxParNameLength = parName.length();
			// workaround, because Minuit2Minimizer::SetVariable() expects
			// that variables are set consecutively. how stupid is that?
			// so we prepare variables here and set values below
			if ((L.parThreshold(i) == 0) or (L.parThreshold(i) < massBinCenter)) {
				if (not minimizer->SetVariable(i, parName, 0, startValStep))
					success = false;
			} else {
				if (not minimizer->SetFixedVariable(i, parName, 0.))  // fix this parameter to 0
					success = false;
				parIsFixed[i] = true;
			}
		}
		const double         sqrtNmbEvts = sqrt((double)nmbEvts);
		std::vector<unsigned int> parIndices  = L.orderedParIndices();
		for (unsigned int i = 0; i < parIndices.size(); ++i) {
			const unsigned int parIndex = parIndices[i];
			double             startVal;
			const std::string       parName = minimizer->VariableName(parIndex);
			if (parName != L.parName(parIndex)) {
				printWarn << "parameter name in minimizer and likelihood is inconsistent "
				          << "(" << parName << " vs. " << L.parName(parIndex) << ")" << std::endl;
				success = false;
			}
			if (startValValid) {
				// get parameter value from fitResult
				assert(startFitResult);
				startVal = startFitResult->fitParameter(parName.c_str());
			} else
				startVal = (useFixedStartValues) ? defaultStartValue
					: random.Uniform(defaultStartValue, sqrtNmbEvts);
			// check if parameter needs to be fixed because of threshold
			if ((L.parThreshold(parIndex) == 0) or (L.parThreshold(parIndex) < massBinCenter)) {
				if (startVal == 0) {
					std::cout << "    read start value 0 for parameter " << parName << ". "
					     << "using default start value." << std::endl;
					startVal = (useFixedStartValues) ? defaultStartValue
						: random.Uniform(defaultStartValue, sqrtNmbEvts);
				}
				std::cout << "    setting parameter [" << std::setw(3) << i << "] "
				     << std::setw(maxParNameLength) << parName << " = " << rpwa::maxPrecisionAlign(startVal) << std::endl;
				if (not minimizer->SetVariableValue(parIndex, startVal))
					success = false;
			} else {
				std::cout << "    fixing parameter  [" << std::setw(3) << i << "] "
				     << std::setw(maxParNameLength) << parName << " = 0" << std::endl;
				if (not minimizer->SetVariableValue(parIndex, 0.))  // fix this parameter to 0
					success = false;
				parIsFixed[parIndex] = true;
			}
			if (not success) {
				printErr << "something went wrong when setting log likelihood parameters. aborting." << std::endl;
				throw;
			}
		}
		// cleanup
		if(startValFile) {
			startValFile->Close();
			delete startValFile;
			startValFile = NULL;
		}
	}

	// ---------------------------------------------------------------------------
	// find minimum of likelihood function
	//bool converged = false;
	//bool hasHesse  = false;
	printInfo << "performing minimization" << std::endl;
	{
		TStopwatch timer;
		timer.Start();
		bool success = minimizer->Minimize();
		timer.Stop();
		if (success)
			printInfo << "minimization finished successfully. " << std::flush;
		else
			printWarn << "minimization failed. " << std::flush;
		std::cout << "used " << std::flush;
		timer.Print();
		//converged = success;
		//printInfo << *minimizer;
		// printInfo << "covariance matrix:" <<endl;
		// for(unsigned int i = 0; i < nmbPar; ++i)
		// 	for(unsigned int j = 0; j < nmbPar; ++j)
		// 		std::cout << "    [" << i << "][" << j << "] = " << minimizer->CovMatrix(i, j) << std::endl;
		if (runHesse) {
			printInfo << "calculating Hessian matrix" << std::endl;
			timer.Start();
			success = minimizer->Hesse();
			timer.Stop();
			if (success)
				printInfo << "successfully calculated Hessian matrix. " << std::flush;
			else
				printWarn << "calculation of Hessian matrix failed. " << std::flush;
			std::cout << "used " << std::flush;
			timer.Print();
			//hasHesse = success;
			//printInfo << *minimizer;
			// printInfo << "covariance matrix:" <<endl;
			// for(unsigned int i = 0; i < nmbPar; ++i)
			// 	for(unsigned int j = 0; j < nmbPar; ++j)
			// 		std::cout << "    [" << i << "][" << j << "] = " << minimizer->CovMatrix(i, j) << std::endl;
		}
	}

	// ---------------------------------------------------------------------------
	// print results
	printInfo << "minimization result:" << std::endl;
	std::vector<unsigned int> parIndices = L.orderedParIndices();
	for (unsigned int i = 0; i< parIndices.size(); ++i) {
		const unsigned int parIndex = parIndices[i];
		std::cout << "    parameter [" << std::setw(3) << i << "] "
		     << std::setw(maxParNameLength) << L.parName(parIndex) << " = ";
		if (parIsFixed[parIndex])
			std::cout << minimizer->X()[parIndex] << " (fixed)" << std::endl;
		else {
			std::cout << std::setw(12) << rpwa::maxPrecisionAlign(minimizer->X()     [parIndex]) << " +- "
			     << std::setw(12) << rpwa::maxPrecisionAlign(minimizer->Errors()[parIndex]);
			if (runMinos) {
				double minosErrLow = 0;
				double minosErrUp  = 0;
				const bool success = minimizer->GetMinosError(parIndex, minosErrLow, minosErrUp);
				if (success)
					std::cout << "    Minos: " << "[" << minosErrLow << ", +" << minosErrUp << "]" << std::endl;
			} else
				std::cout << std::endl;
		}
	}
	printInfo << "function call summary:" << std::endl;
	L.printFuncInfo(std::cout);
#ifdef USE_CUDA
	printInfo << "total CUDA kernel time: "
	          << cuda::likelihoodInterface<cuda::complex<double> >::kernelTime() << " sec" << std::endl;
#endif


	if (minimizer)
		delete minimizer;
	return rpwa::fitResultPtr(new rpwa::fitResult());
}
