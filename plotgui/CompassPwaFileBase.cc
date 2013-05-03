///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
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
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev:: 862                         $: revision of last commit
// $Author:: schmeing                 $: author of last commit
// $Date:: 2012-07-06 13:54:31 +0200 #$: date of last commit
//
// Description:
//      Code file for the CompassPwaFileBase class that provides
//		functionality to read in the next valid line of a txt file
//		generated by Compass pwa and passes it to all CompassPWAFile
//		classes
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#include <sstream>

#include "reportingUtils.hpp"

#include "CompassPwaFileBase.h"

using namespace std;
using namespace rpwa;

bool CompassPwaFileBase::_Debug = false;

// Constructor
CompassPwaFileBase::CompassPwaFileBase():
		_MassBinStart(0),
		_MassBinEnd(0),
		_tBinStart(0),
		_tBinEnd(0){
}

// Returns _MassBinStart;
double CompassPwaFileBase::MassBinStart() const{
	return _MassBinStart;
}

// Returns _MassBinEnd;
double CompassPwaFileBase::MassBinEnd() const{
	return _MassBinEnd;
}

// Returns _tBinStart;
double CompassPwaFileBase::tBinStart() const{
	return _tBinStart;
}

// Returns _tBinEnd;
double CompassPwaFileBase::tBinEnd() const{
	return _tBinEnd;
}
// Prints all important variables of class
std::ostream& CompassPwaFileBase::Print( std::ostream& Out ) const{
	Out << "Mass bin: " << _MassBinStart << ';' << _MassBinEnd << '\n';
	Out << "t' bin: " << _tBinStart << ';' << _tBinEnd << '\n';

	return Out;
}


// Reads the mass and t' bin from file
bool CompassPwaFileBase::ReadBin( istream& File ){
	bool Succesful = true; // Is set to false if an error occurs and returned at the end of the function
	stringstream LineStream;
	char semicolon = 0; // Takes the separating character, which should be a semicolon

	if( CompassPwaFileBase::GetNextValidLine( File, LineStream ) ){
		// Line example between "": "2.50000;2.51000"
		LineStream >> _MassBinStart >> semicolon >> _MassBinEnd;

		if( LineStream.fail() ){
			printErr << "Mass bin not valid \n";
			Succesful = false;
			if( _Debug ){
				printDebug << "MassBin: "<< LineStream.str() << '\n';
			}
		}
		else{
			if( !LineStream.eof() ){
				printWarn << "Mass bin entry longer than expected\n";
				if( _Debug ){
					printDebug << "MassBin: "<< LineStream.str() << '\n';
				}
			}
			if( semicolon != ';' ){
				printWarn << "Mass bin separator not a semicolon\n";
				if( _Debug ){
					printDebug << "Separator: '" << semicolon << "'\n";
				}
			}
		}
	}
	else{
		Succesful = false;
		printErr << "No valid line could be found anymore, but the mass bin was expected\n";
	}

	semicolon = 0;
	if( CompassPwaFileBase::GetNextValidLine( File, LineStream ) ){
		// Line example between "": "0.00000;1.00000"
		LineStream >> _tBinStart >> semicolon >> _tBinEnd;

		if( LineStream.fail() ){
			printErr << "t' bin not valid \n";
			Succesful = false;
			if( _Debug ){
				printDebug << "t 'Bin: "<< LineStream.str() << '\n';
			}
		}
		else{
			if( !LineStream.eof() ){
				printWarn << "t' bin entry longer than expected\n";
				if( _Debug ){
					printDebug << "t 'Bin: "<< LineStream.str() << '\n';
				}
			}
			if( semicolon != ';' ){
				printWarn << "t' bin separator not a semicolon\n";
				if( _Debug ){
					printDebug << "Separator: '" << semicolon << "'\n";
				}
			}
		}
	}
	else{
		Succesful = false;
		printErr << "No valid line could be found anymore, but the t' bin was expected\n";
	}

	return Succesful;
}

// Reads the next line of File which has no "//" as the first two letters and returns true if it was successful
bool CompassPwaFileBase::GetNextValidLine( istream& File, string& Line ){
	while( File.good() ){
		getline( File, Line );

		if( ( '/' != Line[0] ) || ( '/' != Line[1] ) ){
			return true;
		}
	}

	return false;
}

bool CompassPwaFileBase::GetNextValidLine( istream& File, stringstream& Line ){
	string TmpLine;

	if( GetNextValidLine( File, TmpLine ) ){
		Line.str( TmpLine );
		Line.clear();
		Line.seekg(0);
		return true;
	}
	else{
		return false;
	}
}
