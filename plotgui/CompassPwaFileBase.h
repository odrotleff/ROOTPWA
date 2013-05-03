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
// $Rev:: 836                         $: revision of last commit
// $Author:: schmeing                 $: author of last commit
// $Date:: 2011-12-21 12:31:38 +0100 #$: date of last commit
//
// Description:
//      Header file for the CompassPwaFileBase class that provides
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

#ifndef CompassPwaFileBase_H
#define CompassPwaFileBase_H

#include <fstream>
#include <string>
#include <sstream>

namespace rpwa{

	class CompassPwaFileBase{
	private:
		// Variables
		double _MassBinStart;
		double _MassBinEnd;
		double _tBinStart;
		double _tBinEnd;

		static bool _Debug; ///< if set to true, debug messages are printed

		// Functions

	protected:
		// Variables

		// Functions

	public:
		// Constructors + Destructors
		CompassPwaFileBase(); ///< Constructor
		virtual ~CompassPwaFileBase(){}; ///< Virtual destructor

		// Get && Set
		double MassBinStart() const; ///< Returns _MassBinStart;
		double MassBinEnd() const; ///< Returns _MassBinEnd;
		double tBinStart() const; ///< Returns _tBinStart;
		double tBinEnd() const; ///< Returns _tBinEnd;

		static bool Debug() { return _Debug; } ///< returns debug flag
		static void SetDebug(const bool Debug = true) { _Debug = Debug; } ///< sets debug flag

		// Functions
		bool ReadBin( std::istream& File ); ///< Reads the mass and t' bin from file stream
		virtual bool ReadIn( std::istream& File ) = 0; ///< Reads the rest of the information from a specific file and returns 0 if no error occurred or a negative number as the error code
		virtual std::ostream& Print( std::ostream& Out ) const; ///< Prints all important variables of class

		static bool GetNextValidLine( std::istream& File, std::string& Line); ///< Reads the next line of File which has no "//" as the first two letters into Line and returns true if it was successful
		static bool GetNextValidLine( std::istream& File, std::stringstream& Line); ///< Reads the next line of File which has no "//" as the first two letters into Line and returns true if it was successful
	};

} // namespace rpwa

#endif /* CompassPwaFileBase_H */
