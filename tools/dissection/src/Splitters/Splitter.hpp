/*! \file   Splitter.hpp
    \brief  pointer for generic grpah decomposer
    \author Xavier Juvigny, ONERA
    \date   Jul.  2nd 2012
*/

// This file is part of Dissection
// 
// Dissection is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Dissection is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Dissection.  If not, see <http://www.gnu.org/licenses/>.

// ==============================================================
// ==     Definition of pointer function for splitter tool     ==
// ==============================================================
#ifndef _DISSECTION_SPLITTERS_SPLITTER_HPP_
#define _DISSECTION_SPLITTERS_SPLITTER_HPP_
typedef bool (*splitter)(unsigned dim, 
			 const int* ptRows, const int* indCols, 
			 unsigned& nbMaxLevels, unsigned minSize,
			 int* loc2glob, int* glob2loc,
			 int& nbDoms, int*& ptOnDomains,
			 int*& sizeOfDomains, bool checkData,
			 const bool verbose,
			 FILE *fp);
#endif
