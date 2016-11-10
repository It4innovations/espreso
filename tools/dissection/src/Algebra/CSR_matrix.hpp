/*! \file CSR_matrix.hpp
    \brief  Sparse matrix data structure
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Mar. 30th 2012
    \date   Jul. 12th 2015
    \date   Feb. 29th 2016
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

#ifndef _DISSECTION_SPLITTERS_CSR_MATRIX_HPP_
#define _DISSECTION_SPLITTERS_CSR_MATRIX_HPP_

// #define DEBUG_MAPPING_CSR

struct CSR_indirect {
  int n;
  int nnz;
  int *ptRows;
  int *indCols;
  int *indVals0;  // sym for debugging
  int *indVals;   // sym for after the first permutation to exclude diagonals
  int *indVals_unsym;
  bool isSym;
#ifdef DEBUG_MAPPING_CSR
  int *indVals2; // unsym
#endif
};
#endif
