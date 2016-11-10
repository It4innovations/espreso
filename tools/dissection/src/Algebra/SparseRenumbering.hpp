/*! \file SparseRenumbering.cpp
    \brief tridiagonal factorization algorithm with Cuthill-McKee
    \author François-Xavier Roux, ONERA, Laboratoire Jacques-Louis Lions
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
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

#include <vector>
#include "Algebra/CSR_matrix.hpp"

using std::vector;

void CMK_number(const int dim, const int *ptrows, const int *indcols,
		vector<int> &new2old, const bool verbose, FILE *fp);

double frontal_numb(const int dim, const int *ptrows, const int *indcols,
		    const int i0, 
		    vector<int> &list,
		    vector<int> &indic, vector<int> &connect);

int point_front(const int dim, const int *ptrows, const int *indcols,
		vector<int> &new2old, vector<int> &p_front);

int getColorMaskCSR(int *color_mask, const CSR_indirect *csr,
		    const bool verbose,	FILE *fp);
