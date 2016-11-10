/*! \file   TridiagMatrix.hpp
    \brief  management of threads for factorization and Fw/Bw substitution
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Jun. 20th 2014
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

# ifndef _DRIVER_TRIDIAGMATRIX_
# define _DRIVER_TRIDIAGMATRIX_
#include <complex>
#include <vector>
#include "Driver/C_threads_tasks.hpp"
#include "Algebra/SquareBlockMatrix.hpp"

using std::vector;

template<typename T>
class TridiagMatrix
{
public:
  TridiagMatrix(int dim, 
		int nnz, int *ptrows, int *indcols, int *indvals,
		int *remap_eqn,
		bool isSym,
		bool isMapped) :
    _nrow(dim), _nnz(nnz), _isSym(isSym), _isMapped(isMapped)
  {
    _ptrows = new int[dim + 1];
    _indcols = new int[nnz];
    _indvals = new int[nnz];
    _remap_eqn = new int[dim];
    
    for (int i = 0; i < (dim + 1); i++) {
      _ptrows[i] = ptrows[i] + 1;  // C to Fortran
    }
    for (int i = 0; i < nnz; i++) {
      _indcols[i] = indcols[i] + 1; // C to Fortran
      _indvals[i] = indvals[i] + 1; // C to Fortran
    }
    for (int i = 0; i < dim; i++) {
      _remap_eqn[i] = remap_eqn[i];
    }
    _diag =  new SquareBlockMatrix<T>; // to keep pivot information
  }
  /** Destructor.
   */
  ~TridiagMatrix()
  {
    //    FORTRAN_DECL(tridiag_free)(_tridiag_sparse);
    delete [] _ptrows;
    delete [] _indcols;
    delete [] _indvals;
    delete [] _remap_eqn;
    delete _diag;
  }

  int nrow() {
    return _nrow;
  }
  int nnz() {
    return _nnz;
  }
  void setGlobalNonzero(int nnz) {
    _nnz_global = nnz;
  }
  int nnz_global() {
    return _nnz_global;
  }
  int *ptRows() {
    return _ptrows;
  }
  int *indCols() {
    return _indcols;
  } 
  int *indVals() {
    return _indvals;
  }
  int *remap_eqn() {
    return _remap_eqn;
  }
  void setCoef(T *coef) {
    _coef = coef;
  }
  T *getCoef() {
    return _coef;
  }
  SquareBlockMatrix<T>* Diag() {
    return _diag;
  }
  bool isSym() const { return _isSym; }

  bool isMapped() const { return _isMapped; }
  
  void setPtr(void *tridiag_sparse)
  {
    _tridiag_sparse = tridiag_sparse;
  }
  void* &getPtr() {return _tridiag_sparse; }
  const void* &getPtr() const {return _tridiag_sparse; }
  
private:
  // Attributs :
  //  Dissection::Tree _btree;
  int _nrow;
  int _nnz;
  int _nnz_global;
  int *_ptrows;
  int *_indcols;
  int *_indvals;
  int *_remap_eqn;
  T* _coef;
  void *_tridiag_sparse;
  bool _isSym;
  bool _isMapped;
  SquareBlockMatrix<T> *_diag;
};

#endif
