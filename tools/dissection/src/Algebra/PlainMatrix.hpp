/*! \file   PlainMatrix.hpp
    \brief  Common interface for plain matrices (square, rectangulars,...)
    \author Xavier Juvigny, ONERA
    \date   Jan. 25th 2005
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Jun. 11th 2013
    \date   Jul. 12th 2015
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

#ifndef _ALGEBRA_PLAINMATRIX_
# define _ALGEBRA_PLAINMATRIX_
#include <cstdio>
#include <cstring>
#include <new>
//#define DEBUG_MEMORY_PLAINMATRIX
#ifdef DEBUG_MEMORY_PLAINMATRIX
# include <cstdlib>
#include <stdint.h>
#endif

template<typename T>
class PlainMatrix
{
public:

  PlainMatrix(bool isOwner = true)
  {
    _isOwner = isOwner;
    _n = (-1);  // not yet allocated
    _coefs = new T*;
    _coefs_alloc_status = true;
  }

  PlainMatrix(int nbElts)
  {
    _isOwner = true;
    _n = nbElts;
    _coefs = new T*;
    _coefs_alloc_status = true;
    try {
      *_coefs = new T[nbElts];
    } catch (const std::bad_alloc& e) {
      fprintf(stderr, "%s %d : allocation failed : %s", __FILE__, __LINE__,
	      e.what());
    }      
#ifdef DEBUG_MEMORY_PLAINMATRIX
    mem_get();
    fprintf(stderr, "@@PlainMatrix constructor : %.6e %.6e @%x %d\n",  
	    (double)_mem_vrt * 4.0 / (1024.0 *1024.0), 
	    (double)_mem_res * 4.0 / (1024.0 *1024.0), 
	    *_coefs, (sizeof(T) * nbElts));
#endif
  }

  virtual ~PlainMatrix() {
    free();
    if (_coefs_alloc_status) {
      delete _coefs; // same as destructor
      _coefs_alloc_status = false;
    }
  }
#ifdef DEBUG_MEMORY_PLAINMATRIX
  void mem_get()
  {
    int pid = (int)getpid();
    //  cerr << "pid = " << pid << endl;
    char buf[256];
    sprintf(buf, "/proc/%d/statm", pid);
    ifstream fin(buf);
    fin >> _mem_vrt >> _mem_res;
    fin.close();
  }
#endif
  int size() const { return _n; }

  inline T& operator [] (int i) { return (*_coefs)[i]; }
  inline const T& operator [] (int i) const { return (*_coefs)[i]; }
  
  virtual T& operator () (int i, int j) = 0;

  virtual const T& operator () (int i, int j) const = 0;

  T** addrCoefs_pt() { return _coefs; }

  const T** addrCoefs_pt() const { return _coefs; }

  inline T* addrCoefs() { return *_coefs; }

  inline const T* addrCoefs() const { return *_coefs; }

protected:
  void init(int nbElts)
  {
    _isOwner = true;
    if (_n == (-1)) {
      _n = nbElts;
      if (_n > 0) {
	try {
	  *_coefs = new T[nbElts];
	} catch (const std::bad_alloc& e) {
	  fprintf(stderr, "%s %d : allocation failed : %s", __FILE__, __LINE__,
		  e.what());
	}      
#ifdef DEBUG_MEMORY_PLAINMATRIX
	mem_get();
	fprintf(stderr, "@@PlainMatrix init : %.6e %.6e %d\n",  
		(double)_mem_vrt * 4.0 / (1024.0 *1024.0), 
		(double)_mem_res * 4.0 / (1024.0 *1024.0), 
		(sizeof(T) * nbElts));
#endif
      }
    } // _n == (-1)
    else {
#if 0
      if (_n > nbElts) {
	fprintf(stderr,
		"%s %d : PlanMatrix::init(%d) : %d  better deallocate %d\n",
		__FILE__, __LINE__, nbElts, _n, (_n - nbElts));
      }
#endif
      delete [] *_coefs;
      try {
	*_coefs = new T[nbElts];
      } catch (const std::bad_alloc& e) {
	fprintf(stderr, "%s %d : allocation failed : %s", __FILE__, __LINE__,
		e.what());
      }      
      _n = nbElts;
    }
  }

  void init(int nbElts, T* coefs, bool isOwner)
  {
    _isOwner = isOwner;
    _n = nbElts;
    *_coefs = coefs;
  }

  void free() 
  {
    if (_n > 0 && _isOwner) {
#ifdef DEBUG_MEMORY_PLAINMATRIX
      mem_get();
      fprintf(stderr, "@@PlainMatrix free : %.6e %.6e @%x\n",  
	      (double)_mem_vrt * 4.0 / (1024.0 *1024.0), 
	      (double)_mem_res * 4.0 / (1024.0 *1024.0), 
	      *_coefs);
#endif
      delete [] *_coefs;
    }
    _n = (-1);
  }

  T* coefs() 
  {
    return *_coefs;
  }

  const T* coefs() const 
  {
    return *_coefs;
  }

  void ZeroClear()
  {
#if 0    
    const T zero(0.0);
    for (int i = 0; i < _n; i++) {
      (*_coefs)[i] = zero;         // memset()
    }
#else
    memset((void *)*_coefs, 0, sizeof(T) * _n);
#endif
  }
  
  virtual PlainMatrix<T>* clone() const = 0;

  void copy(const PlainMatrix<T> &M) 
  {
    if (_n <= 0) {
      try {
	*_coefs = new T[_n];
      } catch (const std::bad_alloc& e) {
	fprintf(stderr, "%s %d : allocation failed : %s", __FILE__, __LINE__,
		e.what());
      }      
    }
    _isOwner = M._isOwner;
    _n = M._n;
#ifdef DEBUG_MEMORY_PLAINMATRIX
    mem_get();
    fprintf(stderr, "@@PlainMatrix copy : %.6e %.6e @%x %d\n",  
	    (double)_mem_vrt * 4.0 / (1024.0 *1024.0), 
	    (double)_mem_res * 4.0 / (1024.0 *1024.0), 
	    *_coefs, (sizeof(T) * _n));
#endif
#if 0
    for (int i = 0; i < _n; i++) {
      (*_coefs)[i] = (*(M._coefs))[i];  // memcopy
    }
#else
    memcpy((void *)*_coefs, (void *)*(M._coefs), _n * sizeof(T));
#endif
  }

  PlainMatrix(const PlainMatrix& A)
  {
    _isOwner = A._isOwner;
    _n = A._n;
    _coefs = A._coefs;
  }
  PlainMatrix& operator = (const PlainMatrix& A)
  {
    if (this != &A) {
      _isOwner = A._isOwner;
      _n = A._n;
      _coefs = A._coefs;
    }
    return *this;
  }

private:
  bool _isOwner;
  bool _coefs_alloc_status;
  T** _coefs;
  int _n;
#ifdef DEBUG_MEMORY_PLAINMATRIX
  uint64_t _mem_vrt, _mem_res; 
#endif
}; // end class PlainMatrix

#endif
