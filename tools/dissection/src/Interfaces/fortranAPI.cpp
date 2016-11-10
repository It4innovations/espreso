/*! \file   fortranAPI.cpp
    \brief  Fortran style interface 
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

#include <sys/types.h>
#include <stdint.h>
#include <cstdio>
#include <cstdlib>

static int _called = 0;
static int _symbolic = 0;
static int _numeric = 0;
static bool _verbose = true;
#ifdef BLAS_MKL
static int _mkl_num_threads = 1;
#endif
static FILE *_fp;
#ifdef DISSECTION_FORTRAN
#include "Interfaces/fortranAPI.h"
#else
#include "Interfaces/Dissection.hpp"
#endif
#include "Driver/DissectionSolver.hpp"

struct dissection_solver_ptr
{
  int real_or_complex;
  DissectionSolver<double> *rptr;
  DissectionSolver<complex<double>, double> *cptr;
};

DISSECTION_API void DISS_INIT(uint64_t &dslv_, 
			      const int &real_or_complex,
			      const int &nthreads,
			      const int &verbose)
{
  int num_threads;
  dissection_solver_ptr *dslv;

  if (_called == 0) {
    int pid = (int)getpid();
    char fname[256];
    if (verbose > 0) {
      _verbose = true;
    }
    else {
      _verbose = false;
    }
#if 1
    if (_verbose > 0) {
      fprintf(stderr, "pid = %d\n", pid);
      sprintf(fname, "dissection.%04d.log", pid);
      _fp = fopen(fname, "a");
    }
    else {
      _fp = stderr;
    }
#else
    _fp = stderr;
#endif
  }
  if (_verbose > 0) {
    fprintf(_fp, "%s %d : diss_init : called = %d\n", 
	    __FILE__, __LINE__, _called);
  }
  
  _called++;                   // counter for dumping matrix data to debug
#ifdef BLAS_MKL
  if (getenv("MKL_NUM_THREADS")) {
    sscanf(getenv("MKL_NUM_THREADS"), "%d", &_mkl_num_threads);
    if (_verbose > 0) {
      fprintf(_fp,
	      "environmental variable MKL_NUM_THREADS = %d\n",
	      _mkl_num_threads);
    }
  }
  else {
    _mkl_num_threads = mkl_get_max_threads();
  }
  if (_verbose > 0) {
    fprintf(_fp,
	    "MKL_NUM_THREADS = %d\n", _mkl_num_threads);
  }
#endif
  if (nthreads == (-1)) {
    if (getenv("NTHREADS")) {
      sscanf(getenv("NTHREADS"), "%d", &num_threads);
    }
    else {
      num_threads = 1;
    }
  }
  if (nthreads > 0) {
    num_threads = nthreads;
  }
  dslv_ = (uint64_t)new dissection_solver_ptr;
  dslv = (dissection_solver_ptr *)dslv_;
  dslv->real_or_complex = real_or_complex;

  switch(real_or_complex) {
  case DISSECTION_REAL_MATRIX:
    dslv->rptr = new DissectionSolver<double>(num_threads, 
					      (verbose != 0 ? true : false), 
					      _called, _fp);
    break;
  case DISSECTION_COMPLEX_MATRIX:
    dslv->cptr = new DissectionSolver<complex<double>, double>(num_threads, 
							       (verbose != 0 ? true : false), 
							       _called, _fp);
    break;
  default:
    if (_verbose > 0) {
      fprintf(_fp, "%s %d : unknown matrix data type : %d\n", 
	      __FILE__, __LINE__, dslv->real_or_complex);
    }
  }
}

DISSECTION_API void DISS_FREE(uint64_t &dslv_) 
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;

  switch(dslv->real_or_complex) {
  case DISSECTION_REAL_MATRIX:
    delete dslv->rptr;
    break;
  case DISSECTION_COMPLEX_MATRIX:
    delete dslv->cptr;
    break;
  default:
    if (_verbose > 0) {
      fprintf(_fp, "%s %d : unknown matrix data type : %d\n", 
	      __FILE__, __LINE__, dslv->real_or_complex);
    }
  }
  delete dslv;
  dslv_ = (uint64_t)NULL;
  _called--;
  if ((_called == 0) && (_fp != stderr)) {
    fclose(_fp);
  }
#ifdef VECLIB
  unsetenv("VECLIB_MAXIMUM_THREADS");
#endif

}

DISSECTION_API void DISS_NUMERIC_FREE(uint64_t &dslv_)
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;

  switch(dslv->real_or_complex) {
  case DISSECTION_REAL_MATRIX:
    dslv->rptr->NumericFree();
    break;
  case DISSECTION_COMPLEX_MATRIX:
    dslv->cptr->NumericFree();
    break;
  default:
    if (_verbose > 0) {
      fprintf(_fp, "%s %d unknown matrix data type : %d\n", 
	      __FILE__, __LINE__, dslv->real_or_complex);
    }
  }
}
  
DISSECTION_API void DISS_S_FACT(uint64_t &dslv_,
				const int &dim,
				const int *ptRows,
				const int *indCols,
				const int &sym,
				const int &decomposer)
{
  // sym = 1 : symmetric with upper
  //     = 0 : unsymmetric,
  //     = 3 : symmetric with lower
  // decomposer = 0 : SCOTCH 
  //            = 1 : METIS
  //            = 2 : TRIDAIG without nested bisection
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  //  int num_levels;
#if 0
  switch(dslv->real_or_complex) {
  case DISSECTION_REAL_MATRIX:
    fout = dslv->rptr->get_filedescriptor();
    break;
  case DISSECTION_COMPLEX_MATRIX:
    fout = dslv->cptr->get_filedescriptor();
    break;
  default:
    if (_verbose > 0) {
      fprintf(_fp, "%s %d unknown matrix data type : %d\n", 
	      __FILE__, __LINE__, dslv->real_or_complex);
    }
  }
#endif
  switch(dslv->real_or_complex) {
  case DISSECTION_REAL_MATRIX:
    dslv->rptr->SymbolicFact(dim, (int *)ptRows, (int *)indCols,
			     (bool)(sym % 2 == 1), (bool)((sym / 2) == 0),
			     (bool)((sym / 4) == 1),
			     decomposer); // using default parameter
    break;
  case DISSECTION_COMPLEX_MATRIX:
    dslv->cptr->SymbolicFact(dim, (int *)ptRows, (int *)indCols,
			     (bool)(sym % 2 == 1), (bool)((sym / 2) == 0),
			     (bool)((sym / 4) == 1),
			     decomposer); // using default parameter
    break;
  default:
    if (_verbose > 0) {
      fprintf(_fp, "%s %d unknown matrix data type : %d\n", 
	      __FILE__, __LINE__, dslv->real_or_complex);
    }
  }
  if (_verbose > 0) {
    fprintf(_fp, "%s:%d Dissection::SymbolicFact done\n",
	    __FILE__, __LINE__);
  }
  _symbolic++;
}

DISSECTION_API void DISS_N_FACT(uint64_t &dslv_,
				const double *coefs,
				const int &scaling,
				const double &eps_pivot,
				const int &indefinite_flag)
{
  // scaling = 0 : without scaling
  //           1 : 1/sqrt(a_ii) or 1/sqrt(max|a_ij|)
  //           2 : 1/sqrt(a_ii) or Schur complement corresponding to diagonal
  //               kernel_detection_all (for KKT type)
  // eps_pivot = 1.0e-2 : threshold of pivot, ratio of contiguous diagonal
  //                      entries with absolute value
  // indefinite_flag = 1 : indefinite -> kernel_detection_all = false
  // indefinite_flag = 0 : semi-definite -> kernel_detection_all = true

  bool kernel_detection_all = indefinite_flag == 0 ? true : false;
    
  //  FILE *fout;
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;

#ifdef BLAS_MKL
  mkl_set_num_threads(1);
#endif
  #ifdef VECLIB
    setenv("VECLIB_MAXIMUM_THREADS", "1", true);
#endif
  //  dslv->NumericFree(); // for debugging : 20 Nov.2013
  switch(dslv->real_or_complex) {
  case DISSECTION_REAL_MATRIX:
    dslv->rptr->NumericFact(_numeric,
			    (double *)coefs, scaling, 
			    eps_pivot, 
			    kernel_detection_all);
    if (dslv->rptr->getFactorized() == false) {
      fprintf(_fp, "factorization is failed restarted\n");
      dslv->rptr->NumericFree();
      dslv->rptr->NumericFact(_numeric,
			      (double *)coefs, scaling, 
			      eps_pivot, 
			      kernel_detection_all);
      if (dslv->rptr->getFactorized() == false) {
	fprintf(_fp, "factorization is failed again : matrix is dumped\n");
	dslv->rptr->SaveMMMatrix(_called, coefs);
	exit(-1);
      }
    }
    break;
  case DISSECTION_COMPLEX_MATRIX:
    dslv->cptr->NumericFact(_numeric,
			    (complex<double> *)coefs, scaling, 
			    eps_pivot, 
			    kernel_detection_all);
    break;
  }

#ifdef BLAS_MKL
  mkl_set_num_threads(_mkl_num_threads);
#endif
  if (_verbose > 0) {
    fprintf(_fp,
	    "%s %d : Dissection::NumericFact done : %d\n",
	    __FILE__, __LINE__, _numeric);
  }
  _numeric++;
}

DISSECTION_API void DISS_GET_KERN_DIM(uint64_t &dslv_, 
				      int *n0)
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;

  switch(dslv->real_or_complex) {
  case DISSECTION_REAL_MATRIX:
    *n0 = dslv->rptr->kern_dimension();
    break;
  case DISSECTION_COMPLEX_MATRIX:
    if (_verbose > 0) {
      fprintf(_fp, 
	      "%s %d diss_get_kern_dim() for complex is not yet implemented\n",
	      __FILE__, __LINE__);
    }
    *n0 = 0;
    break;
  default:
    if (_verbose > 0) {
      fprintf(_fp, "%s %d unknown matrix data type : %d\n", 
	      __FILE__, __LINE__, dslv->real_or_complex);
    }
  }
}

DISSECTION_API void DISS_GET_NULLPIVOTS(uint64_t &dslv_, 
					int *pivots)
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  switch(dslv->real_or_complex) {
  case DISSECTION_REAL_MATRIX:
    dslv->rptr->GetNullPivotIndices(pivots);
    break;
  case DISSECTION_COMPLEX_MATRIX:
    if (_verbose > 0) {
      fprintf(_fp, 
	      "%s %d diss_get_nullpivots() for complex is not yet implemented\n",
	      __FILE__, __LINE__);
    }
    break;
  default:
    if (_verbose > 0) {
      fprintf(_fp, "%s %d unknown matrix data type : %d\n", 
	      __FILE__, __LINE__, dslv->real_or_complex);
    }
  }
}

DISSECTION_API void DISS_GET_KERN_VECS(uint64_t &dslv_, 
				       double *vec)
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  switch(dslv->real_or_complex) {
  case DISSECTION_REAL_MATRIX:
    dslv->rptr->GetKernelVectors(vec);
    break;
  case DISSECTION_COMPLEX_MATRIX:
    if (_verbose > 0) {
      fprintf(_fp, 
	      "%s %d diss_get_kern_vecs() for complex is not yet implemented\n",
	      __FILE__, __LINE__);
    }
    break;
  default:
    if (_verbose > 0) {
      fprintf(_fp, "%s %d unknown matrix data type : %d\n", 
	      __FILE__, __LINE__, dslv->real_or_complex);
    }
  }

}

DISSECTION_API void DISS_GET_KERNT_VECS(uint64_t &dslv_, 
				       double *vec)
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  switch(dslv->real_or_complex) {
  case DISSECTION_REAL_MATRIX:
    dslv->rptr->GetTransKernelVectors(vec);
    break;
  case DISSECTION_COMPLEX_MATRIX:
    if (_verbose > 0) {
      fprintf(_fp, 
	      "%s %d diss_get_kern_vecs() for complex is not yet implemented\n",
	      __FILE__, __LINE__);
    }
    break;
  default:
    if (_verbose > 0) {
      fprintf(_fp, "%s %d unknown matrix data type : %d\n", 
	      __FILE__, __LINE__, dslv->real_or_complex);
    }
  }

}

DISSECTION_API void DISS_PROJECT(uint64_t &dslv_, 
				 double *x)
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  int n0;
  switch(dslv->real_or_complex) {
  case DISSECTION_REAL_MATRIX:
    n0 = dslv->rptr->kern_dimension();
    if (n0 > 0) {
      dslv->rptr->ProjectionImageSingle(x);
    }
    break;
  case DISSECTION_COMPLEX_MATRIX:
    if (_verbose > 0) {
      fprintf(_fp, 
	      "%s %d diss_project() for complex is not yet implemented\n",
		__FILE__, __LINE__);
    }
    break;
  default:
    if (_verbose > 0) {
      fprintf(_fp, "%s %d unknown matrix data type : %d\n", 
	      __FILE__, __LINE__, dslv->real_or_complex);
    }
  }
}

DISSECTION_API void DISS_SOLVE_1(uint64_t &dslv_, 
				 double *x, 
				 const int &projection,
				 const int &trans)
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  switch(dslv->real_or_complex) {
  case DISSECTION_REAL_MATRIX:
    dslv->rptr->SolveSingle(x, (bool)(projection == 1), (bool)(trans == 1),
			    true); 
                                                     // isTrans
    break;
  case DISSECTION_COMPLEX_MATRIX:    
    fprintf(dslv->cptr->get_filedescriptor(), 
	  "Dissection::SolveSingle : %p\n", dslv->cptr);
    dslv->cptr->SolveSingle((complex<double> *)x, (bool)(projection == 1), 
			    (bool)(trans == 1), true); 
                                                     // isTrans
    break;
  default:
    if (_verbose > 0) {
      fprintf(_fp, "%s %d unknown matrix data type : %d\n", 
	      __FILE__, __LINE__, dslv->real_or_complex);
    }
  }
}

DISSECTION_API void DISS_SOLVE_N(uint64_t &dslv_, 
				 double *x, 
				 const int &nrhs, const int &projection)
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  switch(dslv->real_or_complex) {
  case DISSECTION_REAL_MATRIX:
    dslv->rptr->SolveMulti(x, nrhs, (bool)(projection == 1), false, true);
    break;
  case DISSECTION_COMPLEX_MATRIX:
    dslv->cptr->SolveMulti((complex<double> *)x, nrhs, (bool)(projection == 1),
			   false, true);
    break;
  default:
    if (_verbose > 0) {
      fprintf(_fp, "%s %d unknown matrix data type : %d\n", 
	   __FILE__, __LINE__, dslv->real_or_complex);
    }
  }
}

DISSECTION_API void DISS_MATRIX_PRODUCT(uint64_t &dslv_,
					const double* x, double* y)
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  switch(dslv->real_or_complex) {
  case DISSECTION_REAL_MATRIX:
    dslv->rptr->SpMV(x, y);
    break;
  case DISSECTION_COMPLEX_MATRIX:
    dslv->cptr->SpMV((complex<double> *)x, (complex<double> *)y);    
    break;
  default:
    if (_verbose > 0) {
      fprintf(_fp, "%s %d unknown matrix data type : %d\n", 
	   __FILE__, __LINE__, dslv->real_or_complex);
    }
  }
}
