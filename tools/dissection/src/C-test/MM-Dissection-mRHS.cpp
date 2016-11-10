/*! \file   MM-Dissection-mRHS.cpp
    \brief  test rouinte of dissection solver reading Matrix Market format
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
// along with Disection.  If not, see <http://www.gnu.org/licenses/>.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>

#include "Driver/DissectionSolver.hpp"
#include <mkl_service.h>
#include <pthread.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <complex>
using namespace std;

static int _stat = (-1);

void *thread_child(void *arg)
{
  char buf[256];
  int pid = (int)arg;
  unsigned int mem_tmp, mem_min, mem_max;
  double avg_mem;
  avg_mem = 0.0;
  mem_min = 1U << 32 - 1;
  mem_max = 0U;
  int stat0, stat1;
  stat0 = _stat;
  unsigned int count = 0U;
//  fprintf(stderr, "thread_child forked\n");
  while(_stat != 0) {
     if (_stat == 1) {
       sprintf(buf, "/proc/%d/statm", pid);
       ifstream fin(buf);
       fin >> mem_tmp;
       fin.close();
       if (mem_tmp > mem_max) {
	 mem_max = mem_tmp;
       }
       if (mem_tmp < mem_min) {
	 mem_min = mem_tmp;
       }
       avg_mem += (double)mem_tmp;
       count++;
     }
     stat1 = _stat;
     if ((stat1 == (-1)) && (stat0 == 1)) {
       fprintf(stderr, 
	       "used memory :min: %14.8e  max: %14.8e avg: %14.8e count: %d\n", 
	       (double)mem_min * 4.0 / (1024.0 * 1024.0),
	       (double)mem_max * 4.0 / (1024.0 * 1024.0),
	       (avg_mem / (double)count) * 4.0 / (1024.0 * 1024.0),
	       count);
       count = 0U;
       avg_mem = 0.0;
       mem_min = 1U << 32 - 1;
       mem_max = 0U;
     }
     stat0 = stat1;
     usleep(1000);
  }
// fprintf(stderr, "thread_child join\n count = %ld\n", count);
  pthread_exit(arg);

  return (void *)NULL;
}

template<typename T>
void generate_CSR(std::list<int>* &ind_cols_tmp, std::list<T>* &val_tmp, 
		  int nrow, int nnz,
		  int *irow, int *jcol, T* val)
{
  ind_cols_tmp = new std::list<int>[nrow];
  val_tmp = new std::list<T>[nrow];
  for (int i = 0; i < nnz; i++) {
    const int ii = irow[i];
    const int jj = jcol[i];
    if (ind_cols_tmp[ii].empty()) {
      ind_cols_tmp[ii].push_back(jj);
      val_tmp[ii].push_back(val[i]);
    }
    else {
      if (ind_cols_tmp[ii].back() < jj) {
	ind_cols_tmp[ii].push_back(jj);
	val_tmp[ii].push_back(val[i]);
      }
      else {
	typename std::list<T>::iterator iv = val_tmp[ii].begin();
	std::list<int>::iterator it = ind_cols_tmp[ii].begin();
	for ( ; it != ind_cols_tmp[ii].end(); ++it, ++iv) {
	  if (*it == jj) {
	    fprintf(stderr, "already exits? (%d %d)\n", ii, jj);
	      break;
	  }
	  if (*it > jj) {
	    ind_cols_tmp[ii].insert(it, jj);
	    val_tmp[ii].insert(iv, val[i]);
	    break;
	  }
	}
      }
    }
  }
}

template<typename T>
void copy_CSR(int *indcols, int *ptrows, T* coefs, int nrow, 
	      bool upper_flag, bool isSym,
	      std::list<int>* ind_cols_tmp, std::list<T>* val_tmp)
{
  const T zero(0.0);
  ptrows[0] = 0;
  for (int i = 0; i < nrow; i++) {
    int k;
    int itmp = ind_cols_tmp[i].size();
    if (upper_flag) {
      if (ind_cols_tmp[i].front() == i) {
	ptrows[i + 1] = ptrows[i] + itmp;
	k = ptrows[i];
      }
      else {
	fprintf(stderr, "zero is added to diagonal : %d\n", i);
	ptrows[i + 1] = ptrows[i] + itmp + 1;
	indcols[ptrows[i]] = i;
	coefs[ptrows[i]] = zero;
	k = ptrows[i] + 1;
      }
    }
    else {
      k = ptrows[i];
      if (ind_cols_tmp[i].back() == i || (!isSym)) {
	ptrows[i + 1] = ptrows[i] + itmp;
      }
      else {
	fprintf(stderr, "zero is added to diagonal : %d\n", i);
	ptrows[i + 1] = ptrows[i] + itmp + 1;
	indcols[ptrows[i + 1] - 1] = i;
	coefs[ptrows[i + 1] - 1] = zero;
      }
    }
    std::list<int>::iterator it = ind_cols_tmp[i].begin();
    typename std::list<T>::iterator iv = val_tmp[i].begin();
    for ( ; it != ind_cols_tmp[i].end(); ++it, ++iv, k++) {
      indcols[k] = *it;
      coefs[k] = *iv;
    }
  } // loop : i
}


int main(int argc, char **argv)
{
  int n, itmp;
  char fname[256];
  char buf[1024];
  int nrow, nnz, flag;
  int *ptrows, *indcols, *indvals;
  int *irow, *jcol;
  double *val, *coefs;
  complex<double> *valc, *ccoefs;
  int decomposer;
  int num_threads;
  double eps_pivot;
  int numlevels;
  int minNodes = 64;
  std::list<int>* ind_cols_tmp;
  std::list<double>* val_tmp;
  std::list<complex<double> >* val_tmpc;
  FILE *fp;
  bool isSym, isComplex;
  bool upper_flag = true;

  if (argc < 6) {
    fprintf(stderr, 
	    "MM-dissection [data file] [decomposer] [num_threads] [eps_pivot] [num_levels]\n");
    exit(-1);
  }    
  strcpy(fname, argv[1]);
  decomposer = atoi(argv[2]);
  num_threads = atoi(argv[3]);
  eps_pivot = atof(argv[4]);
  numlevels = atof(argv[5]);
  if (argc >= 7) {
    minNodes = atoi(argv[6]);
  }
  if (argc >= 8) {
    upper_flag = (atoi(argv[7]) == 1);
  }

  // read from the file
  if ((fp = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "fail to open %s\n", fname);
  }
  fgets(buf, 256, fp);
  //
  if (strstr(buf, "symmetric") != NULL) {
   isSym = true;
  }
  else {
    isSym = false;
    upper_flag = false;
  }
  if (strstr(buf, "complex") != NULL) {
   isComplex = true;
  }
  else {
    isComplex = false;
  }

  fprintf(stderr, "symmetric = %s\n", isSym ? "true " : "false");
  while (1) {
    fgets(buf, 256, fp);
    if (buf[0] != '%') {
      sscanf(buf, "%d %d %d", &nrow, &itmp, &nnz);
      break;
    }
  }
  irow = new int[nnz];
  jcol = new int[nnz];

  if (isComplex) {
    valc = new complex<double>[nnz];
    if (upper_flag) {
      for (int i = 0; i < nnz; i++) {
	fscanf(fp, "%d %d %lf %lf", &jcol[i], &irow[i], 
	       &(valc[i].real()), &(valc[i].imag()));
	irow[i]--;
	jcol[i]--;
	if (isSym && irow[i] > jcol[i]) {
	  fprintf(stderr, "exchanged : %d > %d\n", irow[i], jcol[i]);
	  itmp = irow[i];
	  irow[i] = jcol[i];
	  jcol[i] = itmp;
	}
      }
    }
    else {
      for (int i = 0; i < nnz; i++) {
	fscanf(fp, "%d %d %lf %lf", &irow[i], &jcol[i], 
	       &(valc[i].real()), &(valc[i].imag())); 
	irow[i]--;
	jcol[i]--;
	if (isSym && irow[i] < jcol[i]) {
	  fprintf(stderr, "exchanged : %d > %d\n", irow[i], jcol[i]);
	  itmp = irow[i];
	  irow[i] = jcol[i];
	  jcol[i] = itmp;
	}
      }
    }
  }
  else {
    val = new double[nnz];
    if (upper_flag) {
      for (int i = 0; i < nnz; i++) {
	fscanf(fp, "%d %d %lf", &jcol[i], &irow[i], &val[i]); // read lower
	irow[i]--;
	jcol[i]--;
	if (isSym && irow[i] > jcol[i]) {
	  fprintf(stderr, "exchanged : %d > %d\n", irow[i], jcol[i]);
	  itmp = irow[i];
	  irow[i] = jcol[i];
	  jcol[i] = itmp;
	}
    }
    }
    else {
      for (int i = 0; i < nnz; i++) {
	fscanf(fp, "%d %d %lf", &irow[i], &jcol[i], &val[i]); // read lower
	irow[i]--;
	jcol[i]--;
	if (isSym && irow[i] < jcol[i]) {
	  fprintf(stderr, "exchanged : %d > %d\n", irow[i], jcol[i]);
	  itmp = irow[i];
	  irow[i] = jcol[i];
	  jcol[i] = itmp;
	}
      }
    }
  }
  fclose (fp);

  ind_cols_tmp = new std::list<int>[nrow];
  if (isComplex) {
    val_tmpc = new std::list<complex<double> >[nrow];
    generate_CSR<complex<double> >(ind_cols_tmp, val_tmpc, 
				   nrow, nnz,
				   irow, jcol, valc);
  }
  else {
    val_tmp = new std::list<double>[nrow];
    generate_CSR<double>(ind_cols_tmp, val_tmp, 
			 nrow, nnz,
			 irow, jcol, val);
  }
  delete [] irow;
  delete [] jcol;
  if (isComplex) {
    delete [] valc;
  }
  else {
    delete [] val;
  }

  if (upper_flag) {
    for (int i = 0; i < nrow; i++) {
      if (ind_cols_tmp[i].front() != i) {
	nnz++;
      }
    }
  }
  else {
    for (int i = 0; i < nrow; i++) {
      if (ind_cols_tmp[i].back() != i) {
	nnz++;
      }
    }
  }
  ptrows = new int[nrow + 1];
  indcols = new int[nnz];
  indvals = new int[nnz];
  if (isComplex) {
    ccoefs = new complex<double>[nnz];
    copy_CSR<complex<double> >(indcols, ptrows, ccoefs, 
			       nrow, upper_flag, isSym, 
			       ind_cols_tmp, val_tmpc);
  }
  else {
    coefs = new double[nnz];
    copy_CSR<double>(indcols, ptrows, coefs, 
		     nrow, upper_flag, isSym,
		     ind_cols_tmp, val_tmp);
  }

#if 0
  if ((fp = fopen("debug.matrix.data", "w")) != NULL) {
    for (int i = 0; i < nrow; i++) {
      fprintf(fp, "%d : %d :: ", i, (ptrows[i + 1] - ptrows[i]));
      for (int k = ptrows[i]; k < ptrows[i + 1]; k++) {
	fprintf(fp, "%d ", indcols[k]);
      }
      fprintf(fp, "\n");
    }
  }
  fclose(fp);
#endif
  int pid = (int)getpid();
#if 1
  fprintf(stderr, "pid = %d\n", pid);
  sprintf(fname, "dissection.%04d.log", pid);
  fp = fopen(fname, "a");
#else
  fp = stderr;
#endif

  void* results;
  pthread_attr_t th_attr;
  pthread_t thread;
  pthread_attr_init(&th_attr);
  pthread_attr_setdetachstate(&th_attr, PTHREAD_CREATE_JOINABLE);       
  int pthid = pthread_create(&thread, &th_attr, 
			     &thread_child,
			     (void *)pid);
  if (pthid != 0) {
    cout << "bad thread creation ? " << pid << endl;
    exit(0);
  }
  //  DissectionSolver *dslv = new DissectionSolver(num_threads > 1 ? 2 : 1, true, 0, fp);

    clock_t t0_cpu, t1_cpu, t2_cpu, t3_cpu;
    clock_t *tbegin_cpu, *tend_cpu;
    elapsed_t t0_elapsed, t1_elapsed, t2_elapsed, t3_elapsed;
    elapsed_t *tbegin_elapsed, *tend_elapsed;

    int num_levels;
    const int scaling = 1;
    const int pivot_strategy = 0;
    const int dim_augkern = 4;
    int called = 0;
    
    mkl_set_num_threads(1);
    t0_cpu = clock();
    get_realtime(&t0_elapsed);
    if (numlevels <= 0) {
      num_levels = (int)log2((double)nrow / (double)minNodes);
    }
    else {
      num_levels = numlevels;
    }

  if (isComplex) {
    DissectionSolver<complex<double> > *dslv = 
      new DissectionSolver<complex<double> >(num_threads, true, 0, fp);

    dslv->SymbolicFact(nrow, (int *)ptrows, (int *)indcols,
		       isSym,
		       upper_flag,
		       decomposer, num_levels, minNodes); 
    //                  sym, upper
    
    t1_cpu = clock();
    get_realtime(&t1_elapsed);
    
    t2_cpu = clock();
    get_realtime(&t2_elapsed);
    _stat = 1;
    fprintf(stderr, "%s %d : NumericFact()\n", __FILE__, __LINE__);
    dslv->NumericFact((complex<double> *)ccoefs, scaling, 
		      pivot_strategy, eps_pivot, 
		      dim_augkern);
    fprintf(stderr, "%s %d : NumericFact()\n", __FILE__, __LINE__);
    _stat = (-1);
    t3_cpu = clock();
    get_realtime(&t3_elapsed);
    usleep(10000); // sleep
    
    _stat = 1;
    int nrhs = 40;
    complex<double> *x = new complex<double>[nrow];
    complex<double> *y = new complex<double>[nrow * nrhs];
    complex<double> *z = new complex<double>[nrow * nrhs];
    int ntrial = 24;
    tbegin_cpu = new clock_t[ntrial];
    tend_cpu =  new clock_t[ntrial];
    tbegin_elapsed = new elapsed_t[ntrial];
    tend_elapsed = new elapsed_t[ntrial];
    fprintf(stderr, "%s %d : RHS allocated\n", __FILE__, __LINE__);
    _stat = (-1);
  
    for (int k = 0; k < nrhs; k++) {
      const int kshft = k * nrow;
      for (int i = 0; i < nrow; i++) {
	y[kshft + i] = complex<double>((double)((i + k) % 11), 0.0);
      }
      dslv->SpMV(y + kshft, x);

      dslv->SpMV(x, y + kshft);
      for (int i = 0; i < nrow * nrhs; i++) {
	z[i] = y[i];
      }
    }
    _stat = 1;
    for (int k = 1; k <= 20; k++) {
      for (int i = 0; i < nrow * nrhs; i++) {
	y[i] = z[i];
      }
      tbegin_cpu[k - 1] = clock();
      get_realtime(&tbegin_elapsed[k - 1]);  
      if (k == 1) {
	fprintf(stderr, "%s %d : SolveSingle()\n", __FILE__, __LINE__);
	dslv->SolveSingle(y, false, false);
	fprintf(stderr, "%s %d : SolveSingle()\n", __FILE__, __LINE__);
      }
      else {
	const int nnrhs = k;
	fprintf(stderr, "%s %d : SolveMulti() %d\n", __FILE__, __LINE__, nnrhs);
	dslv->SolveMulti(y, nnrhs, false, false); 
	fprintf(stderr, "%s %d : SolveMulti() %d\n", __FILE__, __LINE__, nnrhs);
      }
      tend_cpu[k - 1] = clock();
      get_realtime(&tend_elapsed[k - 1]);
    }
    for (int k = 5; k < 9; k++) {
      for (int i = 0; i < nrow * nrhs; i++) {
	y[i] = z[i];
      }
      tbegin_cpu[k + 15] = clock();
      get_realtime(&tbegin_elapsed[k + 15]);  
      const int nnrhs = k * 5;
      fprintf(stderr, "%s %d : SolveMulti() %d\n", __FILE__, __LINE__, nnrhs);
      dslv->SolveMulti(y, nnrhs, false, false); 
      fprintf(stderr, "%s %d : SolveMulti() %d\n", __FILE__, __LINE__, nnrhs);
      tend_cpu[k + 15] = clock();
      get_realtime(&tend_elapsed[k + 15]);
    }
    _stat = (-1);
#if 1
    for (int k = 0; k < nrhs; k++) {
      const int kshft = k * nrow;
      for (int i = 0; i < nrow; i++) {
	z[i] = complex<double>((double)((i + k) % 11), 0.0);
      }
      dslv->SpMV(z, x);
      dslv->SpMV(x, z);
      
      double norm0, norm1;
      norm0 = 0.0;
      norm1 = 0.0;
      for (int i = 0; i < nrow; i++) {
	norm0 += x[i].real() * x[i].real() + x[i].imag() * x[i].imag();
	complex<double> ztmp;
	ztmp = y[i + kshft] - x[i];
	norm1 += ztmp.real() * ztmp.real() + ztmp.imag() * ztmp.imag();
      }
      fprintf(fp, "## %2d : error    : %18.7e = %18.7e / %18.7e\n", k,
	      sqrt(norm1 / norm0), sqrt(norm1), sqrt(norm0));
      dslv->SpMV(y + kshft, x);
      
      norm0 = 0.0;
      norm1 = 0.0;
      for (int i = 0; i < nrow; i++) {
	norm0 += z[i].real() * z[i].real() + z[i].imag() * z[i].imag();
	complex<double> ztmp;
	ztmp = z[i] - x[i];
	norm1 += ztmp.real() * ztmp.real() + ztmp.imag() * ztmp.imag();
      }
      fprintf(fp, "## %2d : residual : %18.7e = %18.7e / %18.7e\n", k, 
	      sqrt(norm1 / norm0), sqrt(norm1), sqrt(norm0));
    }
#endif
    delete dslv;
    delete [] ccoefs;
    delete [] x;
    delete [] y;
    delete [] z;

  }
  else {
    DissectionSolver<double> *dslv = 
      new DissectionSolver<double>(num_threads, true, 0, fp);
    
    dslv->SymbolicFact(nrow, (int *)ptrows, (int *)indcols,
		       isSym,
		       upper_flag,
		       decomposer, num_levels, minNodes); 
    //                  sym, upper
    
    t1_cpu = clock();
    get_realtime(&t1_elapsed);
    
    t2_cpu = clock();
    get_realtime(&t2_elapsed);
    _stat = 1;
    fprintf(stderr, "%s %d : NumericFact()\n", __FILE__, __LINE__);
    dslv->NumericFact((double *)coefs, scaling, 
		      pivot_strategy, eps_pivot, 
		      dim_augkern);
    fprintf(stderr, "%s %d : NumericFact()\n", __FILE__, __LINE__);
    _stat = (-1);
    t3_cpu = clock();
    get_realtime(&t3_elapsed);
    usleep(10000); // sleep
    int n0;
    n0 = dslv->kern_dimension();
    fprintf(fp, "## kernel dimension = %d\n", n0);
    
    _stat = 1;
    int nrhs = 40;
    double *x = new double[nrow];
    double *y = new double[nrow * nrhs];
    double *z = new double[nrow * nrhs];
    int ntrial = 24;
    tbegin_cpu = new clock_t[ntrial];
    tend_cpu =  new clock_t[ntrial];
    tbegin_elapsed = new elapsed_t[ntrial];
    tend_elapsed = new elapsed_t[ntrial];
    fprintf(stderr, "%s %d : RHS allocated\n", __FILE__, __LINE__);
    _stat = (-1);
    
    if (!isSym && (n0 > 0)) {
      _stat = 1;
      fprintf(stderr, "%s %d : ComputeTransposedKernels()\n", __FILE__, __LINE__);
      dslv->ComputeTransposedKernels(true);
      fprintf(stderr, "%s %d : ComputeTransposedKernels()\n", __FILE__, __LINE__);
      _stat = (-1);
    }
    
    for (int k = 0; k < nrhs; k++) {
      const int kshft = k * nrow;
      for (int i = 0; i < nrow; i++) {
	y[kshft + i] = (double)((i + k) % 11);
      }
      dslv->SpMV(y + kshft, x);
      if (n0 > 0) {
	dslv->ProjectionImageSingle(x);
      }
      dslv->SpMV(x, y + kshft);
      for (int i = 0; i < nrow * nrhs; i++) {
	z[i] = y[i];
      }
    }
    _stat = 1;
    for (int k = 1; k <= 20; k++) {
      for (int i = 0; i < nrow * nrhs; i++) {
	y[i] = z[i];
      }
      tbegin_cpu[k - 1] = clock();
      get_realtime(&tbegin_elapsed[k - 1]);  
      if (k == 1) {
	fprintf(stderr, "%s %d : SolveSingle()\n", __FILE__, __LINE__);
	dslv->SolveSingle(y, false, false);
	fprintf(stderr, "%s %d : SolveSingle()\n", __FILE__, __LINE__);
      }
      else {
	const int nnrhs = k;
	fprintf(stderr, "%s %d : SolveMulti() %d\n", __FILE__, __LINE__, nnrhs);
	dslv->SolveMulti(y, nnrhs, false, false); 
	fprintf(stderr, "%s %d : SolveMulti() %d\n", __FILE__, __LINE__, nnrhs);
      }
      if (n0 > 0) {
	dslv->ProjectionImageSingle(y);
      }
      tend_cpu[k - 1] = clock();
      get_realtime(&tend_elapsed[k - 1]);
    }
    for (int k = 5; k < 9; k++) {
      for (int i = 0; i < nrow * nrhs; i++) {
	y[i] = z[i];
      }
      tbegin_cpu[k + 15] = clock();
      get_realtime(&tbegin_elapsed[k + 15]);  
      const int nnrhs = k * 5;
      fprintf(stderr, "%s %d : SolveMulti() %d\n", __FILE__, __LINE__, nnrhs);
      dslv->SolveMulti(y, nnrhs, false, false); 
      fprintf(stderr, "%s %d : SolveMulti() %d\n", __FILE__, __LINE__, nnrhs);
      if (n0 > 0) {
	dslv->ProjectionImageSingle(y);
      }
      tend_cpu[k + 15] = clock();
      get_realtime(&tend_elapsed[k + 15]);

    }
    _stat = (-1);
#if 1
    for (int k = 0; k < nrhs; k++) {
      const int kshft = k * nrow;
      for (int i = 0; i < nrow; i++) {
	z[i] = (double)((i + k) % 11);
      }
      dslv->SpMV(z, x);
      dslv->SpMV(x, z);
      
      double norm0, norm1;
      norm0 = 0.0;
      norm1 = 0.0;
      for (int i = 0; i < nrow; i++) {
	norm0 += x[i] * x[i];
	double ztmp;
	ztmp = y[i + kshft] - x[i];
	norm1 += ztmp * ztmp;
      }
      fprintf(fp, "## %2d : error    : %18.7e = %18.7e / %18.7e\n", k,
	      sqrt(norm1 / norm0), sqrt(norm1), sqrt(norm0));
      dslv->SpMV(y + kshft, x);
      
      norm0 = 0.0;
      norm1 = 0.0;
      for (int i = 0; i < nrow; i++) {
	norm0 += z[i] * z[i];
	double ztmp;
	ztmp = z[i] - x[i];
	norm1 += ztmp * ztmp;
      }
      fprintf(fp, "## %2d : residual : %18.7e = %18.7e / %18.7e\n", k, 
	      sqrt(norm1 / norm0), sqrt(norm1), sqrt(norm0));
    }
#endif

  delete dslv;
  delete [] coefs;
  delete [] x;
  delete [] y;
  delete [] z;

  }
  
  _stat = 0;
  pthread_attr_destroy(&th_attr);
  pthid = pthread_join(thread, &results);  
  if (pthid != 0) {
    cout << "bad thread join ? " << pthid << endl;
    exit(0);
  }

  fprintf(fp, "## symbolic fact    : cpu time = %.4e elapsed time = %.4e\n", 
	 (double)(t1_cpu - t0_cpu) / (double)CLOCKS_PER_SEC,
	 convert_time(t1_elapsed, t0_elapsed));
 
  fprintf(fp, "## numeric fact     : cpu time = %.4e elapsed time = %.4e\n", 
	 (double)(t3_cpu - t2_cpu) / (double)CLOCKS_PER_SEC,
	 convert_time(t3_elapsed, t2_elapsed));

  for (int k = 1; k <= 20; k++) {
    if (k == 1) {
      fprintf(fp,
	      "## solve single RHS : cpu time = %.4e elapsed time = %.4e\n", 
	      (double)(tend_cpu[0] - tbegin_cpu[0]) / (double)CLOCKS_PER_SEC,
	      convert_time(tend_elapsed[0], tbegin_elapsed[0]));
    }
    else {
      const int nnrhs = k;
      fprintf(fp,
	      "## solve %d mlt-RHS : cpu time = %.4e elapsed time = %.4e\n", 
	      nnrhs,
	      (double)(tend_cpu[k - 1] - tbegin_cpu[k - 1]) / 
	      (double)CLOCKS_PER_SEC,
	      convert_time(tend_elapsed[k - 1], tbegin_elapsed[k - 1]));
    }
  }
  for (int k = 5; k < 9; k++) {
    const int nnrhs = k * 5;
    fprintf(fp,
	    "## solve %d mlt-RHS : cpu time = %.4e elapsed time = %.4e\n", 
	    nnrhs,
	    (double)(tend_cpu[k + 15] - tbegin_cpu[k + 15]) / 
	    (double)CLOCKS_PER_SEC,
	    convert_time(tend_elapsed[k + 15], tbegin_elapsed[k + 15]));
    }

    
#if 0
  double *scalediag;
  scalediag = new double[nrow];
  dslv->GetMatrixScaling(scalediag);
  for (int i = 0; i < nrow; i++) {
    for (int k = ptrows[i]; k < ptrows[i + 1]; k++) {
      if (indcols[k] == i) {
	fprintf(fp, "%i : %g : %g %g\n", i, scalediag[i], 
		1.0 / (scalediag[i] * scalediag[i]), coefs[k]);
	break;
      }
    }
  }
  delete [] scalediag;
#endif

  delete [] ptrows;
  delete [] indcols;
  fclose(fp);

}
