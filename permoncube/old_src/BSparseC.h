
#ifndef BSPARSEC_H_
#define BSPARSEC_H_

#include "utility.h"
#include "Fem.h"
#include "Data.h"
#include "Solid45EqNumbGlobal.h"
#include "StiffnessLocal.h"
//#include "Cluster_g.h"
#include <iostream>
#include "mpi.h"
#include <vector>
#include <algorithm>

#ifdef HFETI_SOLVER
	#include <cilk/cilk.h>
	#include <cilk/cilk_api.h>
	
	#include <tbb/mutex.h>
	#include "tbb/parallel_sort.h"
	#include "tbb/tbb.h"
	#include <math.h>
#endif

class CFem;
class CData;

class CBSparseC {

	public:
    MPI_Comm comm;
    int n_col;

		int *i_eq; 						
		int *j_eq; 						
		double *v_eq;     		
    int *Bi_coo;      		
    int *Bj_coo;      		
    double *Bv_coo;   		
    int n_row_eq;
	 
		double 	*v_bg;    		
		int 		*i_bg;    		
		int 		*j_bg;    		
		int 		n_row_bg;
		int 		nnz_bg;

		double 	*v_bd;    		
		int 		*i_bd;    		
		int 		*j_bd;    		
		int 		n_row_bd;
		int 		nnz_bd;

		double 	*v_bc;        
		int 		*i_bc;        
		int 		*i_bc_crs;    
		int 		*j_bc;        
		int 		n_row_bc;
		int 		nnz_bc;

    int nnz_eq;

    int *multiplicityPD;
    int *multiplicityDD;
    int *i_multiplicityDD;

		vector<int>    BI;
		vector<int>    BJ;
		vector<double> BV;

		vector<int>    BI_l;
		vector<int>    BJ_l;
		vector<double> BV_l;


		vector < vector < int > > lambda_map_sub; 
		vector < int > neigh_domains;

	public:
		CBSparseC(MPI_Comm comm);
		virtual ~CBSparseC();

	public:
//		void createB1(CClust_g & clust);//CFem *fem, CDomain *domainG, int nPrimalDOFs, int *indExtDOFs);
		void createB1(CDomain * domainG, 
                  vector<CFem*> &fem,vector<CData*> &data,
                  int *indExterDOFsOnClust,
                  int n_exterDOFsClust,
                  int *neighbClst,
                  int n_neighbClst,
                  bool flag_redund_lagr_mult);
		static int coo2crs(int m, int nnz, int *row_ind, int *row_ptr);
};

#endif /* BSPARSE_H_ */ 

