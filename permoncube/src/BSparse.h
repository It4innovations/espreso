
#ifndef BSPARSE_H_
#define BSPARSE_H_

#include "utility.h"
#include "Fem.h"
#include "Solid45EqNumbGlobal.h"
#include "StiffnessLocal.h"
#include <iostream>
#include "mpi.h"
#include <vector>
#include <algorithm>

class CFem;

class CBSparse {

	public:
    MPI_Comm comm;
    int n_col;
    int n_bd;

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
    int     n_col_bg;
		int 		nnz_bg;

		double 	*v_bd;    		
		int 		*i_bd;    		
		int 		*j_bd;    		
		int 		n_row_bd;
    int     n_col_bd;
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

		vector<int>    B0_I;
		vector<int>    B0_J;
		vector<double> B0_V;
		int B0_rows; 
		int B0_cols; 
		int B0_nnz;  


		vector<int>    BI_l;
		vector<int>    BJ_l;
		vector<double> BV_l;

		vector<int>    BI_dir;
		vector<int>    BJ_dir;
		vector<double> BV_dir;

		vector<int>    BI_full;
		vector<int>    BJ_full;
		vector<double> BV_full;
		
		int            B_full_cols; 
		int            B_full_rows; 
		int            B_full_nnz; 

		vector < int > lambda_map_sub; 
		vector < int > neigh_domains;
		vector < int > neigh_clusters;

	public:
		CBSparse(MPI_Comm comm);
		virtual ~CBSparse();

	public:
		void createB(CFem *fem, CDomain *domainG, int nPrimalDOFs, int *indExtDOFs);
		static int coo2crs(int m, int nnz, int *row_ind, int *row_ptr);
};

#endif /* BSPARSE_H_ */ 

