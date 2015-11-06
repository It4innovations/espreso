#ifdef WIN32	 
#include "stdafx.h"
#endif

#include <omp.h>
#include "mpi.h"
#include "mkl.h"

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include <map>

using std::vector;
using std::cout;
using std::map; 
using std::make_pair; 

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

#include "SparseMatrix.h"
#include "FEM_Assembler.h"
#include "SparseSolver.h"
#include "TimeEval.h"

#include "utils.h"

#pragma once

class Domain {

public: 
	// Constructor 
	Domain(eslocal domain_index, eslocal use_dynamic_1_no_dynamic_0);
	Domain();

	// Domain specific variables
	eslocal domain_global_index;
	eslocal domain_prim_size;
	eslocal USE_DYNAMIC;
	eslocal USE_KINV;
	eslocal USE_HFETI;
	eslocal DOFS_PER_NODE;
	eslocal isOnACC;


	// Matrices and vectors of the cluster 
	SparseMatrix B0; 
	SparseMatrix B0t;
	SparseMatrix B0_comp;
	SparseMatrix B0t_comp;
	SEQ_VECTOR <eslocal> B0_comp_map_vec;

	SparseMatrix B0Kplus; 
	SparseMatrix B0Kplus_comp; 

	SparseMatrix B0KplusB1_comp;
	SparseMatrix Kplus_R_B1_comp;


	SparseMatrix B1Kplus;
	//SparseMatrix B1KplusB1t;
	
	SparseMatrix B1; 
	SparseMatrix B1t; 
	
	SparseMatrix   B1_comp;
	SparseMatrix   B1t_comp; 
	SEQ_VECTOR< eslocal >  lambda_map_sub;
	map <eslocal, eslocal> my_lamdas_map_indices;
	SEQ_VECTOR< double >B1_scale_vec;

	SparseMatrix B1_comp_dom;
	SparseMatrix B1t_comp_dom; 
	SEQ_VECTOR <eslocal> lambda_map_sub_local;

	SparseSolver Kplus;
	SparseSolver KplusF;
	SEQ_VECTOR <double> f; 
	SEQ_VECTOR <double> vec_c;

	SparseMatrix Kplus_R; 
	SparseMatrix R; 
	SparseMatrix K;
	//SparseMatrix K_non_sym;
	SparseMatrix M; 
	SparseMatrix Prec; 

	SEQ_VECTOR <eslocal>	map_vector_e0;
	SEQ_VECTOR <eslocal>	map_vector;


	SEQ_VECTOR <SEQ_VECTOR <double> > coordinates; 
	//vector <vector <eslocal> > elements;

	SEQ_VECTOR <eslocal> fix_nodes;
	SEQ_VECTOR <eslocal> fix_dofs;

	// variables to export results 
	SEQ_VECTOR <eslocal>	number_of_nodes_in_global0;
	SEQ_VECTOR <eslocal>	map_vector_local2global0;
	SEQ_VECTOR <eslocal>	nodeMulti;
	SEQ_VECTOR <double> ux; 
	SEQ_VECTOR <double> uy; 
	SEQ_VECTOR <double> uz; 

	SEQ_VECTOR <double> up0; 
	SEQ_VECTOR <double> BtLambda_i;
	SEQ_VECTOR <double> norm_vec;
	double norm_c; 
	double norm_f; 

	// temporary variables 
	SEQ_VECTOR <double> compressed_tmp;

	// variables for dynamic 
	double dynamic_timestep;
	double dynamic_beta;
	double dynamic_gama; 

	// CUDA 
	double * cuda_pinned_buff; 
	//float  * cuda_pinned_buff_fl; 
	// END - CUDA 

	// Methods of the class
	void SetDomain(eslocal USE_HFETI, eslocal use_dynamic_1_no_dynamic_0);

	void K_regularization( );
	void K_regularizationFromR ( ); 

	void CreateKplus_R ( ); 

	void multKplusLocal( SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, eslocal x_in_vector_start_index, eslocal y_out_vector_start_index );
	void multKplusLocal( SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out );
	void multKplusLocal( SEQ_VECTOR <double> & x_in_y_out); 

	void get_kernel_from_K();

	//dynamic
	void SetDynamicParameters(double set_dynamic_timestep, double set_dynamic_beta, double set_dynamic_gama);
};


