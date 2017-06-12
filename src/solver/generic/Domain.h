

#include "../generic/SparseMatrix.h"
#include "../specific/sparsesolvers.h"

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
using std::map;
using std::make_pair;

#include "../generic/utils.h"
#include "../../assembler/instance.h"


#pragma once

namespace espreso {
	
class Domain {

public:

	// Constructor
	Domain(const ESPRESOSolver &configuration, Instance *instance_in, eslocal domain_index, eslocal USE_HTFETI_in);

	// Methods of the class
	void SetDomain();

	void multKplusLocal( SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, eslocal x_in_vector_start_index, eslocal y_out_vector_start_index );
	void multKplusLocal( SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out );
	void multKplusLocal( SEQ_VECTOR <double> & x_in_y_out);

	void multKplusLocalCore( SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out );
	void multKplusLocalCore( SEQ_VECTOR <double> & x_in_y_out);

    const ESPRESOSolver &configuration;
	Instance 		    *instance;

	SparseMatrix &K;

	SparseMatrix &Kplus_R;
	SparseMatrix &Kplus_R2;
	SparseMatrix Kplus_Rb;
	SparseMatrix Kplus_Rb2;

	SparseMatrix &_RegMat;

	SEQ_VECTOR <double> &f;

	SparseMatrix B0;
	SparseMatrix B1;

	// Domain specific variables
	eslocal domain_global_index;
	eslocal domain_prim_size;
	eslocal USE_KINV;
	eslocal USE_HFETI;
	eslocal isOnACC;

	eslocal domain_index;
	bool	enable_SP_refinement;


	// Matrices and vectors of the cluster
	SparseMatrix B0t;
	SparseMatrix B0_comp;
	SparseMatrix B0t_comp;
	SEQ_VECTOR <eslocal> B0_comp_map_vec;

	SparseMatrix B0Kplus;
	SparseMatrix B0Kplus_comp;

	SparseMatrix B0KplusB1_comp;
	SparseMatrix Kplus_R_B1_comp;


	SparseMatrix B1Kplus;
	SparseMatrix B1t;
	SparseMatrix B1t_DirPr;
	SEQ_VECTOR <eslocal> B1t_Dir_perm_vec;
	SEQ_VECTOR< eslocal >  lambda_map_sub;
	map <eslocal, eslocal> my_lamdas_map_indices;
	SEQ_VECTOR< double >B1_scale_vec;

	SparseMatrix B1_comp_dom;
	SparseMatrix B1t_comp_dom;
	SEQ_VECTOR <eslocal> lambda_map_sub_local;

	SparseSolverCPU Kplus;

	SparseSolverCPU KplusF;
	SEQ_VECTOR <double> vec_c;
	SEQ_VECTOR <double> vec_lb;



	SparseMatrix R;
	SparseMatrix T;
	SparseMatrix IminusRRt;

	// Matrix and coeficient for regularization

	SparseMatrix M;
	SparseMatrix Prec;

	SEQ_VECTOR <eslocal>	map_vector_e0;
	SEQ_VECTOR <eslocal>	map_vector;

	SEQ_VECTOR <eslocal> 	fix_nodes;
	SEQ_VECTOR <eslocal> 	fix_dofs;

	// variables to export results
	SEQ_VECTOR <eslocal>	number_of_nodes_in_global0;
	SEQ_VECTOR <eslocal>	map_vector_local2global0;
	SEQ_VECTOR <eslocal>	nodeMulti;
	SEQ_VECTOR <double> 	ux;
	SEQ_VECTOR <double> 	uy;
	SEQ_VECTOR <double> 	uz;

	SEQ_VECTOR <double> up0;
	SEQ_VECTOR <double> BtLambda_i;
	SEQ_VECTOR <double> norm_vec;
	double norm_c;
	double norm_f;

	// temporary variables
	SEQ_VECTOR <double> compressed_tmp;
	SEQ_VECTOR <double> compressed_tmp2;

	// CUDA
	double * cuda_pinned_buff;
	float  * cuda_pinned_buff_fl;
	// END - CUDA


};

}

