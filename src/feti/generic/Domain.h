

#include "feti/generic/SparseMatrix.h"
#include "feti/specific/sparsesolvers.h"
#include "feti/specific/densesolvers.h"

#include <map>

using std::vector;
using std::map;

#include "feti/generic/utils.h"
#include "feti/dataholder.h"

#pragma once

namespace espreso {
	
class Domain {

public:

	// Constructor
	Domain(const FETIConfiguration &configuration, DataHolder *instance_in, esint domain_index, esint USE_HTFETI_in);

	// Methods of the class
	void SetDomain();

	void multKplusLocal( SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, esint x_in_vector_start_index, esint y_out_vector_start_index );
	void multKplusLocal( SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out );
	void multKplusLocal( SEQ_VECTOR <double> & x_in_y_out);

	void multKplusLocalCore( SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out );
	void multKplusLocalCore( SEQ_VECTOR <double> & x_in_y_out);

	const FETIConfiguration &configuration;
	DataHolder 		    *instance;

	SparseMatrix &origK;
	SparseMatrix &K;

	SparseMatrix &Kplus_R;
	SparseMatrix &Kplus_R2;
	SparseMatrix Kplus_Rb;
	SparseMatrix Kplus_Rb2;

	SparseMatrix &Kplus_origR;
	SparseMatrix &Kplus_origR2;

	SparseMatrix &_RegMat;

	SEQ_VECTOR <double> &f;

	SparseMatrix B0;
	SparseMatrix B1;

	// Domain specific variables
	esint domain_global_index;
	esint domain_prim_size;
	esint USE_KINV;
	esint USE_HFETI;
	esint isOnACC;

	esint domain_index;
	bool	enable_SP_refinement;


	// Matrices and vectors of the cluster
	SparseMatrix B0t;
	SparseMatrix B0_comp;
	SparseMatrix B0t_comp;
	SEQ_VECTOR <esint> B0_comp_map_vec;

	SparseMatrix B0Kplus;
	SparseMatrix B0Kplus_comp;

	SparseMatrix B0KplusB1_comp;
	SparseMatrix Kplus_R_B1_comp;


	SparseMatrix B1Kplus;
	SparseMatrix B1t;
	SparseMatrix B1t_DirPr;
	SEQ_VECTOR <esint> B1t_Dir_perm_vec;
	SEQ_VECTOR< esint >  lambda_map_sub;
	map <esint, esint> my_lamdas_map_indices;
	SEQ_VECTOR< double >B1_scale_vec;

	SparseMatrix B1_comp_dom;
	SparseMatrix B1t_comp_dom;
	SEQ_VECTOR <esint> lambda_map_sub_local;

//	SparseSolverAcc Kplus;

#ifdef BEM4I_TO_BE_REMOVED
	DenseSolverCPU Kplus;
#else
	SparseSolverCPU Kplus;
#endif

	SparseSolverCPU KplusF;
	SEQ_VECTOR <double> vec_c;
	SEQ_VECTOR <double> vec_lb;



	SparseMatrix R;
	SparseMatrix T;
	SparseMatrix IminusRRt;

	// Matrix and coeficient for regularization

	SparseMatrix M;
	SparseMatrix Prec;

	SEQ_VECTOR <esint>	map_vector_e0;
	SEQ_VECTOR <esint>	map_vector;

	SEQ_VECTOR <esint> 	fix_nodes;
	SEQ_VECTOR <esint> 	fix_dofs;

	// variables to export results
	SEQ_VECTOR <esint>	number_of_nodes_in_global0;
	SEQ_VECTOR <esint>	map_vector_local2global0;
	SEQ_VECTOR <esint>	nodeMulti;
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

