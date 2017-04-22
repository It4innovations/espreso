

#include "../generic/SparseMatrix.h"
#include "../specific/sparsesolvers.h"
#include "../specific/densesolvers.h"

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
//	Domain(const ESPRESOSolver &configuration);
	Domain(const ESPRESOSolver &configuration, Instance *instance_in, eslocal domain_index, eslocal USE_HTFETI_in);
//	Domain(const ESPRESOSolver &configuration, eslocal domain_index, eslocal use_dynamic_1_no_dynamic_0);

//    ClusterBase(const ESPRESOSolver &configuration, Instance *instance_in):
//     	configuration(configuration),
//			instance(instance_in),


    const ESPRESOSolver &configuration;
	Instance *instance;
	// ************************

//	Instance(size_t domains, const std::vector<int> &neighbours);
//	Instance(Instance &other, Matrices &share);
//	~Instance();
//
//	size_t domains;
//	std::vector<size_t> 	&domainDOFCount;
//	std::vector<Property> 	&properties;
//	std::vector<int> 		neighbours;
//
//	std::vector<SparseMatrix> 			&K, &N1, &N2, &RegMat;

	SparseMatrix &K;

	SparseMatrix &Kplus_R;
	SparseMatrix &Kplus_R2;
	SparseMatrix Kplus_Rb;
	SparseMatrix Kplus_Rb2;

	SparseMatrix &_RegMat;


//	std::vector<SparseMatrix> 			&M;
//	std::vector<std::vector<double> > 	&R, &f;

	SEQ_VECTOR <double> &f;

//	// matrices for Hybrid FETI constraints
//	std::vector<SparseMatrix> 			&B0;
//	std::vector<std::vector<esglobal> > &B0subdomainsMap; // TODO: not needed

	SparseMatrix B0;

//	// matrices for FETI constraints
//	std::vector<SparseMatrix> &B1;

	SparseMatrix B1;

//	std::vector<std::vector<esglobal> > &B1subdomainsMap; // TODO: not needed
//	std::vector<std::vector<esglobal> > &B1clustersMap; // TODO: get it directly
//
//	std::vector<std::vector<double> > 	&B1c, &LB, &B1duplicity;
//
//	std::vector<SparseMatrix> 			&inequality;
//	std::vector<std::vector<double> > 	&inequalityC;
//
//	// blocks types of B1
//	enum CONSTRAINT {
//		DIRICHLET,
//		EQUALITY_CONSTRAINTS,
//		INEQUALITY_CONSTRAINTS,
//	};

//	std::vector<size_t> &block;
//
//
//	std::vector<std::vector<double> > primalSolution;
//	std::vector<std::vector<double> > dualSolution;
//
//	std::vector<Solution*> solutions;


	// ************************




	// Domain specific variables
	eslocal domain_global_index;
	eslocal domain_prim_size;
	eslocal USE_DYNAMIC;
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
	//SparseMatrix B1KplusB1t;
	SparseMatrix B1t;
	SparseMatrix B1t_DirPr;
	SEQ_VECTOR <eslocal> B1t_Dir_perm_vec;
	//SparseMatrix   B1_comp;
	//SparseMatrix   B1t_comp;
	SEQ_VECTOR< eslocal >  lambda_map_sub;
	map <eslocal, eslocal> my_lamdas_map_indices;
	SEQ_VECTOR< double >B1_scale_vec;

	SparseMatrix B1_comp_dom;
	SparseMatrix B1t_comp_dom;
	SEQ_VECTOR <eslocal> lambda_map_sub_local;

//	SparseSolverAcc Kplus;

	SparseSolverCPU Kplus;
	//DenseSolverCPU Kplus;


	SparseSolverCPU KplusF;
	SEQ_VECTOR <double> vec_c;
	SEQ_VECTOR <double> vec_lb;



	SparseMatrix R;
	SparseMatrix T;

	// Matrix and coeficient for regularization

	//SparseMatrix K_non_sym;
	SparseMatrix M;
	SparseMatrix Prec;

	SEQ_VECTOR <eslocal>	map_vector_e0;
	SEQ_VECTOR <eslocal>	map_vector;


	//SEQ_VECTOR <SEQ_VECTOR <double> > coordinates;
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
	SEQ_VECTOR <double> compressed_tmp2;

	// CUDA
	double * cuda_pinned_buff;
	float  * cuda_pinned_buff_fl;
	// END - CUDA

	// Methods of the class
	void SetDomain();

	void multKplusLocal( SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, eslocal x_in_vector_start_index, eslocal y_out_vector_start_index );
	void multKplusLocal( SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out );
	void multKplusLocal( SEQ_VECTOR <double> & x_in_y_out);


};

}

