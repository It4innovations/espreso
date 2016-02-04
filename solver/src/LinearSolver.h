/*
 * LinearSolver.h
 *
 *  Created on: Sep 24, 2015
 *      Author: lriha
 */

#ifndef SOLVER_SRC_LINEARSOLVER_H_
#define SOLVER_SRC_LINEARSOLVER_H_

#include "esconfig.h"
#include "esmesh.h"
#include "essolver.h"


//#include <omp.h>
//#include "mpi.h"
//#include "mkl.h"
//
//#include <string>
//#include <sstream>
//#include <iostream>
//#include <vector>
//#include <fstream>
//#include <algorithm>
//#include <math.h>
//#include <iomanip>
//#include <map>
//
//#include "stdlib.h"
//#include "stdio.h"
//#include "string.h"
//
//#include <ctime>
//#include <stack>
//#include <time.h>
//std::stack<clock_t> tictoc_stack;
//
//using std::vector;
//using std::cout;
//using std::map;
//using std::make_pair;
//
//#include <cilk/cilk.h>
//#include <cilk/cilk_api.h>
//
//#include "SparseMatrix.h"
//#include "FEM_Assembler.h"
//#include "SparseDSS.h"
//#include "SparseSolverCPU.h"
//#include "TimeEval.h"
//#include "Domain.h"
//#include "Cluster.h"
//#include "IterSolver.h"
//#include "ConfigFile.h"
//
//#include "utils.h"

class LinearSolver {
public:

	LinearSolver();

	virtual ~LinearSolver();

	void setup( eslocal rank, eslocal size, bool SINGULAR );

	void init(
			const mesh::Mesh &mesh,

			std::vector < SparseMatrix >	& K_mat,
			std::vector < SparseMatrix >	& T_mat,
			std::vector < SparseMatrix >	& B1_mat,
			std::vector < SparseMatrix >	& B0_mat,

			std::vector < std::vector <eslocal> >	& lambda_map_sub_B1,
			std::vector < std::vector <eslocal> >	& lambda_map_sub_B0,
			std::vector < std::vector <eslocal> >	& lambda_map_sub_clst,
			std::vector < std::vector <double> >	& B1_duplicity,

			std::vector < std::vector <double > >	& f_vec,
			std::vector < std::vector <double > >	& vec_c,

			const std::vector < std::vector <eslocal > >	& fix_nodes,

			std::vector < eslocal > & neigh_clusters
	);

	void init(
			std::vector < SparseMatrix >	& K_mat,
			std::vector < SparseMatrix >	& B1_mat,
			std::vector < SparseMatrix >	& B0_mat,

			std::vector < std::vector <eslocal> >	& lambda_map_sub_B1,
			std::vector < std::vector <eslocal> >	& lambda_map_sub_B0,
			std::vector < std::vector <eslocal> >	& lambda_map_sub_clst,
			std::vector < std::vector <double> >	& B1_duplicity,

			std::vector < std::vector <double > >	& f_vec,
			std::vector < std::vector <double > >	& vec_c,

			std::vector < eslocal > & neigh_clusters
	);

	void Preprocessing( std::vector < std::vector < eslocal > > & lambda_map_sub );

	void Solve( std::vector < std::vector <double > >	& f_vec, vector < vector < double > > & prim_solution );

	void Postprocessing ();

	void finilize();

	void CheckSolution( vector < vector < double > > & prim_solution );

	void set_B1(
			std::vector < SparseMatrix >				& B1_mat,
			std::vector < std::vector <double> >        & B1_duplicity);

	void set_B0(
			std::vector < SparseMatrix >				& B0_mat );

	void set_R(
			const mesh::Mesh &mesh
	);

  	void set_R_from_K();

	eslocal DOFS_PER_NODE;

private:

	eslocal MPI_rank;
	eslocal MPI_size;
	eslocal number_of_subdomains_per_cluster;

	bool 	SINGULAR;
	bool 	KEEP_FACTORS;
  	bool 	R_from_mesh;

	TimeEval timeEvalMain; //(string("ESPRESO Solver Overal Timing"));

	Cluster cluster;

	IterSolver solver;



};

#endif /* SOLVER_SRC_LINEARSOLVER_H_ */
