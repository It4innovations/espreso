/*
 * LinearSolver.h
 *
 *  Created on: Sep 24, 2015
 *      Author: lriha
 */

#ifndef SOLVER_GENERIC_LINEARSOLVER_H_
#define SOLVER_GENERIC_LINEARSOLVER_H_

#include "../essolver.h"
#include "../specific/itersolvers.h"
#include "esconfig.h"
#include "esmesh.h"

namespace espreso {

class LinearSolver {
public:

	LinearSolver();

	virtual ~LinearSolver();

	void setup( eslocal rank, eslocal size, bool SINGULAR );

	void init(
			const Mesh &mesh,

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

			const std::vector < int > & neigh_clusters
	);

	void Preprocessing( std::vector < std::vector < eslocal > > & lambda_map_sub );

	void Solve( std::vector < std::vector <double > >	& f_vec, std::vector < std::vector < double > > & prim_solution );

	void Postprocessing ();

	void finilize();

	void CheckSolution( std::vector < std::vector < double > > & prim_solution );

	void set_B1(
			std::vector < SparseMatrix >				& B1_mat,
			std::vector < std::vector <double> >        & B1_duplicity);

	void set_B0(
			std::vector < SparseMatrix >				& B0_mat );

	void set_R(
			const Mesh &mesh
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

}

#endif /* SOLVER_GENERIC_LINEARSOLVER_H_ */
