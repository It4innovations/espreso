/*
 * LinearSolver.h
 *
 *  Created on: Sep 24, 2015
 *      Author: lriha
 */

#ifndef SOLVER_GENERIC_LINEARSOLVER_H_
#define SOLVER_GENERIC_LINEARSOLVER_H_

#include "../specific/itersolvers.h"

#include "../../assembler/constraints/constraints.h"
#include "../../assembler/old_physics/assembler.h"
#include "../../assembler/instance.h"


namespace espreso {

struct Instance;

class LinearSolver {
public:

	LinearSolver(Instance *instance, const ESPRESOSolver &configuration);
	LinearSolver(const ESPRESOSolver &configuration, OldPhysics &physics, Constraints &constraints); // to be removed

	void init();

	void update(Matrices matrices);
	void run();

	virtual ~LinearSolver();

	void setup();

	void init(const std::vector<int> &neighbours);

//	void init(
//			const Mesh &mesh,
//
//			std::vector < SparseMatrix >	& K_mat,
//			std::vector < SparseMatrix >	& T_mat,
//			std::vector < SparseMatrix >	& B1_mat,
//			std::vector < SparseMatrix >	& B0_mat,
//
//			std::vector < std::vector <eslocal> >	& lambda_map_sub_B1,
//			std::vector < std::vector <eslocal> >	& lambda_map_sub_B0,
//			std::vector < std::vector <eslocal> >	& lambda_map_sub_clst,
//			std::vector < std::vector <double> >	& B1_duplicity,
//
//			std::vector < std::vector <double > >	& f_vec,
//			std::vector < std::vector <double > >	& vec_c,
//
//			const std::vector < std::vector <eslocal > >	& fix_nodes,
//
//			const std::vector < int > & neigh_clusters
//	);

	void Preprocessing( std::vector < std::vector < eslocal > > & lambda_map_sub );

	void Solve( std::vector < std::vector <double > >	& f_vec, std::vector < std::vector < double > > & prim_solution );
	void Solve( std::vector < std::vector <double > >	& f_vec, std::vector < std::vector < double > > & prim_solution, std::vector < std::vector < double > > & dual_solution );

	void Postprocessing ();

	void finilize();

	void CheckSolution( std::vector < std::vector < double > > & prim_solution );

	void set_B1(
			const std::vector < SparseMatrix >				& B1_mat,
			const std::vector < std::vector <double> >        & B1_duplicity);

	void set_B0(
			const std::vector < SparseMatrix >				& B0_mat );

//	void set_R(
//			const Mesh &mesh
//	);

	void set_R_from_K();

	Instance *instance;
	const ESPRESOSolver &configuration;

	// TODO: to be removed
	OldPhysics *physics;
	Constraints *constraints;


	TimeEval timeEvalMain; //(string("ESPRESO Solver Overal Timing"));

private:
	eslocal number_of_subdomains_per_cluster;

	bool 	SINGULAR;
	bool 	KEEP_FACTORS;

	Cluster *cluster;
	IterSolver *solver;



};

}

#endif /* SOLVER_GENERIC_LINEARSOLVER_H_ */
