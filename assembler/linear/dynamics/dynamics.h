
#ifndef ASSEMBLER_LINEAR_DYNAMICS_DYNAMICS_H_
#define ASSEMBLER_LINEAR_DYNAMICS_DYNAMICS_H_

#include "../elasticity/linearelasticity.h"

namespace espreso {

template <class TInput>
class TransientElasticity: public LinearElasticity<TInput> {

public:
	TransientElasticity(TInput &input): LinearElasticity<TInput>(input),
	_beta(0.25), _gama(0.5), _timestep(1e-6), _deltaT(_timestep), _time(0) { };
	virtual ~TransientElasticity() {};

	virtual void init();
	virtual void pre_solve_update();
	virtual void post_solve_update();
	virtual void solve(std::vector<std::vector<double> > &solution);

protected:
	virtual double timeConstant() { return 1 / (_beta * _timestep * _timestep); }

	double _beta;
	double _gama;
	double _timestep;
	double _deltaT;

	size_t _time;

	std::vector<double> _constantA;
	std::vector<std::vector<double> > _u, _v, _w;    // old vectors
	std::vector<std::vector<double> > _un, _vn, _wn; // new vectors
	std::vector<std::vector<double> > _b, _tmp;
//
//private:
//
//    eslocal MPI_rank;
//    eslocal MPI_size;
//
//	// Matrices for Mesh generator and Assembler
//	std::vector < SparseCSRMatrix<eslocal> > K_mat;
//	std::vector < SparseCSRMatrix<eslocal> > M_mat;
//	std::vector < SparseIJVMatrix<eslocal> > B1_mat;
//	std::vector < SparseIJVMatrix<eslocal> > B0_mat;
//
//	// Matrices for Linear Solver
//	std::vector < SparseMatrix >			K_mat_ls;
//	std::vector < SparseMatrix >			M_mat_ls;
//	std::vector < SparseMatrix >			B1_mat_ls;
//	std::vector < SparseMatrix >			B0_mat_ls;
//
//	// Suporting data required for communication layer - all generated by generator of Matrix B
//	std::vector < std::vector <eslocal> >	lambda_map_sub_B1;
//	std::vector < std::vector <eslocal> >	lambda_map_sub_B0;
//	std::vector < std::vector <eslocal> >	lambda_map_sub_clst;
//	std::vector < std::vector <double> >	B1_duplicity;
//	std::vector < std::vector <eslocal> >	l2g_vec;				// l2g vector per cluster - not for entire problem
//
//
//	std::vector < std::vector <double > >	f_vec; 					// right hand side
//	std::vector < eslocal > domain_prim_size;						// primal sizes of all domains per cluster
//
//
//	std::vector < std::vector <eslocal > >	fix_nodes;
//
//	std::vector < eslocal > neigh_clusters;							// list of neighboring cluster - defined by their MPI rank
//
//
//	size_t partsCount; 											//number of domains per cluster
//	int DOFS_PER_NODE;											// number of degrees of freedom per one node
//
//	//Time measurement instances
//	TimeEval timeEvalMain;										// Time measurement instancefor overal solver runtime
//
//	// ESPRESO Linear Solver instance
//	LinearSolver lin_solver; // espreso solver class to be determined
//
//
//
//	// Specific data for transient elasticity
//	double dynamic_beta;
//	double dynamic_gama;
//	double dynamic_timestep;
//	double time_const;
//
//	eslocal timeStep;
//
//	std::vector < std::vector <double> > vec_u;
//	std::vector < std::vector <double> > vec_v;
//	std::vector < std::vector <double> > vec_w;
//
//	std::vector < std::vector <double> > vec_u_n;
//	std::vector < std::vector <double> > vec_v_n;
//	std::vector < std::vector <double> > vec_w_n;
//
//	std::vector < std::vector <double> > vec_b;
//	std::vector < std::vector <double> > vec_t_tmp;
//
//	std::vector <double> const_a;

};

}

#include "dynamics.hpp"

#endif /* ASSEMBLER_LINEAR_DYNAMICS_DYNAMICS_H_ */
