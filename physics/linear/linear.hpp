
#include "linear.h"

namespace physics {

template <MatrixComposer TMatrixComposer>
void Linear<TMatrixComposer>::init()
{
	this->_timeStatistics.SetName("Linear Elasticity Solver Overall Timing");
	this->_timeStatistics.totalTime.AddStartWithBarrier();
	std::cout.precision(15);

	TimeEvent timeKasm("Create K and RHS");
	timeKasm.AddStart();

	_K.resize(subdomains());
	_M.resize(subdomains());
	_f.resize(subdomains());
	for (size_t s = 0; s < subdomains(); s++) {
		// TODO: set dynamics
		KMf(s, false);

		if (this->_verbose && this->_mesh.rank() == 0) {
			std::cout << s << " " ;
		}
	}
	if (this->_verbose && this->_mesh.rank() == 0) {
		std::cout << std::endl;
	}

	timeKasm.AddEndWithBarrier();
	this->_timeStatistics.AddEvent(timeKasm);

	TimeEvent timeLocalB("Create local B");
	timeLocalB.AddStart();

	_localB.resize(subdomains());
	_globalB.resize(subdomains());

	_lambda_map_sub_B1.resize(subdomains());
	_lambda_map_sub_B0.resize(subdomains());
	_B1_duplicity.resize(subdomains());
	_vec_c.resize(subdomains());
	localB();

	timeLocalB.AddEndWithBarrier();
	this->_timeStatistics.AddEvent(timeLocalB);

	TimeEvent timeGlobalB("Create global B");
	timeGlobalB.AddStart();


	globalB();
//	 for (int i = 0; i < vec_c.size(); i++)
//		 vec_c[i].resize(B1_duplicity[i].size(), 0.0);

	timeGlobalB.AddEndWithBarrier();
	this->_timeStatistics.AddEvent(timeGlobalB);


	TimeEvent timeBforces("Fill right hand side");
	timeBforces.AddStart();

	RHS();

	timeBforces.AddEndWithBarrier();
	this->_timeStatistics.AddEvent(timeBforces);

	TimeEvent timeLSconv(string("Linear Solver - preprocessing"));
	timeLSconv.AddStart();

	lin_solver.DOFS_PER_NODE = this->DOFs();
	lin_solver.setup(this->_mesh.rank(), this->_mesh.size(), true);

	initSolver();
//
//	 if (BEM) {
//
//		 lin_solver.init(
//
//			_instance.surf_mesh(),
//
//			K_mat_ls,
//
//			B1_mat_ls,
//			B0_mat_ls,
//
//			lambda_map_sub_B1,
//			lambda_map_sub_B0,
//			lambda_map_sub_clst,
//			B1_duplicity,
//
//			f_vec,
//			vec_c,
//
//			fix_nodes,
//			l2g_vec,
//
//			neigh_clusters
//
//		);
//	 } else {
//		 lin_solver.init(
//
//			_instance.mesh(),
//
//			K_mat_ls,
//
//			B1_mat_ls,
//			B0_mat_ls,
//
//			lambda_map_sub_B1,
//			lambda_map_sub_B0,
//			lambda_map_sub_clst,
//			B1_duplicity,
//
//			f_vec,
//			vec_c,
//
//			fix_nodes,
//			l2g_vec,
//
//			neigh_clusters
//
//		);
//
//	 }
//
//	 timeLSconv.AddEndWithBarrier();
//	 timeEvalMain.AddEvent(timeLSconv);
}

template <MatrixComposer TMatrixComposer>
void Linear<TMatrixComposer>::pre_solve_update()
{

}

template <MatrixComposer TMatrixComposer>
void Linear<TMatrixComposer>::post_solve_update()
{

}

template <MatrixComposer TMatrixComposer>
void Linear<TMatrixComposer>::solve()
{

}

template <MatrixComposer TMatrixComposer>
void Linear<TMatrixComposer>::finalize()
{

}


}
