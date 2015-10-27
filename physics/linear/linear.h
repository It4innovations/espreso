
#ifndef PHYSICS_LINEAR_LINEAR_H_
#define PHYSICS_LINEAR_LINEAR_H_

#include "../physics.h"

namespace physics {

template <MatrixComposer TMatrixComposer>
class Linear: public Physics<TMatrixComposer> {

public:
	virtual ~Linear() {};

	void init();
	void pre_solve_update();
	void post_solve_update();
	void solve();
	void finalize();

protected:
	Linear(const mesh::Mesh &mesh): Physics<TMatrixComposer>(mesh) {};

	// FEM specific
	virtual size_t DOFs() = 0;
	virtual void inertia(std::vector<double> &inertia) = 0;
	virtual void C(DenseMatrix &C) = 0;
	virtual double CP() = 0;
	virtual double rho() = 0;

	// Matrices for Linear Solver
	std::vector<SparseMatrix> _K, _M, _localB, _globalB;
	std::vector<SparseIJVMatrix<eslocal> > _B0, _B1;
	// RHS
	std::vector<std::vector<double> > _f;

	// Result
	vector<vector<double> > _prim_solution;

	// Description ??
	std::vector<std::vector<eslocal> > _lambda_map_sub_B1;
	std::vector<std::vector<eslocal> > _lambda_map_sub_B0;
	std::vector<std::vector<eslocal> > _lambda_map_sub_clst;
	std::vector<std::vector<double> > _B1_duplicity;
	std::vector<std::vector<double> > _vec_c;
	std::vector<eslocal> _neighClusters;

	LinearSolver _lin_solver;

private:
	size_t subdomains();
	void KMf(size_t part, bool dynamics);
	void KeMefe(
			DenseMatrix &Ke, DenseMatrix &Me, std::vector<double> &fe,
			DenseMatrix &Ce, const mesh::Element *e, size_t part, bool dynamics);
	void integrate(
			DenseMatrix &Ke, DenseMatrix &Me, std::vector<double> &fe,
			SparseVVPMatrix<eslocal> &K, SparseVVPMatrix<eslocal> &M, std::vector<double> &f,
			const mesh::Element *e, bool dynamics);

	void localB();
	void globalB();
	void RHS();
	void initSolver();
	void saveResult();
};


}

#include "linear.hpp"
#include "fem.hpp"
#include "bem.hpp"


#endif /* PHYSICS_LINEAR_LINEAR_H_ */
