
#ifndef ASSEMBLER_LINEAR_LINEAR_H_
#define ASSEMBLER_LINEAR_LINEAR_H_

#include "../gluing/gluing.h"

namespace assembler {

template <MatrixComposer TMatrixComposer>
class Linear: public Gluing<TMatrixComposer> {

public:
	virtual ~Linear() {};

	void init();
	void pre_solve_update();
	void post_solve_update();
	void solve();
	void finalize();

protected:
	Linear(const mesh::Mesh &mesh): Gluing<TMatrixComposer>(mesh) {};

	// FEM specific
	virtual size_t DOFs() = 0;
	virtual void inertia(std::vector<double> &inertia) = 0;
	virtual void C(DenseMatrix &C) = 0;
	virtual double CP() = 0;
	virtual double rho() = 0;

	// Matrices for Linear Solver
	std::vector<SparseMatrix> _K, _M;
	// RHS
	std::vector<std::vector<double> > _f;

	// Result
	vector<vector<double> > _prim_solution;

	LinearSolver _lin_solver;

private:
	void KMf(size_t part, bool dynamics);
	void KeMefe(
			DenseMatrix &Ke, DenseMatrix &Me, std::vector<double> &fe,
			DenseMatrix &Ce, const mesh::Element *e, size_t part, bool dynamics);
	void integrate(
			DenseMatrix &Ke, DenseMatrix &Me, std::vector<double> &fe,
			SparseVVPMatrix<eslocal> &K, SparseVVPMatrix<eslocal> &M, std::vector<double> &f,
			const mesh::Element *e, bool dynamics);

	void RHS();
	void initSolver();
	void saveResult();
};


}

#include "linear.hpp"
#include "fem.hpp"
#include "bem.hpp"


#endif /* ASSEMBLER_LINEAR_LINEAR_H_ */
