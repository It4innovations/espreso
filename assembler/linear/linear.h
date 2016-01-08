
#ifndef ASSEMBLER_LINEAR_LINEAR_H_
#define ASSEMBLER_LINEAR_LINEAR_H_

#include "../gluing/gluing.h"

namespace assembler {

template <class TInput>
class Linear: public Gluing<TInput> {

public:
	virtual ~Linear() {};

	void init();
	void pre_solve_update();
	void post_solve_update();
	void solve(std::vector<std::vector<double> > &solution);
	void finalize();

	void fillAPIHolder(APIHolder *holder);

protected:
	Linear(TInput &input): Gluing<TInput>(input) {};

	// FEM specific
	virtual void inertia(std::vector<double> &inertia) = 0;
	virtual void C(DenseMatrix &C) = 0;
	virtual double CP() = 0;
	virtual double rho() = 0;

	// Matrices for Linear Solver
	std::vector<SparseMatrix> _K, _T, _M;
	// RHS
	std::vector<std::vector<double> > _f;

	LinearSolver _lin_solver;

private:
	void KMf(size_t part, bool dynamics);
	void T(size_t part);
	void RHS();
	void initSolver();

	void KeMefe(
			DenseMatrix &Ke, DenseMatrix &Me, std::vector<double> &fe,
			DenseMatrix &Ce, const mesh::Element *e, size_t part, bool dynamics);
	void integrate(
			DenseMatrix &Ke, DenseMatrix &Me, std::vector<double> &fe,
			SparseVVPMatrix<eslocal> &K, SparseVVPMatrix<eslocal> &M, std::vector<double> &f,
			const mesh::Element *e, bool dynamics);


};


}

#include "linear.hpp"
#include "fem.hpp"
#include "bem.hpp"
#include "api.hpp"


#endif /* ASSEMBLER_LINEAR_LINEAR_H_ */
