
#ifndef ASSEMBLER_ASSEMBLER_H_
#define ASSEMBLER_ASSEMBLER_H_

#include "esoutput.h"
#include "essolver.h"
#include "esmesh.h"
#include "esbem.h"
#include "../libespreso/feti4i.h"

namespace assembler {

struct FEM {
	FEM(mesh::Mesh &mesh): mesh(mesh) { };

	mesh::Mesh &mesh;
};

struct BEM {
	BEM(mesh::Mesh &mesh, mesh::SurfaceMesh &surface): mesh(mesh), surface(surface) { };

	mesh::Mesh &mesh;
	mesh::SurfaceMesh &surface;
};


// Design for testing
struct APIHolder {
	SparseCSRMatrix<eslocal> *K;
	std::vector<double> *rhs;
	std::vector<eslocal> *dirichlet_indices;
	std::vector<double> *dirichlet_values;
	std::vector<eslocal> *l2g;
	std::vector<eslocal> *neighbourRanks;
	eslocal indexing;

	~APIHolder() {
		delete K;
		delete rhs;
		delete dirichlet_indices;
		delete dirichlet_values;
		delete l2g;
		delete neighbourRanks;
	}
};

struct API {

	API(SparseCSRMatrix<eslocal> &K,
		std::vector<double> &rhs,
		std::vector<eslocal> &dirichlet_indices,
		std::vector<double> &dirichlet_values,
		std::vector<eslocal> &l2g,
		std::vector<eslocal> &neighbourRanks)
	:K(K), rhs(rhs), dirichlet_indices(dirichlet_indices), dirichlet_values(dirichlet_values),
	 l2g(l2g), neighbourRanks(neighbourRanks), indexing(0) { };

	API(APIHolder &holder): K(*holder.K), rhs(*holder.rhs),
		dirichlet_indices(*holder.dirichlet_indices), dirichlet_values(*holder.dirichlet_values),
		l2g(*holder.l2g), neighbourRanks(*holder.neighbourRanks), indexing(holder.indexing) { };

	SparseCSRMatrix<eslocal> &K;
	std::vector<double> &rhs;
	std::vector<eslocal> &dirichlet_indices;
	std::vector<double> &dirichlet_values;
	std::vector<eslocal> &l2g;
	std::vector<eslocal> &neighbourRanks;
	eslocal indexing;
};

struct API2 {
	API2() {};
	API2(APIHolder &holder):
		K(holder.K),
		rhs_size(holder.rhs->size()), rhs(holder.rhs->data()),
		dirichlet_size(holder.dirichlet_indices->size()),
		dirichlet_indices(holder.dirichlet_indices->data()),
		dirichlet_values(holder.dirichlet_values->data()),
		l2g_size(holder.l2g->size()), l2g(holder.l2g->data()),
		neighbours_size(holder.neighbourRanks->size()), neighbours(holder.neighbourRanks->data()),
		indexing(holder.indexing) { };

	SparseCSRMatrix<eslocal> *K;

	eslocal rhs_size;
	double *rhs;

	eslocal dirichlet_size;
	eslocal *dirichlet_indices;
	double *dirichlet_values;

	eslocal l2g_size;
	eslocal *l2g;

	eslocal neighbours_size;
	eslocal *neighbours;

	eslocal indexing;
};


class AssemblerBase {

public:
	virtual void init() = 0;
	virtual void pre_solve_update() = 0;
	virtual void post_solve_update() = 0;
	virtual void solve(std::vector<std::vector<double> > &solution) = 0;
	virtual void finalize() = 0;

	virtual size_t DOFs() = 0;
	virtual void fillAPIHolder(APIHolder *holder) = 0;

	virtual ~AssemblerBase() {};
};

template <class TInput>
class Assembler: public AssemblerBase {

protected:
	Assembler(TInput &input): _input(input), _verbose(true) {};

	virtual size_t subdomains();
	virtual size_t rank();
	virtual size_t size();

	TInput _input;

	bool _verbose;
	TimeEval _timeStatistics;

};

}

#include "assembler.hpp"

#endif /* ASSEMBLER_ASSEMBLER_H_ */
