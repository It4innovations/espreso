
#ifndef ASSEMBLER_ASSEMBLER_H_
#define ASSEMBLER_ASSEMBLER_H_

#include "esbasis.h"
#include "esoutput.h"
#include "esmesh.h"
#include "esbem.h"
#include "../../libespreso/feti4i.h"
#include "../solver/essolver.h"

namespace espreso {

struct FEM {
	FEM(Mesh &mesh): mesh(mesh) { };

	Mesh &mesh;
};

struct BEM {
	BEM(Mesh &mesh, Mesh &surface): mesh(mesh), surface(surface) { };

	Mesh &mesh;
	Mesh &surface;
};

struct API {
	API(APIMesh *mesh): mesh(mesh), K(NULL) { };

	APIMesh *mesh;
	SparseCSRMatrix<eslocal> *K;

	eslocal size;
	double *rhs;
	eslocal *l2g;

	eslocal dirichlet_size;
	eslocal *dirichlet_indices;
	double *dirichlet_values;

	int neighbours_size;
	int *neighbours;

	eslocal indexing;
};


class AssemblerBase {

public:
	virtual void init() = 0;
	virtual void pre_solve_update() {};
	virtual void post_solve_update() {};
	virtual void solve(std::vector<std::vector<double> > &solution) = 0;
	virtual void finalize() = 0;

	virtual ~AssemblerBase() {};
};

template <class TInput>
class Assembler: public AssemblerBase {

protected:
	Assembler(TInput &input): _input(input), _verbose(true), _timeStatistics("Solver Overall Timing") {};

	virtual size_t subdomains();
	virtual size_t DOFs();

	TInput _input;

	bool _verbose;
	TimeEval _timeStatistics;

};

}

#include "assembler.hpp"

#endif /* ASSEMBLER_ASSEMBLER_H_ */
