
#ifndef ASSEMBLER_ASSEMBLER_H_
#define ASSEMBLER_ASSEMBLER_H_

#include "essolver.h"
#include "esmesh.h"
#include "esbem.h"

namespace assembler {

enum MatrixComposer {
	FEM,
	BEM,
	ELMER
};

template <MatrixComposer TMatrixComposer>
class Assembler {

public:
	virtual void init() = 0;
	virtual void pre_solve_update() = 0;
	virtual void post_solve_update() = 0;
	virtual void solve() = 0;
	virtual void finalize() = 0;

	virtual ~Assembler() {};

protected:
	Assembler(const mesh::Mesh &mesh);

	virtual size_t subdomains();

	const mesh::Mesh &_mesh;
	const mesh::SurfaceMesh _surface;

	bool _verbose;
	TimeEval _timeStatistics;

};

}

#include "assembler.hpp"

#endif /* ASSEMBLER_ASSEMBLER_H_ */
