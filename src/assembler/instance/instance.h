
#ifndef SRC_ASSEMBLER_INSTANCE_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_INSTANCE_H_

#include "esbasis.h"
#include "esmesh.h"
#include "essolver.h"

namespace espreso {

struct Instance {

	Instance(const Mesh &mesh): _mesh(mesh), _timeStatistics("Solver Overall Timing") {};

	virtual void init() = 0;

	virtual void pre_solve_update(std::vector<std::vector<double> > &solution) {};
	virtual void solve(std::vector<std::vector<double> > &solution) = 0;
	virtual void post_solve_update(std::vector<std::vector<double> > &solution) {};

	virtual void finalize() = 0;

	virtual ~Instance() {};

protected:
	const Mesh &_mesh;
	TimeEval _timeStatistics;
};

}


#endif /* SRC_ASSEMBLER_INSTANCE_INSTANCE_H_ */
