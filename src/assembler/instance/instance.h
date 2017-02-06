
#ifndef SRC_ASSEMBLER_INSTANCE_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_INSTANCE_H_

#include "../../basis/logging/timeeval.h"

namespace espreso {

class Mesh;
class OldPhysics;
class Constraints;

struct OldInstance {

	OldInstance(const Mesh &mesh): _mesh(mesh), _timeStatistics("Solver Overall Timing") {};

	virtual void init() = 0;
	virtual void solve(std::vector<std::vector<double> > &solution) = 0;
	virtual void finalize() = 0;

	virtual const OldPhysics& physics() const = 0;
	virtual const Constraints& constraints() const = 0;
	virtual OldPhysics& physics() = 0;
	virtual Constraints& constraints() = 0;

	virtual ~OldInstance() {};

protected:
	const Mesh &_mesh;
	TimeEval _timeStatistics;
};

}


#endif /* SRC_ASSEMBLER_INSTANCE_INSTANCE_H_ */
