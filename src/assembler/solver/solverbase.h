
#ifndef SRC_ASSEMBLER_SOLVER_SOLVERBASE_H_
#define SRC_ASSEMBLER_SOLVER_SOLVERBASE_H_

#include <cstddef>
#include <vector>
#include <string>

#include "../instance.h"

namespace espreso {

class Mesh;
struct Step;
class TimeEval;

class SolverBase
{
public:
	SolverBase(const std::string &name, const std::string &physicsName, Mesh *mesh, double duration);

	virtual void run(Step &step) =0;

	const std::string& name() const { return _name; }
	double duration() const { return _duration; }

	virtual ~SolverBase();

protected:
	std::string _name;
	Mesh *_mesh;
	double _duration;

	TimeEval *_timeStatistics;
};

}

#endif /* SRC_ASSEMBLER_SOLVER_SOLVERBASE_H_ */
