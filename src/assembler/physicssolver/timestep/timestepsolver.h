
#ifndef SRC_ASSEMBLER_PHYSICSSOLVER_TIMESTEP_TIMESTEPSOLVER_H_
#define SRC_ASSEMBLER_PHYSICSSOLVER_TIMESTEP_TIMESTEPSOLVER_H_

#include <string>

namespace espreso {

class LinearSolver;
class LoadStepSolver;
class Assembler;
struct Step;

class TimeStepSolver {

	friend class LoadStepSolver;

public:
	TimeStepSolver(const std::string &description, Assembler &assembler);
	virtual ~TimeStepSolver() {}

	virtual void solve(Step &step, LoadStepSolver &loadStepSolver) =0;

	std::string description() const;

protected:
	std::string _description;
	Assembler &_assembler;
};

}



#endif /* SRC_ASSEMBLER_PHYSICSSOLVER_TIMESTEP_TIMESTEPSOLVER_H_ */
