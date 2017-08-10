
#ifndef SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_LOADSTEPSOLVER_H_
#define SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_LOADSTEPSOLVER_H_

#include <string>

namespace espreso {

struct Step;
class TimeStepSolver;
class Assembler;
enum Matrices: int;

class LoadStepSolver {

	friend class LinearTimeStep;
	friend class NewtonRhapson;

public:
	LoadStepSolver(const std::string &description, TimeStepSolver &timeStepSolver, double duration);
	virtual ~LoadStepSolver() {}

	void run(Step &step);

	std::string description() const;
	double duration() const;

protected:
	virtual Matrices updateStructuralMatrices(Step &step, Matrices matrices) =0;
	virtual Matrices reassembleStructuralMatrices(Step &step, Matrices matrices) =0;

	virtual void initLoadStep(Step &step);
	virtual bool hasNextTimeStep(Step &step);
	virtual void runNextTimeStep(Step &step) =0;
	virtual void processTimeStep(Step &step) =0;
	virtual void finalizeLoadStep(Step &step);

	std::string _description;
	TimeStepSolver &_timeStepSolver;
	Assembler &_assembler;
	double _duration;

	double _startTime;
	double _precision;
};

}


#endif /* SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_LOADSTEPSOLVER_H_ */
