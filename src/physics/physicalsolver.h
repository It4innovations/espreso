
#ifndef SRC_PHYSICS_PHYSICALSOLVER_H_
#define SRC_PHYSICS_PHYSICALSOLVER_H_

namespace espreso {

class LoadStepSolver;
class SubStepSolver;
class LinearSystem;

struct PhysicalSolver {
	static void run();

	~PhysicalSolver();
protected:
	PhysicalSolver();

	void clear();

	template <typename TPhysics>
	static void runSingle(PhysicalSolver &solver, TPhysics &configuration);
	template <typename TPhysics>
	static void runCoupled(PhysicalSolver &first, PhysicalSolver &second, TPhysics &configuration);

	LoadStepSolver *loadStepSolver;
	SubStepSolver *subStepSolver;
	LinearSystem *system;
};

}



#endif /* SRC_PHYSICS_PHYSICALSOLVER_H_ */
