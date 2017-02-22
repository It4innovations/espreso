
#ifndef APP_FACTORY_FACTORY_H_
#define APP_FACTORY_FACTORY_H_

#include <vector>

namespace espreso {

struct Solver;
struct Physics;
struct Instance;
class LinearSolver;
namespace store { class ResultStore; }


struct GlobalConfiguration;
struct Results;
struct OldInstance;
struct Mesh;

struct Factory {

	Factory(const GlobalConfiguration &configuration);
	~Factory();

	void solve();
	void finalize();

	void check(const Results &configuration);

	double norm() const;

	std::vector<Solver*> loadSteps;
	store::ResultStore* store;

	OldInstance *instance;
	Mesh *mesh;

private:
	void meshPreprocessing();

	std::vector<Solver*> _solvers;
	std::vector<Physics*> _physics;
	std::vector<Instance*> _instances;
	std::vector<LinearSolver*> _linearSolvers;

	bool _newAssembler;

	std::vector<std::vector<double> > _solution;
};

}



#endif /* APP_FACTORY_FACTORY_H_ */
