
#ifndef APP_FACTORY_FACTORY_H_
#define APP_FACTORY_FACTORY_H_

#include <vector>

namespace espreso {

struct Solver;
struct NewPhysics;
struct NewInstance;
class LinearSolver;
namespace store { class ResultStore; }


struct GlobalConfiguration;
struct Results;
struct Instance;
struct Mesh;

struct Factory {

	Factory(const GlobalConfiguration &configuration);
	~Factory();

	void solve();
	void check(const Results &configuration);

	double norm() const;

	Solver *solver;
	std::vector<NewPhysics*> physics;
	std::vector<NewInstance*> instances;
	std::vector<LinearSolver*> linearSolvers;
	store::ResultStore* store;

	Instance *instance;
	Mesh *mesh;

private:
	std::vector<std::vector<double> > _solution;
};

}



#endif /* APP_FACTORY_FACTORY_H_ */
