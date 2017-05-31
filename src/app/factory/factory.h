
#ifndef APP_FACTORY_FACTORY_H_
#define APP_FACTORY_FACTORY_H_

#include <vector>
#include "async/Dispatcher.h"

namespace espreso {

struct SolverBase;
struct Physics;
struct Instance;
class LinearSolver;
namespace output {
class AsyncStore;
class ResultStoreList;
}


struct GlobalConfiguration;
struct OutputConfiguration;
struct Results;
struct OldInstance;
struct Mesh;

struct Factory {

	Factory(const GlobalConfiguration &configuration);
	~Factory();

	void solve();
	void finalize();

	void check(const Results &configuration);

	std::vector<SolverBase*> loadSteps;
	output::ResultStoreList* store;

	OldInstance *instance;
	Mesh *mesh;

private:
	void meshPreprocessing(const OutputConfiguration &configuration);

	std::vector<SolverBase*> _solvers;
	std::vector<Physics*> _physics;
	std::vector<Instance*> _instances;
	std::vector<LinearSolver*> _linearSolvers;

	bool _newAssembler;

	std::vector<std::vector<double> > _solution;

	/**
	 * We always create the async store (even if the output is not enabled).
	 * This is a drawback of the ASYNC library but required to get synchronization
	 * right.
	 */
	output::AsyncStore* _asyncStore;

	/** The dispatcher for the I/O ranks */
	async::Dispatcher _dispatcher;
};

}



#endif /* APP_FACTORY_FACTORY_H_ */
