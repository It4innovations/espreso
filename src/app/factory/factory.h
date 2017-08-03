
#ifndef APP_FACTORY_FACTORY_H_
#define APP_FACTORY_FACTORY_H_

#include <vector>
#include "async/Dispatcher.h"

namespace espreso {

struct Instance;
struct Physics;
class LinearSolver;
class Assembler;
class TimeStepSolver;
class LoadStepSolver;
class Mesh;
namespace output {
class Store;
class AsyncStore;
class ResultStoreList;
}

struct GlobalConfiguration;
struct OutputConfiguration;
struct LoadStepSettingsBase;

class FactoryLoader {

public:
	virtual ~FactoryLoader();

	void preprocessMesh();
	virtual size_t loadSteps() const =0;

	virtual LoadStepSolver* getLoadStepSolver(size_t step, Mesh *mesh, output::Store *store) =0;

	template<class TLoadStepSettings>
	const TLoadStepSettings& getLoadStepsSettings(size_t step, const std::map<size_t, TLoadStepSettings*> &setting) const
	{
		if (setting.find(step + 1) == setting.end()) {
			printError("Missing setting for LOAD STEP " + std::to_string(step + 1));
		}
		return *setting.find(step + 1)->second;
	}

	LinearSolver* getLinearSolver(const LoadStepSettingsBase &settings, Instance *instance) const;

protected:
	void printError(const std::string &error) const;

	std::vector<Instance*> _instances;
	std::vector<Physics*> _physics;
	std::vector<LinearSolver*> _linearSolvers;
	std::vector<Assembler*> _assemblers;
	std::vector<TimeStepSolver*> _timeStepSolvers;
	std::vector<LoadStepSolver*> _loadStepSolvers;
};


class Factory {

public:
	Factory(const GlobalConfiguration &configuration);
	~Factory();

	void solve();
	void finalize();

protected:
	void initAsync(const OutputConfiguration &configuration);
	void loadPhysics(const GlobalConfiguration &configuration);
	void setOutput(const OutputConfiguration &configuration);

	FactoryLoader* createFactoryLoader(const GlobalConfiguration &configuration);

	Mesh *_mesh;
	output::ResultStoreList* _storeList;
	FactoryLoader *_loader;
	std::vector<LoadStepSolver*> _loadSteps;

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
