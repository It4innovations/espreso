
#ifndef SRC_OUTPUT_RESULTSTORE_ASYNCSTORE_H_
#define SRC_OUTPUT_RESULTSTORE_ASYNCSTORE_H_

#include "store.h"
#include "asyncstoreexecutor.h"
#include "../../assembler/step.h"
#include "async/Module.h"

namespace espreso {

class Mesh;

namespace output {

class AsyncStore : public Store, private async::Module<AsyncStoreExecutor, OutputConfiguration, Step> {

public:
	AsyncStore(const OutputConfiguration &configuration): Store(configuration), _finalized(false) {};
	virtual ~AsyncStore() {};

	/** Called by ASYNC on all ranks, needs to set at least the executor */
	void setUp() { std::cout << "setUp" << std::endl; setExecutor(_executor); };

	void init(const Mesh *mesh, const std::string &path);

	virtual void storeSettings(const Step &step) {};
	virtual void storeSettings(size_t steps) {};
	virtual void storeSettings(const std::vector<size_t> &steps) {};

	virtual void storeFETIData(const Step &step, const Instance &instance) {};

	virtual void storeSolution(const Step &step, const std::vector<Solution*> &solution) {};
	virtual void finalize();

	/** Called by ASYNC on all ranks, needs to finalize the executor */
	void tearDown() { _executor.finalize(); }

private:
	/** The executor doing the real work */
	AsyncStoreExecutor _executor;

	bool _finalized;
};

}
}

#endif /* SRC_OUTPUT_RESULTSTORE_ASYNCSTORE_H_ */
