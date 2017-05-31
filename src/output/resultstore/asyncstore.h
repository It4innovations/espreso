
#ifndef SRC_OUTPUT_RESULTSTORE_ASYNCSTORE_H_
#define SRC_OUTPUT_RESULTSTORE_ASYNCSTORE_H_

#include "resultstore.h"
#include "asyncstoreexecutor.h"
#include "../../assembler/step.h"
#include "async/Module.h"

namespace espreso {

class Mesh;

namespace output {

class AsyncStore : public ResultStore, private async::Module<AsyncStoreExecutor, OutputConfiguration, Param> {

public:
	AsyncStore(const OutputConfiguration &configuration, const Mesh *mesh, MeshInfo::InfoMode mode = MeshInfo::EMPTY)
		: ResultStore(configuration, mesh, mode), _finalized(false), _bufferSize(0), _headerStore(0L) {};
	virtual ~AsyncStore() { delete _headerStore; };

	/** Called by ASYNC on all ranks, needs to set at least the executor */
	void setUp() { setExecutor(_executor); };

	void init(const Mesh *mesh);

	virtual std::string store(const std::string &name, const RegionData &regionData);

	virtual std::string linkClusters(const std::string &root, const std::string &name, const RegionData &regionData);
	virtual void linkSteps(const std::string &name, const std::vector<std::pair<Step, std::vector<std::string> > > &steps);

	virtual void finalize();

	/** Called by ASYNC on all ranks, needs to finalize the executor */
	void tearDown() { _executor.finalize(); }

private:
	/** The executor doing the real work */
	AsyncStoreExecutor _executor;

	bool _finalized;

	/** Size of the packed data */
	size_t _bufferSize;

	/** A result store that stores header files */
	ResultStore* _headerStore;
};

}
}

#endif /* SRC_OUTPUT_RESULTSTORE_ASYNCSTORE_H_ */
