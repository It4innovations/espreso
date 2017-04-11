
#ifndef SRC_OUTPUT_RESULTSTORE_ASYNCSTOREEXECUTOR_H_
#define SRC_OUTPUT_RESULTSTORE_ASYNCSTOREEXECUTOR_H_

#include "../store.h"
#include "async/ExecInfo.h"

namespace espreso {
namespace output {

class AsyncStoreExecutor {

public:
	AsyncStoreExecutor() : _store(0L) {}
	~AsyncStoreExecutor() { delete _store; }

	void execInit(const async::ExecInfo &info, const OutputConfiguration &config);

	void exec(const async::ExecInfo &info, const Step &step);

	void finalize()
	{
		if (_store)
			_store->finalize();
	}

private:
	/** The real store */
	Store* _store;
};

}
}

#endif /* SRC_OUTPUT_RESULTSTORE_ASYNCSTOREEXECUTOR_H_ */
