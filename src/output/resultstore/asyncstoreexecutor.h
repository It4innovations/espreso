
#ifndef SRC_OUTPUT_RESULTSTORE_ASYNCSTOREEXECUTOR_H_
#define SRC_OUTPUT_RESULTSTORE_ASYNCSTOREEXECUTOR_H_

#include "../resultstore.h"
#include "async/ExecInfo.h"

namespace espreso {

struct Param
{ };

class AsyncStoreExecutor {

public:
	AsyncStoreExecutor() : _store(0L) {}
	~AsyncStoreExecutor() { delete _store; }

	void execInit(const async::ExecInfo &info, const OutputConfiguration &config);

	void exec(const async::ExecInfo &info, const Param &param);

	void finalize()
	{
		if (_store)
			_store->finalize();
	}

private:
	/** The real store */
	ResultStore* _store;
};

}

#endif /* SRC_OUTPUT_RESULTSTORE_ASYNCSTOREEXECUTOR_H_ */
