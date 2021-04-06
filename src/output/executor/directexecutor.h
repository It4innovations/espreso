
#ifndef SRC_OUTPUT_RESULT_EXECUTOR_DIRECTEXECUTOR_H_
#define SRC_OUTPUT_RESULT_EXECUTOR_DIRECTEXECUTOR_H_

#include "executor.h"

namespace espreso {

class DirectExecutor: public ResultStoreExecutor {

public:
	void updateMesh()
	{
		for (size_t i = 0; i < _resultStore.size(); i++) {
			_resultStore[i]->updateMesh();
		}
	}

	void updateMonitors()
	{
		for (size_t i = 0; i < _resultStore.size(); i++) {
			_resultStore[i]->updateMonitors();
		}
	}

	void updateSolution()
	{
		for (size_t i = 0; i < _resultStore.size(); i++) {
			_resultStore[i]->updateSolution();
		}
	}

	DirectExecutor(const Mesh &mesh): ResultStoreExecutor(mesh) {}
};

}



#endif /* SRC_OUTPUT_RESULT_EXECUTOR_DIRECTEXECUTOR_H_ */
