
#include "output/visualization/visualization.h"
#include "executor.h"

#include "basis/utilities/utils.h"

#include "output/monitors/monitoring.h"

using namespace espreso;

ResultStoreExecutor::~ResultStoreExecutor()
{
	for (size_t i = 0; i < _resultStore.size(); i++) {
		delete _resultStore[i];
	}
}

void ResultStoreExecutor::addResultStore(ResultStoreBase *resultStore)
{
	_resultStore.push_back(resultStore);
}

bool ResultStoreExecutor::storeStep()
{
	return Monitoring::storeStep() || Visualization::storeStep();
}


