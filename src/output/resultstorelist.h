
#ifndef SRC_OUTPUT_RESULTSTORELIST_H_
#define SRC_OUTPUT_RESULTSTORELIST_H_

#include <algorithm>

#include "store.h"

namespace espreso {

class ResultStoreList: public Store {

public:
	ResultStoreList(const OutputConfiguration &output): Store(output) { };
	~ResultStoreList() { std::for_each(_results.begin(), _results.end(), [] (Store *rs) { delete rs; } ); }

	void add(Store *rs) { _results.push_back(rs); }

	virtual void updateMesh()
	{
		std::for_each(_results.begin(), _results.end(), [&] (Store *rs) { rs->updateMesh(); } );
	}

	virtual void storeSettings(const Step &step)
	{
		std::for_each(_results.begin(), _results.end(), [&] (Store *rs) { rs->storeSettings(step); } );
	}

	virtual void storeFETIData(const Step &step, const Instance &instance)
	{
		std::for_each(_results.begin(), _results.end(), [&] (Store *rs) { rs->storeFETIData(step, instance); } );
	}

	virtual void storeSolution(const Step &step, const std::vector<Solution*> &solution, const std::vector<std::pair<ElementType, Property> > &properties)
	{
		std::for_each(_results.begin(), _results.end(), [&] (Store *rs) { rs->storeSolution(step, solution, properties); } );
	}

	virtual void storeSubSolution(const Step &step, const std::vector<Solution*> &solution, const std::vector<std::pair<ElementType, Property> > &properties)
	{
		std::for_each(_results.begin(), _results.end(), [&] (Store *rs) { rs->storeSubSolution(step, solution, properties); } );
	}

	virtual void finalize()
	{
		std::for_each(_results.begin(), _results.end(), [&] (Store *rs) { rs->finalize(); } );
	}

protected:
	std::vector<Store*> _results;
};

}

#endif /* SRC_OUTPUT_RESULTSTORELIST_H_ */
