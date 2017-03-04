
#ifndef SRC_OUTPUT_RESULTSTORELIST_H_
#define SRC_OUTPUT_RESULTSTORELIST_H_

#include "resultstore.h"

namespace espreso {
namespace output {

class ResultStoreList: public Store {

public:
	ResultStoreList(const OutputConfiguration &output): Store(output) { };
	~ResultStoreList() { std::for_each(_results.begin(), _results.end(), [] (ResultStore *rs) { delete rs; } ); }

	void add(ResultStore *rs) { _results.push_back(rs); }

	virtual void storeSettings(const Step &step)
	{
		std::for_each(_results.begin(), _results.end(), [&] (ResultStore *rs) { rs->storeSettings(step); } );
	}

	virtual void storeSettings(size_t steps)
	{
		std::for_each(_results.begin(), _results.end(), [&] (ResultStore *rs) { rs->storeSettings(steps); } );
	}

	virtual void storeSettings(const std::vector<size_t> &steps)
	{
		std::for_each(_results.begin(), _results.end(), [&] (ResultStore *rs) { rs->storeSettings(steps); } );
	}

	virtual void storeSolution(const Step &step, const std::vector<Solution*> &solution)
	{
		std::for_each(_results.begin(), _results.end(), [&] (ResultStore *rs) { rs->storeSolution(step, solution); } );
	}

	virtual void finalize()
	{
		std::for_each(_results.begin(), _results.end(), [&] (ResultStore *rs) { rs->finalize(); } );
	}

protected:
	std::vector<ResultStore*> _results;
};

}
}

#endif /* SRC_OUTPUT_RESULTSTORELIST_H_ */
