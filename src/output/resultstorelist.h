
#ifndef SRC_OUTPUT_RESULTSTORELIST_H_
#define SRC_OUTPUT_RESULTSTORELIST_H_

#include "resultstore.h"

namespace espreso {
namespace store {

class ResultStoreList: public ResultStore {

public:
	ResultStoreList(const OutputConfiguration &output): ResultStore(output, NULL, "") { };
	~ResultStoreList() { std::for_each(_results.begin(), _results.end(), [] (ResultStore *rs) { delete rs; } ); }

	void add(ResultStore *rs) { _results.push_back(rs); }

	void storeGeometry(size_t timeStep = -1)
	{
		std::for_each(_results.begin(), _results.end(), [&] (ResultStore *rs) { rs->storeGeometry(timeStep); } );
	}

	void storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType)
	{
		std::for_each(_results.begin(), _results.end(), [&] (ResultStore *rs) { rs->storeProperty(name, properties, eType); } );
	}

	void storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType)
	{
		std::for_each(_results.begin(), _results.end(), [&] (ResultStore *rs) { rs->storeValues(name, dimension, values, eType); } );
	}

	void finalize()
	{
		std::for_each(_results.begin(), _results.end(), [&] (ResultStore *rs) { rs->finalize(); } );
	}

protected:
	std::vector<ResultStore*> _results;
};

}
}

#endif /* SRC_OUTPUT_RESULTSTORELIST_H_ */
