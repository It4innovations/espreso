
#ifndef SRC_ANALYSIS_COMPOSER_NODES_UNIFORM_FETI_H_
#define SRC_ANALYSIS_COMPOSER_NODES_UNIFORM_FETI_H_

#include "math2/utils/utils_feti.h"

#include <vector>
#include <functional>

namespace espreso {

class UniformNodesFETIPattern {

public:
	struct RegionInfo {
		esint size;
		std::vector<esint> K, f;
	};

	UniformNodesFETIPattern(): _DOFs(1), _DOFMap(nullptr) {}

	void initUpper();
	void initFull();

protected:
	int _DOFs;

	serializededata<esint, DI> *_DOFMap;

	std::vector<RegionInfo> elements; // RegionInfo per domain
	std::vector<std::vector<RegionInfo> > bregion; // RegionInfo per domain per boundary region

	void init();

//	std::vector<std::vector<esint> > _KPermutation, _RHSPermutation;
//	std::vector<esint> _domainDOFsSize, _domainDirichletSize;
//	std::vector<int> _BEMDomain;
};

}

#endif /* SRC_ANALYSIS_COMPOSER_NODES_UNIFORM_FETI_H_ */
