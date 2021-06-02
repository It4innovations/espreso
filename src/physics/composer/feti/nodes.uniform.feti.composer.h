
#ifndef SRC_PHYSICS_COMPOSER_FETI_NODES_UNIFORM_FETI_COMPOSER_H_
#define SRC_PHYSICS_COMPOSER_FETI_NODES_UNIFORM_FETI_COMPOSER_H_

#include "feti.composer.h"

namespace espreso {

struct FETIAssemblerData;
class NodesUniformDistributedComposer;

class NodesUniformFETIComposer: public FETIComposer {

public:
	NodesUniformFETIComposer(Kernel *kernel, FETIAssemblerData *data, int DOFs);

	void init();

protected:
	void synchronize(const Builder &builder);

	void _initDOFMap();
	void _buildPatterns();
	void _buildKFEMPattern(esint domain);
	void _buildKBEMPattern(esint domain);
	void _buildDirichlet();
	void _buildMortars();
	void _buildInequality();

	int _DOFs;
};

}



#endif /* SRC_PHYSICS_COMPOSER_FETI_NODES_UNIFORM_FETI_COMPOSER_H_ */
