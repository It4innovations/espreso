
#ifndef SRC_PHYSICS_COMPOSER_DISTRIBUTED_FACES_EDGES_UNIFORM_DISTRIBUTED_COMPOSER_H_
#define SRC_PHYSICS_COMPOSER_DISTRIBUTED_FACES_EDGES_UNIFORM_DISTRIBUTED_COMPOSER_H_

#include "distributed.composer.h"

namespace espreso {

struct DistributedAssemblerData;

class FacesEdgesUniformDistributedComposer: public DistributedComposer {

public:
	FacesEdgesUniformDistributedComposer(Kernel *kernel, DistributedAssemblerData *data, int fDOFs, int eDOFs);

	void init();

protected:
	void _initDOFMap();
	void _buildPatterns();
	void _buildDirichlet();

	int _fDOFs, _eDOFs;
};

}

#endif /* SRC_PHYSICS_COMPOSER_DISTRIBUTED_FACES_EDGES_UNIFORM_DISTRIBUTED_COMPOSER_H_ */
