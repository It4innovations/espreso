
#ifndef SRC_PHYSICS_COMPOSER_DISTRIBUTED_NODES_UNIFORM_DISTRIBUTED_COMPOSER_H_
#define SRC_PHYSICS_COMPOSER_DISTRIBUTED_NODES_UNIFORM_DISTRIBUTED_COMPOSER_H_

#include "distributed.composer.opt.h"

namespace espreso {

struct DistributedAssemblerData;

class NodesUniformDistributedComposer: public DistributedComposerOpt {

	struct FaceID {
		esint e1, e2;

		bool operator<(const FaceID &other) {
			if (e1 == other.e1) {
				return e2 < other.e2;
			}
			return e1 < other.e1;
		}
	};

	struct EdgeID {
		esint e, n1, n2;
	};

public:
	NodesUniformDistributedComposer(Kernel *kernel, ModuleOpt *opt, DistributedAssemblerData *data, int DOFs);

	int esize(esint interval);
	int bsize(esint region, esint interval);

	void init();

protected:
	void _initDOFMap();
	void _buildPatterns();
	void _buildDirichlet();

	int _DOFs;
};

}



#endif /* SRC_PHYSICS_COMPOSER_DISTRIBUTED_NODES_UNIFORM_DISTRIBUTED_COMPOSER_H_ */
