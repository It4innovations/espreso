
#ifndef SRC_PHYSICS_COMPOSER_DISTRIBUTED_DISTRIBUTEDCOMPOSER_H_
#define SRC_PHYSICS_COMPOSER_DISTRIBUTED_DISTRIBUTEDCOMPOSER_H_

#include "physics/composer/composer.h"
#include <vector>

namespace espreso {

struct DistributedAssemblerData;
struct SparseVector;
template <typename TEBoundaries, typename TEData> class serializededata;

class DistributedComposer: public Composer {

public:
	DistributedComposer(Kernel *kernel, ModuleOpt *opt, DistributedAssemblerData *data);
	~DistributedComposer();

	void assemble(const Builder &builder);

protected:
	DistributedAssemblerData *_data;

	serializededata<esint, esint> *_DOFMap;

	std::vector<esint> _tKOffsets, _tRHSOffsets;
	std::vector<esint> _KPermutation, _RHSPermutation;
	std::vector<esint> _nDistribution;
};

}


#endif /* SRC_PHYSICS_COMPOSER_DISTRIBUTED_DISTRIBUTEDCOMPOSER_H_ */
