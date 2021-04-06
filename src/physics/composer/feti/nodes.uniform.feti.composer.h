
#ifndef SRC_PHYSICS_COMPOSER_FETI_NODES_UNIFORM_FETI_COMPOSER_H_
#define SRC_PHYSICS_COMPOSER_FETI_NODES_UNIFORM_FETI_COMPOSER_H_

#include "feti.composer.h"

namespace espreso {

struct FETIAssemblerData;

class NodesUniformFETIComposer: public FETIComposer {

public:
	NodesUniformFETIComposer(FETIConfiguration &configuration, Kernel *kernel, FETIAssemblerData *data, int DOFs);

	void init();

	void buildB0FromCorners(MatrixIJVFETI &B0);

protected:
	void synchronize(const Builder &builder);

	void computeCornerNodes();
	void computeFixPoints();
	void computeFixPointsOnSurface();

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
