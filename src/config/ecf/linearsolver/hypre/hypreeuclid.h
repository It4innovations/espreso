
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREEUCLID_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREEUCLID_H_

#include "config/description.h"

namespace espreso {

struct HYPREEuclidConfiguration: public ECFDescription {

	int levels;
	int stats;
	int memory_stats;
	double sparse_tol;
	int row_scale;
	int set_bj;
	double ilut_tol;

	
	HYPREEuclidConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREEUCLID_H_ */
