
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPARASAILS_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPARASAILS_H_

#include "config/description.h"

namespace espreso {

struct HYPREParaSailsConfiguration: public ECFDescription {

	double threshold;
	int n_levels;
	double filter;

	enum class SYMMETRY {
		NON_INF,
		SPD,
		NON_DEF_SPD
	};
	SYMMETRY symmetry;

	double loadbal;
	int reuse;
	int logging;
	
	HYPREParaSailsConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPARASAILS_H_ */
