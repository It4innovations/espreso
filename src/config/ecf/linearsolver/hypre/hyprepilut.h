
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPILUT_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPILUT_H_

#include "config/description.h"

namespace espreso {

struct HYPREPilutConfiguration: public ECFDescription {

	int max_iter;
	double drop_tol;
	int row_size;
	
	HYPREPilutConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPILUT_H_ */
