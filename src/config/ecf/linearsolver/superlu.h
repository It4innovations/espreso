
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_SUPERLU_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_SUPERLU_H_

#include "config/description.h"

namespace espreso {

struct SuperLUConfiguration: public ECFDescription {

    int np_row;
    int np_col;

    SuperLUConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_LINEARSOLVER_SUPERLU_H_ */
