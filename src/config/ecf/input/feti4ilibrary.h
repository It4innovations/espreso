
#ifndef SRC_CONFIG_ECF_INPUT_FETI4ILIBRARY_H_
#define SRC_CONFIG_ECF_INPUT_FETI4ILIBRARY_H_

#include "config/ecf/linearsolver/feti.h"

namespace espreso {

struct FETI4ILibraryConfiguration: public ECFDescription {

    size_t domains;
    FETIConfiguration solver;

    FETI4ILibraryConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_FETI4ILIBRARY_H_ */
