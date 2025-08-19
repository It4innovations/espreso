
#ifndef SRC_CONFIG_ECF_OPERATIONS_SCHUR_CSX_DNY_PASTIX_H_
#define SRC_CONFIG_ECF_OPERATIONS_SCHUR_CSX_DNY_PASTIX_H_

#include "config/description.h"

namespace espreso {

struct SchurCsxDnyPastixConfig: public ECFDescription {

    enum struct AUTOBOOL {
        AUTO,
        TRUE,
        FALSE
    };

    AUTOBOOL use_gpu;

    SchurCsxDnyPastixConfig();

};

}

#endif /* #ifndef SRC_CONFIG_ECF_OPERATIONS_SCHUR_CSX_DNY_PASTIX_H_ */
