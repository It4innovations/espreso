
#ifndef SRC_PHYSICS_COMPOSER_DISTRIBUTED_DISTRIBUTED_COMPOSER_OPT_H_
#define SRC_PHYSICS_COMPOSER_DISTRIBUTED_DISTRIBUTED_COMPOSER_OPT_H_

#include "distributed.composer.h"

namespace espreso {

class DistributedComposerOpt: public DistributedComposer {

public:
    using DistributedComposer::DistributedComposer;

    void assemble(const Builder &builder);
};

}

#endif /* SRC_PHYSICS_COMPOSER_DISTRIBUTED_DISTRIBUTED_COMPOSER_OPT_H_ */
