
#ifndef SRC_PHYSICS_COMPOSER_FETI_NODES_UNIFORM_API_COMPOSER_H_
#define SRC_PHYSICS_COMPOSER_FETI_NODES_UNIFORM_API_COMPOSER_H_

#include "nodes.uniform.feti.composer.h"

namespace espreso {

struct FETISolverData;

class NodesUniformAPIComposer: public NodesUniformFETIComposer {

public:
    NodesUniformAPIComposer(const FETIConfiguration &configuration, int DOFs);

    void fill(FETISolverData &data);

    const serializededata<esint, DI>* DOFMap() const { return _DOFMap; };

protected:
    void _setDecomposition(FETISolverData &data);
    void _setPattern(FETISolverData &data, esint domain);
};

}



#endif /* SRC_PHYSICS_COMPOSER_FETI_NODES_UNIFORM_API_COMPOSER_H_ */
