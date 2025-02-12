
#ifndef SRC_PHYSICS_COMPOSER_FETI_FETICOMPOSER_H_
#define SRC_PHYSICS_COMPOSER_FETI_FETICOMPOSER_H_

#include "physics/composer/composer.h"
#include "math/domainindices.h"
#include <vector>

namespace espreso {

class FETIConfiguration;
class MatrixIJVFETI;
struct FETIAssemblerData;
struct FETIVector;
template <typename TEBoundaries, typename TEData> class serializededata;

class FETIComposer: public Composer {

public:
    FETIComposer(const FETIConfiguration &configuration, Kernel *kernel, ModuleOpt *opt, FETIAssemblerData *data);
    ~FETIComposer();

    void assemble(const Builder &builder);

    virtual void buildB0FromCorners(MatrixIJVFETI &B0) =0;

protected:
    virtual void synchronize(const Builder &builder) =0;

    const FETIConfiguration &_configuration;
    FETIAssemblerData *_data;

    serializededata<esint, DI> *_DOFMap;

//    FETISolverConfiguration &_configuration;

    std::vector<std::vector<esint> > _KPermutation, _RHSPermutation;
    std::vector<esint> _domainDOFsSize, _domainDirichletSize;
    std::vector<int> _BEMDomain;
};

}


#endif /* SRC_PHYSICS_COMPOSER_FETI_FETICOMPOSER_H_ */
