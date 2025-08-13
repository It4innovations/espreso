
#ifndef SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_CONSTRAINS_H_
#define SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_CONSTRAINS_H_

#include "analysis/math/vector_distributed.h"
#include "feti/feti.h"
#include "equalityconstrains.h"
#include "mortar.h"
#include "fixedwall.h"
#include "fixedsphere.h"
#include "fixedtube.h"

namespace espreso {

template <typename T>
struct Constrains {

    void set(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet);
    void update(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet);

    EqualityConstrains<T> eq;
    MortarContact<T> mortar;
    FixedWall<T> fw;
    FixedSphere<T> fs;
    FixedTube<T> ft;
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_CONSTRAINS_H_ */
