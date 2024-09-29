
#ifndef SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_EQUALITYCONSTRAINS_H_
#define SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_EQUALITYCONSTRAINS_H_

#include "analysis/math/vector_distributed.h"
#include "feti/feti.h"

namespace espreso {

template <typename T>
struct EqualityConstrains {

    void set(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet);
    void update(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet);

    void enforce(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet);

protected:
    std::vector<size_t> doffset;
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_EQUALITYCONSTRAINS_H_ */
