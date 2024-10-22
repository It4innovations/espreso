
#ifndef SRC_ANALYSIS_PATTERN_APPLY_H_
#define SRC_ANALYSIS_PATTERN_APPLY_H_

#include "math/primitives/vector_dense.h"
#include "analysis/math/vector_distributed.h"

namespace espreso {

template <typename Matrix, typename T>
struct ApplyMatrix {
    virtual ~ApplyMatrix() {}

    virtual void apply(Matrix &m, Vector_Distributed<Vector_Dense, T> &y, const T &alpha, const T &beta, const Vector_Distributed<Vector_Dense, T> &x) =0;
};


}

#endif /* SRC_ANALYSIS_PATTERN_APPLY_H_ */
