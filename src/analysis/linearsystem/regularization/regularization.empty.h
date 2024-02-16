
#ifndef SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_EMPTY_H_
#define SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_EMPTY_H_

#include "feti/feti.h"

namespace espreso {

template <typename T>
struct RegularizationEmpty {
    static void set(FETI<T> &feti);
    static void update(FETI<T> &feti);
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_EMPTY_H_ */
