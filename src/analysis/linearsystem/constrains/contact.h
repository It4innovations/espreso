
#ifndef SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_CONTACT_H_
#define SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_CONTACT_H_

#include "analysis/math/vector_distributed.h"
#include "esinfo/ecfinfo.h"
#include "feti/feti.h"

namespace espreso {

template <typename T>
struct Contact {

    void set(const step::Step &step, FETI<T> &feti);
    void update(const step::Step &step, FETI<T> &feti);
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_CONTACT_H_ */
