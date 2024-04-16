
#ifndef SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_CLUSTER_FACES_H_
#define SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_CLUSTER_FACES_H_

#include "esinfo/ecfinfo.h"
#include "feti/feti.h"

namespace espreso {

template <typename T>
struct ClusterFacesGluing {

    void set(const step::Step &step, FETI<T> &feti);
    void update(const step::Step &step, FETI<T> &feti);
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_CLUSTER_FACES_H_ */
