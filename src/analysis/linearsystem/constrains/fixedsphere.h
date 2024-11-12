
#ifndef SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_FIXEDSPHERE_H_
#define SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_FIXEDSPHERE_H_

#include "analysis/math/vector_distributed.h"
#include "esinfo/ecfinfo.h"
#include "feti/feti.h"

namespace espreso {

template <typename T>
struct FixedSphere {

    FixedSphere(): interval(1), cmapsize(0) {}

    void set(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet);
    void update(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet);

    std::vector<std::vector<esint> > cindex;
    std::vector<size_t> dsize;
    size_t interval, cmapsize;

private:
    void _store(const Point &center, double radius);
    void _store(const std::vector<Point> &coords, const std::vector<std::pair<Point, double> > &gaps);
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_FIXEDSPHERE_H_ */
