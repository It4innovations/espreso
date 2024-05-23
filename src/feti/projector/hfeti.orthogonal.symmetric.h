
#ifndef SRC_FETI_PROJECTOR_HFETI_ORTHOGONAL_SYMMETRIC_H_
#define SRC_FETI_PROJECTOR_HFETI_ORTHOGONAL_SYMMETRIC_H_

#include "projector.h"

#include <map>

namespace espreso {

/*
 * R: KxR : single block
 * B: LxK : from primal to dual
 *
 * e = Rt * f : R
 * G = Rt * Bt: RxL: from dual to kernels
 *
 * y = Q * x = Gt * inv(GGt) * G * x
 *
 *     Gx = G * x         :: x      -> Gx     :      L -> R
 * iGGtGx = inv(GGt) * Gx :: Gx     -> iGGtGx : totalR -> R
 *      y = Gt * iGGtGx   :: iGGtGx -> y      :      R -> L
 */

template <typename T>
struct HFETIOrthogonalSymmetric: public Projector<T> {
    HFETIOrthogonalSymmetric(FETI<T> &feti);
    ~HFETIOrthogonalSymmetric();

    void update(const step::Step &step);

protected:
    void _computeDualGraph();
    void _setG();
    void _setGGt();
    void _updateG();
    void _updateGGt();

    using Projector<T>::feti;
    using Projector<T>::kernel;
    using Projector<T>::e;
    using Projector<T>::G;
    using Projector<T>::Gt;
    using Projector<T>::GGt;
    using Projector<T>::invGGt;
    using Projector<T>::invL;
    using Projector<T>::invU;

    using Projector<T>::Gx;
    using Projector<T>::iGGtGx;

    size_t domainOffset;
    size_t GGtDataOffset, GGtDataSize, GGtNnz;

    struct ClusterInfo {
        int cluster, koffset, kernels;

        ClusterInfo() = default;
        ClusterInfo(int cluster, int koffset, int kernels): cluster(cluster), koffset(koffset), kernels(kernels) {}

        bool operator< (const ClusterInfo &other) const { return cluster <  other.cluster; }
        bool operator<=(const ClusterInfo &other) const { return cluster <= other.cluster; }
        bool operator!=(const ClusterInfo &other) const { return cluster != other.cluster; }
    };
    struct NeighborClusterInfoInfo: ClusterInfo {
        struct CIndices { int offset, count; };
        std::vector<CIndices> cindices;
        int ncols;

        NeighborClusterInfoInfo() = default;

        NeighborClusterInfoInfo& operator=(const ClusterInfo &other) {
            this->cluster = other.cluster;
            this->koffset = other.koffset;
            this->kernels = other.kernels;
            return *this;
        }
    };
    std::vector<ClusterInfo> cinfo;
    std::vector<std::vector<ClusterInfo> > dualGraph;
    std::vector<NeighborClusterInfoInfo> neighInfo;
};

}


#endif /* SRC_FETI_PROJECTOR_HFETI_ORTHOGONAL_SYMMETRIC_H_ */
