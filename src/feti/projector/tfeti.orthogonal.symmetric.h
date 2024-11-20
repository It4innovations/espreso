
#ifndef SRC_FETI_PROJECTOR_TFETI_ORTHOGONAL_SYMMETRIC_H_
#define SRC_FETI_PROJECTOR_TFETI_ORTHOGONAL_SYMMETRIC_H_

#include "projector.h"
#include "dualgraph.h"

#include <map>

namespace espreso {

/*
 * R: KxR : block diagonal
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
struct TFETIOrthogonalSymmetric: public Projector<T> {
    TFETIOrthogonalSymmetric(FETI<T> &feti);
    ~TFETIOrthogonalSymmetric();

    void set(const step::Step &step);
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

    struct DomainInfo {
        int domain, koffset, kernels;

        DomainInfo() = default;
        DomainInfo(int domain, int koffset, int kernels): domain(domain), koffset(koffset), kernels(kernels) {}

        bool operator< (const DomainInfo &other) const { return domain <  other.domain; }
        bool operator<=(const DomainInfo &other) const { return domain <= other.domain; }
        bool operator!=(const DomainInfo &other) const { return domain != other.domain; }
    };

    struct NeighborDomainInfo: DomainInfo {
        struct CIndices { int offset, count; };
        std::vector<CIndices> cindices;
        int ncols;

        NeighborDomainInfo() = default;

        NeighborDomainInfo& operator=(const DomainInfo &other) {
            this->domain = other.domain;
            this->koffset = other.koffset;
            this->kernels = other.kernels;
            return *this;
        }
    };

    std::vector<DomainInfo> dinfo;
    std::vector<std::vector<DomainInfo> > dualGraph;
    std::map<int, NeighborDomainInfo> upinfo;
    std::vector<std::map<int, NeighborDomainInfo> > downinfo;

    DualGraph dual;
};

}



#endif /* SRC_FETI_PROJECTOR_TFETI_ORTHOGONAL_SYMMETRIC_H_ */
