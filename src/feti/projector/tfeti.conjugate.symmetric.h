
#ifndef SRC_FETI_PROJECTOR_TFETI_CONJUGATE_SYMMETRIC_H_
#define SRC_FETI_PROJECTOR_TFETI_CONJUGATE_SYMMETRIC_H_

#include "projector.h"

#include <map>

namespace espreso {

/*
 * R: KxR : block diagonal
 * B: LxK : from primal to dual
 * F: KxK : dual operator
 *
 * e = Rt * f : R
 * G = Rt * Bt: RxL: from dual to kernels
 *
 * y = Q * x = Gt * inv(GFGt) * G * F * x
 *
 *     Gx = G * x         :: x      -> Gx     :      L -> R
 * iGGtGx = inv(GGt) * Gx :: Gx     -> iGGtGx : totalR -> R
 *      y = Gt * iGGtGx   :: iGGtGx -> y      :      R -> L
 */

template <typename T>
struct TFETIConjugateSymmetric: public Projector<T> {
    TFETIConjugateSymmetric(FETI<T> &feti);
    ~TFETIConjugateSymmetric();

    void info();
    void update(const step::Step &step);

    void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);
    void apply_e(const Vector_Kernel<T> &x, Vector_Dual<T> &y);
    void apply_Ra(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y);

protected:
    void _computeDualGraph();
    void _setG();
    void _setGFGt();
    void _updateG();
    void _updateGFGt();

    void _applyG(const Vector_Dual<T> &in, Vector_Kernel<T> &out);
    void _applyInvGFGt(const Vector_Kernel<T> &in, Vector_Dense<T> &out);
    void _applyGt(const Vector_Dense<T> &in, const T &alpha, Vector_Dual<T> &out);
    void _applyR(const Vector_Dense<T> &in, std::vector<Vector_Dense<T> > &out);

    void _print(const step::Step &step);

    using Projector<T>::feti;
    using Projector<T>::e;

    Matrix_CSR<T> G, Gt, GFGt;
    Matrix_Dense<T> invGFGt;

    Vector_Kernel<T> Gx; // we need whole vector
    Vector_Dense<T> iGFGtGx; // only local part is sufficient

    size_t domainOffset;
    size_t GFGtDataOffset, GFGtDataSize, GFGtNnz;

    struct DomainInfo {
        esint domain, koffset, kernels;

        DomainInfo() = default;
        DomainInfo(esint domain, esint koffset, esint kernels): domain(domain), koffset(koffset), kernels(kernels) {}

        bool operator< (const DomainInfo &other) const { return domain <  other.domain; }
        bool operator<=(const DomainInfo &other) const { return domain <= other.domain; }
        bool operator!=(const DomainInfo &other) const { return domain != other.domain; }
    };

    struct NeighborDomainInfo: DomainInfo {
        struct CIndices { esint offset, count; };
        std::vector<CIndices> cindices;
        esint ncols;

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
    std::map<esint, NeighborDomainInfo> upinfo;
    std::vector<std::map<esint, NeighborDomainInfo> > downinfo;
};

}

#endif /* SRC_FETI_PROJECTOR_TFETI_CONJUGATE_SYMMETRIC_H_ */
