
#ifndef SRC_FETI_PROJECTOR_TFETI_ORTHOGONAL_SYMMETRIC_WITHFACTORS_H_
#define SRC_FETI_PROJECTOR_TFETI_ORTHOGONAL_SYMMETRIC_WITHFACTORS_H_

#include "projector.h"

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
struct TFETIOrthogonalSymmetricWithFactors: public Projector<T> {
    TFETIOrthogonalSymmetricWithFactors(FETI<T> &feti);
    ~TFETIOrthogonalSymmetricWithFactors();

    void info();
    void update(const step::Step &step);

    void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);
    void apply_e(const Vector_Kernel<T> &x, Vector_Dual<T> &y);
    void apply_R(const Vector_Kernel<T> &x, std::vector<Vector_Dense<T> > &y);
    void apply_Ra(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y);
    void apply_invU(const Vector_Kernel<T> &x, Vector_Kernel<T> &y);
    void apply_invL(const Vector_Kernel<T> &x, Vector_Kernel<T> &y);
    void apply_GtinvU(const Vector_Kernel<T> &x, Vector_Dual<T> &y);
    void apply_invLG(const Vector_Dual<T> &x, Vector_Kernel<T> &y);

protected:
    void _computeDualGraph();
    void _setG();
    void _setGGt();
    void _updateG();
    void _updateGGt();

    void _applyG(const Vector_Dual<T> &in, Vector_Kernel<T> &out);
    void _applyInvGGt(const Vector_Kernel<T> &in, Vector_Dense<T> &out);
    void _applyInvL(const Vector_Kernel<T> &in, Vector_Dense<T> &out);
    void _applyInvU(const Vector_Kernel<T> &in, Vector_Dense<T> &out);
    void _applyGt(const Vector_Dense<T> &in, const T &alpha, Vector_Dual<T> &out);
    void _applyR(const Vector_Dense<T> &in, std::vector<Vector_Dense<T> > &out);

    void _print(const step::Step &step);

    using Projector<T>::feti;
    using Projector<T>::e;

    Matrix_CSR<T> G, Gt, GGt;
    Matrix_Dense<T> invL, invU;

    Vector_Kernel<T> Gx; // we need whole vector
    Vector_Dense<T> iGGtGx; // only local part is sufficient

    size_t domainOffset;
    size_t GGtDataOffset, GGtDataSize, GGtNnz;

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



#endif /* SRC_FETI_PROJECTOR_TFETI_ORTHOGONAL_SYMMETRIC_WITHFACTORS_H_ */
