
#ifndef SRC_FETI_PROJECTOR_HFETI_ORTHOGONAL_SYMMETRIC_H_
#define SRC_FETI_PROJECTOR_HFETI_ORTHOGONAL_SYMMETRIC_H_

#include "projector.h"
#include "dualgraph.h"

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

    void set(const step::Step &step);
    void update(const step::Step &step);

    void orthonormalizeKernels(const step::Step &step);

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

    size_t GGtDataOffset, GGtDataSize;

    DualGraph dual;
};

}


#endif /* SRC_FETI_PROJECTOR_HFETI_ORTHOGONAL_SYMMETRIC_H_ */
