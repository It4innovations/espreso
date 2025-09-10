
#ifndef SRC_FETI_PROJECTOR_HFETI_CONJUGATE_SYMMETRIC_H_
#define SRC_FETI_PROJECTOR_HFETI_CONJUGATE_SYMMETRIC_H_

#include "projector.h"
#include "dualgraph.h"

namespace espreso {

template <typename T>
struct HFETIConjugateSymmetric: public Projector<T> {
    HFETIConjugateSymmetric(FETI<T> &feti);
    ~HFETIConjugateSymmetric();

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
    std::vector<int> nonzeros, distributed;
};

}




#endif /* SRC_FETI_PROJECTOR_HFETI_CONJUGATE_SYMMETRIC_H_ */
