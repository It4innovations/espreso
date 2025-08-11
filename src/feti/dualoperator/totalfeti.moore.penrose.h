
#ifndef SRC_FETI_DUALOPERATOR_TOTALFETI_MOOREPENROSE_H_
#define SRC_FETI_DUALOPERATOR_TOTALFETI_MOOREPENROSE_H_

#include "dualoperator.h"
#include "math/wrappers/math.spsolver.h"

namespace espreso {


template <typename T>
class TotalFETIMoorePenrose: public DualOperator<T> {
public:
    TotalFETIMoorePenrose(FETI<T> &feti);
    ~TotalFETIMoorePenrose();

    void info();
    void set(const step::Step &step);
    void update(const step::Step &step);

    // y = B * K+ * Bt * x
    void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);
    void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y);

    // y = K+(f - Bt * x)
    void toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y);

protected:
    using DualOperator<T>::feti;
    using DualOperator<T>::d;

    void _apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);

    std::vector<Vector_Dense<T> > Btx, KplusBtx;
};

}

#endif /* SRC_FETI_DUALOPERATOR_TOTALFETI_MOOREPENROSE_H_ */
