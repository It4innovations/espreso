
#ifndef SRC_FETI_DUALOPERATOR_TOTALFETI_MOOREPENROSE_H_
#define SRC_FETI_DUALOPERATOR_TOTALFETI_MOOREPENROSE_H_

#include "dualoperator.h"
#include "math/wrappers/math.spsolver.h"

namespace espreso {


template <typename T>
class TotalFETIMoorePenrose: public DualOperator<T> {
public:
    TotalFETIMoorePenrose(FETI<T> &feti);
    virtual ~TotalFETIMoorePenrose();

    void info() override;
    void set(const step::Step &step) override;
    void update(const step::Step &step) override;

    // y = B * K+ * Bt * x
    void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) override;
    void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y) override;

    // y = K+(f - Bt * x)
    void toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y) override;

protected:
    using DualOperator<T>::feti;
    using DualOperator<T>::d;

    void _apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);

    std::vector<Vector_Dense<T> > Btx, KplusBtx;
};

}

#endif /* SRC_FETI_DUALOPERATOR_TOTALFETI_MOOREPENROSE_H_ */
