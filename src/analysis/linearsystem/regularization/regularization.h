
#ifndef SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_H_
#define SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_H_

#include "esinfo/ecfinfo.h"
#include "feti/feti.h"
#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

template <typename T>
struct Regularization {

    void set(const step::Step &step, FETI<T> &feti);
    void update(const step::Step &step, FETI<T> &feti, Vector_Distributed<Vector_Dense, T> *solution);

protected:
    template <typename Settings, typename Configuration>
    void analyze(FETI<T> &feti, Settings &settings, Configuration &configuration);

    void empty(FETI<T> &feti);
    void algebraic(FETI<T> &feti, int defect, int sc_size);
    void svd(FETI<T> &feti);

    void set    (FETI<T> &feti, HeatTransferLoadStepConfiguration &configuration);
    void update (FETI<T> &feti, HeatTransferLoadStepConfiguration &configuration);

    void set    (FETI<T> &feti, StructuralMechanicsLoadStepConfiguration &configuration);
    void update (FETI<T> &feti, StructuralMechanicsLoadStepConfiguration &configuration, Vector_Distributed<Vector_Dense, T> *solution);

    bool regMat, R1, R2, onSurface;

    std::vector<std::vector<esint> > fixPoints, fixCols, permutation;
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_H_ */
