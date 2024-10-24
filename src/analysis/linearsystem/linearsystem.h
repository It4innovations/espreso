
#ifndef SRC_ANALYSIS_LINEARSOLVER_LINEARSOLVER_H_
#define SRC_ANALYSIS_LINEARSOLVER_LINEARSOLVER_H_

#include "analysis/math/matrix_base.h"
#include "analysis/math/vector_base.h"
#include "analysis/pattern/pattern.h"

namespace espreso {

namespace step { struct Step; }

template <typename T>
struct LinearSystemSolver {

    virtual Pattern<T>* getPattern(HeatTransferLoadStepConfiguration &configuration       , int multiplicity) =0;
    virtual Pattern<T>* getPattern(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity) =0;

    virtual ~LinearSystemSolver() {}

    virtual void set(step::Step &step) =0;
    virtual void update(step::Step &step) =0;
    virtual bool solve(step::Step &step) =0;

    Pattern<T> *assember, *solver;

    Matrix_Base<T> *A;
    Vector_Distributed<Vector_Dense, T> *x, *b;
    Vector_Distributed<Vector_Sparse, T> *dirichlet;

    virtual T rhs_without_dirichlet_norm() =0;
};

}

#endif /* SRC_ANALYSIS_LINEARSOLVER_LINEARSOLVER_H_ */
