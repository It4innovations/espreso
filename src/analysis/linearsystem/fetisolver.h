
#ifndef SRC_ANALYSIS_LINEARSYSTEM_FETISOLVER_H_
#define SRC_ANALYSIS_LINEARSYSTEM_FETISOLVER_H_

#include "analysis/math/matrix_feti.h"
#include "analysis/math/vector_distributed.h"
#include "analysis/math/vector_feti.h"
#include "analysis/pattern/pattern.h"
#include "linearsystem.h"
#include "constrains/constrains.h"
#include "regularization/regularization.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "feti/feti.h"

namespace espreso {

template <typename T>
struct FETILinearSystemSolver: LinearSystemSolver<T> {

    Pattern<T>* getPattern(int DOFs)                                                                  { return new PatternUniformFETI<T>(DOFs); }
    Pattern<T>* getPattern(HeatTransferLoadStepConfiguration &configuration       , int multiplicity) { return new PatternUniformFETI<T>(configuration, multiplicity); }
    Pattern<T>* getPattern(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity) { return new PatternUniformFETI<T>(configuration, multiplicity); }

    FETILinearSystemSolver(FETIConfiguration &configuration);
    ~FETILinearSystemSolver();

    void set(step::Step &step);
    void update(step::Step &step);
    bool solve(step::Step &step);
    bool postSolve(step::Step &step);

    T rhs_without_dirichlet_norm();

private:
    Matrix_FETI<T> A;
    struct {
        Vector_FETI<Vector_Dense, T> feti;
        Vector_Distributed<Vector_Dense, T> physics;
    } x, b, BtL, ineqBtL;
    Vector_Distributed<Vector_Sparse, T> dirichlet;

    Constrains<T> constrains;
    Regularization<T> regularization;

    FETI<T> feti;
    bool bem;
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_FETISOLVER_H_ */
