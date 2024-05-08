
#ifndef SRC_ANALYSIS_LINEARSYSTEM_FETISOLVER_H_
#define SRC_ANALYSIS_LINEARSYSTEM_FETISOLVER_H_

#include "analysis/math/matrix_feti.h"
#include "analysis/math/vector_distributed.h"
#include "analysis/math/vector_feti.h"
#include "linearsystem.h"
#include "constrains/constrains.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "feti/feti.h"

namespace espreso {

template <typename T>
struct FETILinearSystemSolver: LinearSystemSolver<T> {

    FETILinearSystemSolver(PhysicsConfiguration &physics, LoadStepSolverConfiguration &loadStep);
    ~FETILinearSystemSolver();

    void set(step::Step &step);
    void update(step::Step &step);
    bool solve(step::Step &step);

private:
    PhysicsConfiguration &physics;
    LoadStepSolverConfiguration &loadStep;

    Matrix_FETI<T> A;
    Vector_FETI<Vector_Dense, T> x, b;
    Vector_Distributed<Vector_Sparse, T> dirichlet;

    Constrains<T> constrains;

    FETI<T> feti;
    bool bem;
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_FETISOLVER_H_ */
