
#ifndef SRC_ANALYSIS_LINEARSYSTEM_DIRECTSOLVER_H_
#define SRC_ANALYSIS_LINEARSYSTEM_DIRECTSOLVER_H_

#include "analysis/math/matrix_distributed.h"
#include "analysis/math/vector_distributed.h"
#include "linearsystem.h"


namespace espreso {

template <typename T>
struct DirectLinearSystemSolver: LinearSystemSolver<T> {

    DirectLinearSystemSolver()
    {
        LinearSystemSolver<T>::A = &A;
        LinearSystemSolver<T>::x = &x;
        LinearSystemSolver<T>::b = &b;
        LinearSystemSolver<T>::dirichlet = &dirichlet;
    }

    ~DirectLinearSystemSolver()
    {

    }

protected:
    void setDirichlet();

    Matrix_Distributed<T> A;
    Vector_Distributed<Vector_Dense, T> x, b;
    Vector_Distributed<Vector_Sparse, T> dirichlet;
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_DIRECTSOLVER_H_ */
