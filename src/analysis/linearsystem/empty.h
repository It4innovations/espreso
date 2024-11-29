
#ifndef SRC_ANALYSIS_LINEARSYSTEM_EMPTY_H_
#define SRC_ANALYSIS_LINEARSYSTEM_EMPTY_H_

#include "linearsystem.h"
#include "analysis/pattern/pattern.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/sysutils.h"

namespace espreso {

template <typename T>
struct EmptySystemSolver: LinearSystemSolver<T> {

    Pattern<T>* getPattern(int DOFs)                                                                  { return new PatternUniformDirect<T>(DOFs); }
    Pattern<T>* getPattern(HeatTransferLoadStepConfiguration &configuration       , int multiplicity) { return new PatternUniformDirect<T>(configuration, multiplicity); }
    Pattern<T>* getPattern(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity) { return new PatternUniformDirect<T>(configuration, multiplicity); }

    EmptySystemSolver()
    {
        LinearSystemSolver<T>::A = &A;
        LinearSystemSolver<T>::x = &x;
        LinearSystemSolver<T>::b = &b;
        LinearSystemSolver<T>::dirichlet = &dirichlet;
    }

    ~EmptySystemSolver()
    {

    }

    void set(step::Step &step)
    {

    }

    void update(step::Step &step)
    {
        if (info::ecf->output.print_matrices) {
            eslog::storedata(" STORE: system/{A, b, dirichlet}\n");
            math::store(this->A, utils::filename(utils::debugDirectory(step) + "/system", "A").c_str());
            math::store(this->b, utils::filename(utils::debugDirectory(step) + "/system", "b").c_str());
            math::store(this->dirichlet, utils::filename(utils::debugDirectory(step) + "/system", "dirichlet").c_str());
        }
    }

    bool solve(step::Step &step)
    {
        return true;
    }

    bool postSolve(step::Step &step)
    {
        return true;
    }

    T rhs_without_dirichlet_norm()
    {
        return T{};
    }

protected:
//    void setDirichlet() {}

    Matrix_Distributed<T> A;
    Vector_Distributed<Vector_Dense, T> x, b;
    Vector_Distributed<Vector_Sparse, T> dirichlet;
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_EMPTY_H_ */
