
#ifndef SRC_ANALYSIS_LINEARSYSTEM_SUITESPARSE_H_
#define SRC_ANALYSIS_LINEARSYSTEM_SUITESPARSE_H_

#include "directsolver.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "wrappers/suitesparse/w.suitesparse.direct.h"

namespace espreso {

template <typename T>
struct SuiteSparseLinearSystemSolver: DirectLinearSystemSolver<T> {

    Pattern<T>* getPattern(HeatTransferLoadStepConfiguration &configuration       , int multiplicity) { return new PatternUniformDirect<T>(configuration, multiplicity); }
    Pattern<T>* getPattern(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity) { return new PatternUniformDirect<T>(configuration, multiplicity); }

    SuiteSparseLinearSystemSolver(SuiteSparseConfiguration &configuration): suitesparse(configuration) {}
    ~SuiteSparseLinearSystemSolver() {}

    void set(step::Step &step)
    {
        suitesparse.set(this->A);
    }

    void update(step::Step &step)
    {
        if (info::ecf->output.print_eigen_values) {
            this->A.printEigenValues("A[SOL]", 8);
        }

        if (this->A.updated || this->b.updated || this->dirichlet.updated) {
            this->setDirichlet();
            suitesparse.update(this->A);
            this->A.updated = this->b.updated = true; // set dirichlet update both
        }
        this->dirichlet.updated = false;

        if (info::ecf->output.print_matrices) {
            eslog::storedata(" STORE: system/{A, b, dirichlet}\n");
            math::store(this->A, utils::filename(utils::debugDirectory(step) + "/system", "K").c_str());
            math::store(this->b, utils::filename(utils::debugDirectory(step) + "/system", "b").c_str());
            math::store(this->dirichlet, utils::filename(utils::debugDirectory(step) + "/system", "BC").c_str());
        }
    }

    bool solve(step::Step &step)
    {
        if (suitesparse.solve(this->b, this->x)) {
            this->x.scatter();

            if (false) {
                double x1 = 1e15, x2 = 1 / x1;
                for (int i = 0; i < this->x.cluster.size; ++i) {
                    if (std::fabs(this->x.cluster.vals[i]) < x2) {
                        this->x.cluster.vals[i] = 0;
                    } else {
                        this->x.cluster.vals[i] = std::ceil(x1 * this->x.cluster.vals[i]) * x2;
                    }
                }
            }

            if (info::ecf->output.print_matrices) {
                eslog::storedata(" STORE: system/{x}\n");
                math::store(this->x, utils::filename(utils::debugDirectory(step) + "/system", "x").c_str());
            }
            return true;
        }
        return false;
    }

    T rhs_without_dirichlet_norm()
    {
        return this->b.norm();
    }

private:
    SuiteSparseDirectSolver<T> suitesparse;
};

}


#endif /* SRC_ANALYSIS_LINEARSYSTEM_SUITESPARSE_H_ */
