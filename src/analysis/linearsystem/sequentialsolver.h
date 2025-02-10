
#ifndef SRC_ANALYSIS_LINEARSYSTEM_SEQUENTIALSOLVER_H_
#define SRC_ANALYSIS_LINEARSYSTEM_SEQUENTIALSOLVER_H_

#include "directsolver.h"
#include "analysis/pattern/pattern.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "mesh/store/nodestore.h"
#include "math/wrappers/math.spsolver.h"

namespace espreso {

template <typename T>
struct SequentialLinearSystemSolver: DirectLinearSystemSolver<T> {

    Pattern<T>* getPattern(int DOFs)                                                                  { return new PatternUniformDirect<T>(DOFs); }
    Pattern<T>* getPattern(HeatTransferLoadStepConfiguration &configuration       , int multiplicity) { return new PatternUniformDirect<T>(configuration, multiplicity); }
    Pattern<T>* getPattern(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity) { return new PatternUniformDirect<T>(configuration, multiplicity); }

    SequentialLinearSystemSolver(SequentialConfiguration &configuration)
    {
        if (info::mpi::size > 1) {
            eslog::globalerror("Sequential solver does not support distributed runs.\n");
        }
    }

    ~SequentialLinearSystemSolver() {}

    bool isSymmetric(Matrix_Type type)
    {
        return type == Matrix_Type::REAL_SYMMETRIC_INDEFINITE
            || type == Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE
            || type == Matrix_Type::COMPLEX_SYMMETRIC
            || type == Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE
            || type == Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE;
    }

    void copyPattern() // shallow copy for asymmetric
    {
        if (isSymmetric(this->A.type)) {
            for (esint i = 0; i < this->A.cluster.nrows; i++) {
                for (esint c = this->A.cluster.rows[i] - Indexing::CSR; c < this->A.cluster.rows[i + 1] - Indexing::CSR; ++c) {
                    if (i <= this->A.cluster.cols[c] - Indexing::CSR) {
                        ++_A.nnz;
                    }
                }
            }
            _A.resize(this->A.cluster.nrows, this->A.cluster.ncols, _A.nnz);
            _A.rows[0] = Indexing::CSR;
            for (esint i = 0, offset = 0; i < this->A.cluster.nrows; i++) {
                for (esint c = this->A.cluster.rows[i] - Indexing::CSR; c < this->A.cluster.rows[i + 1] - Indexing::CSR; ++c) {
                    if (i <= this->A.cluster.cols[c] - Indexing::CSR) {
                        _A.cols[offset++] = this->A.cluster.cols[c];
                    }
                }
                _A.rows[i + 1] = offset + Indexing::CSR;
            }
            _A.shape = Matrix_Shape::UPPER;
        } else {
            _A.shallowCopy(this->A.cluster);
        }
        _A.type = this->A.type;
    }

    void copyValues()
    {
        if (isSymmetric(this->A.type)) {
            for (esint i = 0, offset = 0; i < this->A.cluster.nrows; i++) {
                for (esint c = this->A.cluster.rows[i] - Indexing::CSR; c < this->A.cluster.rows[i + 1] - Indexing::CSR; ++c) {
                    if (i <= this->A.cluster.cols[c] - Indexing::CSR) {
                        _A.vals[offset++] = this->A.cluster.vals[c];
                    }
                }
            }
        }
    }

    void set(step::Step &step)
    {
        copyPattern();
        solver.commit(_A);
        solver.symbolicFactorization();
    }

    void update(step::Step &step)
    {
        if (info::ecf->output.print_eigen_values) {
            this->A.printEigenValues("A[SOL]", 8);
        }

        if (this->A.updated || this->b.updated || this->dirichlet.updated) {
            this->setDirichlet();
            copyValues();
            solver.numericalFactorization();
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
        solver.solve(this->b.cluster, this->x.cluster);
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

    bool postSolve(step::Step &step)
    {
        solver.numericalFactorization();
        solver.solve(this->B.cluster, this->X.cluster);
        return true;
    }

    T rhs_without_dirichlet_norm()
    {
        return this->b.norm();
    }

private:
    DirectSparseSolver<T, esint> solver;
    Matrix_CSR<T, esint> _A;
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_SEQUENTIALSOLVER_H_ */
