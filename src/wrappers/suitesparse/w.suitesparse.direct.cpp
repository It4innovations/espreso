
#include "w.suitesparse.direct.h"
#include "esinfo/mpiinfo.h"
#include "math/primitives/matrix_csr.h"

#include "wrappers/mpi/communication.h"
#include "esinfo/eslog.hpp"

#include <vector>
#include <complex>

#ifdef HAVE_SUITESPARSE

#include "w.suitesparse.cholmod.h"

namespace espreso {

template<typename T>
struct SuiteSparseDataHolder {
    cholmod_common cm_common;
    cholmod_sparse * cm_matrix = nullptr;
    cholmod_factor * cm_factor = nullptr;
    int stage = 0; // 0 = completely uninitialized, 1 = symbolic factorization done, 2 = numeric factorization done

    Matrix_CSR<T> A;
};

}

#endif

namespace espreso {

template <typename T>
void SuiteSparseDirectSolver<T>::check()
{
#ifndef HAVE_SUITESPARSE
    eslog::globalerror("ESPRESO run-time error: cannot call SuiteSparse solver (the library with the solver is not linked).\n");
#else
    if (info::mpi::size > 1) {
        eslog::globalerror("ESPRESO run-time error: SuiteSparse does not support MPI parallelization.\n");
    }
#endif
}

bool _isSymmetric(Matrix_Type type)
{
    return type == Matrix_Type::REAL_SYMMETRIC_INDEFINITE
        || type == Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE
        || type == Matrix_Type::COMPLEX_SYMMETRIC
        || type == Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE
        || type == Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE;
}

template<typename T>
static void _info(const Matrix_Distributed<T> &A)
{
    eslog::info(" = LINEAR SOLVER :: SUITE SPARSE                                                             = \n");
    switch (A.cluster.type) {
    case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:
        eslog::info(" = MATRIX TYPE ::                                           REAL SYMMETRIC POSITIVE DEFINITE = \n");
        break;
    case Matrix_Type::REAL_SYMMETRIC_INDEFINITE:
        eslog::info(" = MATRIX TYPE ::                                                  REAL SYMMETRIC INDEFINITE = \n");
        break;
    case Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC:
        eslog::info(" = MATRIX TYPE ::                                  REAL NONSYMMETRIC (STRUCTURALY SYMMETRIC) = \n");
        break;
    case Matrix_Type::REAL_NONSYMMETRIC:
        eslog::info(" = MATRIX TYPE ::                                                          REAL NONSYMMETRIC = \n");
        break;

    case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE:
        eslog::info(" = MATRIX TYPE ::                                        COMPLEX_HERMITIAN_POSITIVE_DEFINITE = \n");
        break;
    case Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE:
        eslog::info(" = MATRIX TYPE ::                                               COMPLEX_HERMITIAN_INDEFINITE = \n");
        break;
    case Matrix_Type::COMPLEX_SYMMETRIC:
        eslog::info(" = MATRIX TYPE ::                                                          COMPLEX_SYMMETRIC = \n");
        break;
    case Matrix_Type::COMPLEX_STRUCTURALLY_SYMMETRIC:
        eslog::info(" = MATRIX TYPE ::                                             COMPLEX_STRUCTURALLY_SYMMETRIC = \n");
        break;
    case Matrix_Type::COMPLEX_NONSYMMETRIC:
        eslog::info(" = MATRIX TYPE ::                                                       COMPLEX_NONSYMMETRIC = \n");
        break;
    }
}

template <typename T>
bool SuiteSparseDirectSolver<T>::set(const Matrix_Distributed<T> &A)
{
#ifdef HAVE_SUITESPARSE
    if(A.cluster.nrows != A.cluster.ncols) {
        eslog::error("commit: matrix has to be square\n");
    }
    if(
            A.type != Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE &&
            A.type != Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE) {
        eslog::error("SuiteSparse: not implemented support for non SPD or HPD matrices\n");
    }
    _info(A); // print the info before call the solver
    double start = eslog::time();

    external = new SuiteSparseDataHolder<T>();

    // pick only upper triangle (since composer does not set correct dirichlet in symmetric matrices)
    if (_isSymmetric(A.cluster.type)) {
        esint nhalo = A.decomposition->halo.size();
        for (esint i = nhalo; i < A.cluster.nrows; i++) {
            for (esint c = A.cluster.rows[i] - Indexing::CSR; c < A.cluster.rows[i + 1] - Indexing::CSR; ++c) {
                if (A.decomposition->begin + i - nhalo <= A.cluster.cols[c] - Indexing::CSR) {
                    ++external->A.nnz;
                }
            }
        }
        external->A.resize(A.cluster.nrows - nhalo, A.cluster.ncols, external->A.nnz);
        external->A.rows[0] = 0;
        for (esint i = nhalo, offset = 0; i < A.cluster.nrows; i++) {
            for (esint c = A.cluster.rows[i] - Indexing::CSR; c < A.cluster.rows[i + 1] - Indexing::CSR; ++c) {
                if (A.decomposition->begin + i - nhalo <= A.cluster.cols[c] - Indexing::CSR) {
                    external->A.cols[offset++] = A.cluster.cols[c] - Indexing::CSR;
                }
            }
            external->A.rows[i - nhalo + 1] = offset;
        }
    } else {
        esint nhalo = A.decomposition->halo.size();
        for (esint i = nhalo; i < A.cluster.nrows; i++) {
            for (esint c = A.cluster.rows[i] - Indexing::CSR; c < A.cluster.rows[i + 1] - Indexing::CSR; ++c) {
                ++external->A.nnz;
            }
        }
        external->A.resize(A.cluster.nrows - nhalo, A.cluster.ncols, external->A.nnz);
        external->A.rows[0] = 0;
        for (esint i = nhalo, offset = 0; i < A.cluster.nrows; i++) {
            for (esint c = A.cluster.rows[i] - Indexing::CSR; c < A.cluster.rows[i + 1] - Indexing::CSR; ++c) {
                external->A.cols[offset++] = A.cluster.cols[c] - Indexing::CSR;
            }
            external->A.rows[i - nhalo + 1] = offset;
        }
    }

    if (external->cm_matrix == nullptr) {
        external->cm_matrix = new cholmod_sparse();
    } else {
        delete external->cm_matrix;
        external->cm_matrix = new cholmod_sparse();
    }

    external->cm_matrix->nrow = external->A.ncols;
    external->cm_matrix->ncol = external->A.nrows;
    external->cm_matrix->nzmax = external->A.nnz;
    external->cm_matrix->nz = nullptr;
    external->cm_matrix->z = nullptr;
    external->cm_matrix->stype = -1; // UPPER in CSR, but LOWER in CSC
    external->cm_matrix->itype = _getCholmodItype<int>();
    external->cm_matrix->xtype = _getCholmodXtype<T>();
    external->cm_matrix->dtype = _getCholmodDtype<T>();
    external->cm_matrix->sorted = 1;
    external->cm_matrix->packed = 1;

    external->cm_matrix->p = external->A.rows;
    external->cm_matrix->i = external->A.cols;
    external->cm_matrix->x = external->A.vals;

    external->cm_factor = _analyze<int>(external->cm_matrix, external->cm_common);

    external->stage = 1;
    eslog::solver(" = PREPARE PERSISTENT DATA (SYMBOLIC FACTORIZATION)                               %8.3f s = \n", eslog::time() - start);
    return true;
#endif
    return false;
}

template <typename T>
bool SuiteSparseDirectSolver<T>::update(const Matrix_Distributed<T> &A)
{
#ifdef HAVE_SUITESPARSE
    double start = eslog::time();
    if (_isSymmetric(A.cluster.type)) {
        esint nhalo = A.decomposition->halo.size();
        for (esint i = nhalo, offset = 0; i < A.cluster.nrows; i++) {
            for (esint c = A.cluster.rows[i] - Indexing::CSR; c < A.cluster.rows[i + 1] - Indexing::CSR; ++c) {
                if (A.decomposition->begin + i - nhalo <= A.cluster.cols[c] - Indexing::CSR) {
                    external->A.vals[offset++] = A.cluster.vals[c];
                }
            }
        }
    } else {
        esint nhalo = A.decomposition->halo.size();
        for (esint i = nhalo, offset = 0; i < A.cluster.nrows; i++) {
            for (esint c = A.cluster.rows[i] - Indexing::CSR; c < A.cluster.rows[i + 1] - Indexing::CSR; ++c) {
                external->A.vals[offset++] = A.cluster.vals[c];
            }
        }
    }
    if (external->stage < 1) eslog::error("update: invalid order of operations in SuiteSparse solver\n");

    _factorize<int>(external->cm_factor, external->cm_matrix, external->cm_common);

    eslog::solver("       - NUMERICAL FACTORIZATION                                            %8.3f s -  \n", eslog::time() - start);
    external->stage = 2;
    return true;
#endif
    return false;
}

template <typename T>
bool SuiteSparseDirectSolver<T>::solve(const Vector_Distributed<Vector_Dense, T> &b, Vector_Distributed<Vector_Dense, T> &x)
{
#ifdef HAVE_SUITESPARSE
    double start = eslog::time();

    if(external->stage < 2) eslog::error("solve: invalid order of operations in SuiteSparse solver\n");

    cholmod_dense cm_rhs;
    cm_rhs.nrow = cm_rhs.d = cm_rhs.nzmax = b.cluster.size;
    cm_rhs.ncol = 1;
    cm_rhs.x = b.cluster.vals + b.decomposition->halo.size();
    cm_rhs.xtype = _getCholmodXtype<T>();
    cm_rhs.dtype = _getCholmodDtype<T>();

    cholmod_dense * cm_sol = _solve<int>(CHOLMOD_A, external->cm_factor, &cm_rhs, external->cm_common);

    std::copy_n(reinterpret_cast<T*>(cm_sol->x), cm_sol->nrow, x.cluster.vals + x.decomposition->halo.size());

    _free<int>(cm_sol, external->cm_common);
    eslog::solver("       - SOLVER TIME                                                        %8.3f s -  \n", eslog::time() - start);
    return true;
#endif
    return false;
}

template <typename T>
void SuiteSparseDirectSolver<T>::clear()
{
#ifdef HAVE_SUITESPARSE
    delete external;
#endif
}

template class SuiteSparseDirectSolver<double>;

}
