
#include "w.mkl.pdss.h"
#include "esinfo/mpiinfo.h"
#include "math/primitives/matrix_csr.h"

#include "wrappers/mpi/communication.h"
#include "esinfo/eslog.hpp"

#include <vector>
#include <complex>

#ifdef HAVE_MKLPDSS
#include "mkl_cluster_sparse_solver.h"

namespace espreso {

template<typename T>
struct MKLPDSSDataHolder {
    void *pt[64];

    esint maxfct;
    esint mnum;

    esint mtype;
    esint phase;
    esint perm;
    esint n;
    esint nrhs;
    esint iparm[64];
    esint msglvl;
    int   comm;
    esint error;

    Matrix_CSR<T> A;
    Vector_Dense<T> b, x;
};

}

#endif

namespace espreso {

template<typename T>
void MKLPDSS<T>::check()
{
#ifndef HAVE_MKLPDSS
    eslog::globalerror("ESPRESO run-time error: cannot call MKL PDSS solver (the library with the solver is not linked).\n");
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
bool MKLPDSS<T>::call(esint phase)
{
#ifdef HAVE_MKLPDSS
    external->phase = phase;
    cluster_sparse_solver(
            external->pt, &external->maxfct, &external->mnum,
            &external->mtype,
            &external->phase,
            &external->n, external->A.vals, external->A.rows, external->A.cols,
            &external->perm, &external->nrhs, external->iparm, &external->msglvl,
            external->b.vals, external->x.vals,
            &external->comm, &external->error);

    switch (external->error) {
    case   0: break;
    case  -1: eslog::error("MKL PDSS: input inconsistent.\n"); break;
    case  -2: eslog::error("MKL PDSS: not enough memory.\n"); break;
    case  -3: eslog::error("MKL PDSS: reordering problem.\n"); break;
    case  -4: eslog::error("MKL PDSS: zero pivot, numerical factorization or iterative refinement problem.\n"); break;
    case  -5: eslog::error("MKL PDSS: unclassified (internal) error.\n"); break;
    case  -6: eslog::error("MKL PDSS: reordering failed.\n"); break;
    case  -7: eslog::error("MKL PDSS: diagonal matrix is singular.\n"); break;
    case  -8: eslog::error("MKL PDSS: 32-bit integer overflow problem.\n"); break;
    case  -9: eslog::error("MKL PDSS: not enough memory for OOC.\n"); break;
    case -10: eslog::error("MKL PDSS: error opening OOC files.\n"); break;
    case -11: eslog::error("MKL PDSS: read/write error with OOC files.\n"); break;
    }
    return external->error == 0;
#endif
    return false;
}

template<typename T>
static void _info(const Matrix_Distributed<T> &A)
{
    eslog::info(" = LINEAR SOLVER :: MKL                                TYPE :: PARALLEL DIRECT SPARSE SOLVER = \n");
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

template<typename T>
bool MKLPDSS<T>::set(const Matrix_Distributed<T> &A)
{
#ifdef HAVE_MKLPDSS
    _info(A); // print the info before call the solver
    double start = eslog::time();

    external = new MKLPDSSDataHolder<T>();
    external->n = A.decomposition->totalSize;

    external->maxfct = 1; // dummy
    external->mnum = 1; // dummy

    external->nrhs = 1;

    std::fill(external->iparm, external->iparm + 64, 0);
    external->iparm[0] = 1; // Use filled values.
    // Fill-in reducing ordering for the input matrix.
    external->iparm[1] = info::mpi::size > 1 ? 10 : 3; // MPI or parallel
    // Matrix input format.
    external->iparm[34] = 1 - Indexing::CSR; // 0- or 1- based indexing
    external->iparm[39] = 2; // distributed A, x, rhs
    external->iparm[40] = A.decomposition->begin + Indexing::CSR;
    external->iparm[41] = A.decomposition->end - 1 + Indexing::CSR;

    external->msglvl = 0;
    external->comm = MPI_Comm_c2f(info::mpi::comm);

    switch (A.cluster.type) {
    case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:    external->mtype =  2; break;
    case Matrix_Type::REAL_SYMMETRIC_INDEFINITE:           external->mtype = -2; break;
    case Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC:         external->mtype =  1; break;
    case Matrix_Type::REAL_NONSYMMETRIC:                   external->mtype = 11; break;
    case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE: external->mtype =  4; break;
    case Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE:        external->mtype = -4; break;
    case Matrix_Type::COMPLEX_SYMMETRIC:                   external->mtype =  6; break;
    case Matrix_Type::COMPLEX_STRUCTURALLY_SYMMETRIC:      external->mtype =  3; break;
    case Matrix_Type::COMPLEX_NONSYMMETRIC:                external->mtype = 13; break;
    }

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
        external->A.rows[0] = Indexing::CSR;
        for (esint i = nhalo, offset = 0; i < A.cluster.nrows; i++) {
            for (esint c = A.cluster.rows[i] - Indexing::CSR; c < A.cluster.rows[i + 1] - Indexing::CSR; ++c) {
                if (A.decomposition->begin + i - nhalo <= A.cluster.cols[c] - Indexing::CSR) {
                    external->A.cols[offset++] = A.cluster.cols[c];
                }
            }
            external->A.rows[i - nhalo + 1] = offset + Indexing::CSR;
        }
    } else {
        esint nhalo = A.decomposition->halo.size();
        for (esint i = nhalo; i < A.cluster.nrows; i++) {
            for (esint c = A.cluster.rows[i] - Indexing::CSR; c < A.cluster.rows[i + 1] - Indexing::CSR; ++c) {
                ++external->A.nnz;
            }
        }
        external->A.resize(A.cluster.nrows - nhalo, A.cluster.ncols, external->A.nnz);
        external->A.rows[0] = Indexing::CSR;
        for (esint i = nhalo, offset = 0; i < A.cluster.nrows; i++) {
            for (esint c = A.cluster.rows[i] - Indexing::CSR; c < A.cluster.rows[i + 1] - Indexing::CSR; ++c) {
                external->A.cols[offset++] = A.cluster.cols[c];
            }
            external->A.rows[i - nhalo + 1] = offset + Indexing::CSR;
        }
    }
    bool status = call(11);
    eslog::solver(" = PREPARE PERSISTENT DATA (SYMBOLIC FACTORIZATION)                               %8.3f s = \n", eslog::time() - start);
    return status;
#endif
    return false;
}

template<typename T>
bool MKLPDSS<T>::update(const Matrix_Distributed<T> &A)
{
#ifdef HAVE_MKLPDSS
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
    bool status = call(22);
    eslog::solver("       - NUMERICAL FACTORIZATION                                            %8.3f s -  \n", eslog::time() - start);
    return status;
#endif
    return false;
}

template<typename T>
bool MKLPDSS<T>::solve(const Vector_Distributed<Vector_Dense, T> &b, Vector_Distributed<Vector_Dense, T> &x)
{
#ifdef HAVE_MKLPDSS
    external->b.vals = b.cluster.vals + b.decomposition->halo.size();
    external->x.vals = x.cluster.vals + x.decomposition->halo.size();
    double start = eslog::time();

    bool status = call(33); // solve at once
    eslog::solver("       - SOLVER TIME                                                        %8.3f s -  \n", eslog::time() - start);
    return status;
#endif
    return false;
}

template<typename T>
bool MKLPDSS<T>::solve(const Vector_Distributed<Matrix_Dense, T> &B, Vector_Distributed<Matrix_Dense, T> &X)
{
#ifdef HAVE_MKLPDSS
    bool status = true;
    double start = eslog::time();
    for (int r = 0; r < B.cluster.nrows; ++r) {
        external->b.vals = B.cluster.vals + B.decomposition->halo.size() + r * B.cluster.ncols;
        external->x.vals = X.cluster.vals + X.decomposition->halo.size() + r * B.cluster.ncols;

        status |= call(33); // solve at once
    }
    eslog::solver("       - SOLVER TIME                                                        %8.3f s -  \n", eslog::time() - start);
    return status;
#endif
    return false;
}

template<typename T>
void MKLPDSS<T>::clear()
{
#ifdef HAVE_MKLPDSS
    if (external) delete external;
#endif
}

template struct MKLPDSS<double>;
template struct MKLPDSS<std::complex<double> >;

}
