
#include "dirichlet.implicit.h"
#include "feti/common/applyB.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "basis/utilities/sysutils.h"
#include "basis/utilities/stacktimer.h"

#include <algorithm>

namespace espreso {

template <typename T>
DirichletImplicit<T>::DirichletImplicit(FETI<T> &feti)
: Preconditioner<T>(feti)
{
    Btx.resize(feti.K.size());
    KBtx.resize(feti.K.size());
    indices.resize(feti.K.size());
    permutation.resize(feti.K.size());
    A11.resize(feti.K.size());
    A12.resize(feti.K.size());
    A22.resize(feti.K.size());
    solver.resize(feti.K.size());

    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        Btx[d].resize(feti.K[d].nrows);
        KBtx[d].resize(feti.K[d].nrows);
    }

    eslog::checkpointln("FETI: SET DIRICHLET PRECONDITIONER");
}

template <typename T>
DirichletImplicit<T>::~DirichletImplicit()
{

}

template <typename T>
void DirichletImplicit<T>::info()
{
    if (feti.configuration.exhaustive_info) {
        eslog::info(" = DIRICHLET PRECONDITIONER PROPERTIES                                                       = \n");
        eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
    }
}

template <typename T>
void DirichletImplicit<T>::set(const step::Step &step)
{
stacktimer::enable();
stacktimer::push("DirichletImplicit::set");
    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        permutation[d].resize(feti.B1[d].ncols);
        for (int i = 0; i < feti.B1[d].nnz; ++i) {
            permutation[d][feti.B1[d].cols[i]] = 1;
        }
        indices[d].clear();
        for (size_t i = 0; i < permutation[d].size(); ++i) {
            if (permutation[d][i] == 1) {
                indices[d].push_back(i);
            }
        }
        for (size_t i = 0, j = 0, k = feti.B1[d].ncols - indices[d].size(); i < permutation[d].size(); ++i) {
            permutation[d][i] = (permutation[d][i] == 1 ? k++ : j++);
        }

        Matrix_CSR<T> pK(feti.K[d]);
        math::permute(pK, permutation[d]);

        int size_surface = indices[d].size();
        int size = feti.K[d].nrows;
        int size_A11 = size - size_surface;

        math::spblas::submatrix(pK, A11[d], 0, size_A11, 0, size_A11);

        solver[d].commit(A11[d]);
        solver[d].symbolicFactorization();
    }
stacktimer::pop();
stacktimer::disable();
}

template <typename T>
void DirichletImplicit<T>::update(const step::Step &step)
{
stacktimer::enable();
stacktimer::push("DirichletImplicit::update");
    if (feti.updated.K || feti.updated.B) {
        #pragma omp parallel for
        for (size_t d = 0; d < feti.K.size(); ++d) {
            Matrix_CSR<T> pK(feti.K[d]);
            math::permute(pK, permutation[d]);

            int size_surface = indices[d].size();
            int size = feti.K[d].nrows;
            int size_A11 = size - size_surface;

            math::spblas::submatrix(pK, A11[d],        0, size_A11,        0, size_A11);
            math::spblas::submatrix(pK, A12[d],        0, size_A11, size_A11,     size);
            math::spblas::submatrix(pK, A22[d], size_A11,     size, size_A11,     size);
            solver[d].numericalFactorization();
        }
    }
stacktimer::pop();
stacktimer::disable();
    _print(step);
    eslog::checkpointln("FETI: COMPUTE DIRICHLET PRECONDITIONER");
}

template <typename T>
void DirichletImplicit<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
stacktimer::enable();
stacktimer::push("DirichletImplicit::apply (vector)");
    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        applyBt(feti, d, x, Btx[d]);
        for (size_t i = 0; i < indices[d].size(); ++i) {
            Btx[d].vals[i] = Btx[d].vals[indices[d][i]];
        }

        Vector_Dense<T> vx, vy;
        vx.resize(A12[d].ncols);
        vy.resize(A12[d].nrows);
        math::copy(vx, Btx[d]);

        math::set(vy, T{0});
        math::spblas::apply(vy, T{1}, A12[d], vx);
        solver[d].solve(vy, vx);

        math::set(KBtx[d], T{0});
        math::spblas::applyT(KBtx[d], T{-1}, A12[d], vx);
        for (int r = 0; r < A22[d].nrows; ++r) { // A22 has upper triangle only
            for (int c = A22[d].rows[r]; c < A22[d].rows[r + 1]; ++c) {
                KBtx[d].vals[r] += A22[d].vals[c] * Btx[d].vals[A22[d].cols[c]];
                if (r != A22[d].cols[c]) {
                    KBtx[d].vals[A22[d].cols[c]] += A22[d].vals[c] * Btx[d].vals[r];
                }
            }
        }

        for (size_t i = indices[d].size(); i > 0; --i) {
            KBtx[d].vals[indices[d][i - 1]] = KBtx[d].vals[i - 1];
        }
    }
    applyB(feti, KBtx, y);
    y.synchronize();
stacktimer::pop();
stacktimer::disable();
}

template <typename T>
void DirichletImplicit<T>::_print(const step::Step &step)
{
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: feti/preconditioner/{DirichletImplicit}\n");
        // nothing to store
    }
}

template struct DirichletImplicit<double>;
template struct DirichletImplicit<std::complex<double> >;

}
