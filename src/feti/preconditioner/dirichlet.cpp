
#include "dirichlet.h"
#include "feti/common/applyB.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "basis/utilities/sysutils.h"

#include "math/operations/schur_csx_dny.h"

#include <algorithm>

namespace espreso {

template <typename T>
Dirichlet<T>::Dirichlet(FETI<T> &feti)
: Preconditioner<T>(feti)
{
    Btx.resize(feti.K.size());
    KBtx.resize(feti.K.size());
    sc.resize(feti.K.size());
    indices.resize(feti.K.size());

    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        Btx[d].resize(feti.K[d].nrows);
        KBtx[d].resize(feti.K[d].nrows);
        sc[d].shape = Matrix_Shape::UPPER;
    }

    eslog::checkpointln("FETI: SET DIRICHLET PRECONDITIONER");
}

template <typename T>
Dirichlet<T>::~Dirichlet()
{

}

template <typename T>
void Dirichlet<T>::info()
{
    if (feti.configuration.exhaustive_info) {
        eslog::info(" = DIRICHLET PRECONDITIONER PROPERTIES                                                       = \n");
        eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
    }
}

template <typename T>
void Dirichlet<T>::update(const step::Step &step)
{
    if (feti.updated.K || feti.updated.B) {
        #pragma omp parallel for
        for (size_t d = 0; d < feti.K.size(); ++d) {
            std::vector<int> permutation(feti.B1[d].ncols);
            for (int i = 0; i < feti.B1[d].nnz; ++i) {
                permutation[feti.B1[d].cols[i]] = 1;
            }
            indices[d].clear();
            for (size_t i = 0; i < permutation.size(); ++i) {
                if (permutation[i] == 1) {
                    indices[d].push_back(i);
                }
            }
            if (DirectSparseSolver<T>::provideSC()) {
                DirectSparseSolver<T> Ksolver;
                Ksolver.commit(feti.K[d]);
                sc[d].resize(indices[d].size(), indices[d].size());
                Ksolver.getSC(sc[d], permutation);
            } else {
                for (size_t i = 0, j = 0, k = feti.B1[d].ncols - indices[d].size(); i < permutation.size(); ++i) {
                    permutation[i] = (permutation[i] == 1 ? k++ : j++);
                }
                sc[d].resize(indices[d].size(), indices[d].size());
                _manual(feti.K[d], sc[d], permutation);
            }
        }
    }
    _print(step);
    eslog::checkpointln("FETI: COMPUTE DIRICHLET PRECONDITIONER");
}

template <typename T>
void Dirichlet<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        applyBt(feti, d, x, Btx[d]);
        for (size_t i = 0; i < indices[d].size(); ++i) {
            Btx[d].vals[i] = Btx[d].vals[indices[d][i]];
        }
        math::blas::apply(KBtx[d], T{1}, sc[d], T{0}, Btx[d]);
        for (size_t i = indices[d].size(); i > 0; --i) {
            KBtx[d].vals[indices[d][i - 1]] = KBtx[d].vals[i - 1];
        }
    }
    applyB(feti, KBtx, y);
    y.synchronize();
}

template <typename T>
void Dirichlet<T>::_manual(Matrix_CSR<T> &K, Matrix_Dense<T> &sc, std::vector<int> &permutation)
{
    Matrix_CSR<T> pK(K);
    math::permute(pK, permutation);

    int size_sc = sc.nrows;

    MatrixCsxView_new<T,int> pK_new = MatrixCsxView_new<T,int>::from_old(pK);

    MatrixDenseData_new<T> sc_new;
    sc_new.set(sc.nrows, sc.ncols, 'R', AllocatorCPU_new::get_singleton());
    sc_new.prop = pK_new.prop;
    sc_new.alloc();

    using op_schur_t = math::operations::schur_csx_dny<T,int>;
    std::unique_ptr<op_schur_t> op_sc = op_schur_t::make(op_schur_t::implementation_selector::sparse_solver);
    op_sc->set_matrix(&pK_new, size_sc);
    op_sc->set_sc(&sc_new);
    op_sc->set_need_solve_A11(false);
    op_sc->preprocess();
    op_sc->perform();

    size_t idx = 0;
    for(size_t r = 0; r < sc_new.nrows; r++) {
        size_t start = ((sc.shape == Matrix_Shape::UPPER) ?            r :   0);
        size_t end   = ((sc.shape == Matrix_Shape::UPPER) ? sc_new.ncols : r+1);
        for(size_t c = start; c < end; c++) {
            sc.vals[idx] = sc_new.vals[r * sc_new.ld + c];
            idx++;
        }
    }

    // Matrix_CSR<T> A11_sp;
    // Matrix_CSR<T> A21_sp; // = A12c_sp
    // Matrix_Dense<T> A22t_dn;
    // Matrix_Dense<T> A12t_dn;
    // Matrix_Dense<T> A11iA12_dn;
    // math::spblas::submatrix(pK, A11_sp ,        0, size_A11,        0, size_A11);
    // math::spblas::submatrix(pK, A21_sp,         0, size_A11, size_A11,     size, false,  true); // = A12c_sp
    // math::spblas::submatrix(pK, A22t_dn, size_A11,     size, size_A11,     size,  true, false, true);
    // math::spblas::submatrix(pK, A12t_dn,        0, size_A11, size_A11,     size,  true, false, true);

    // DirectSparseSolver<T, int> solver;
    // solver.commit(A11_sp);
    // solver.symbolicFactorization();
    // solver.numericalFactorization();
    // solver.solve(A12t_dn, A11iA12_dn);

    // SpBLAS<Matrix_CSR, T> A21(A21_sp, true);
    // A21.apply(A22t_dn, T{-1}, T{1}, A11iA12_dn, true);

    // sc.shape = K.shape;
    // sc.resize(A22t_dn);
    // for(int r = 0, i = 0; r < sc.nrows; ++r) {
    //     for(int c = sc.shape == Matrix_Shape::FULL ? 0 : r; c < sc.ncols; ++c, ++i) {
    //         sc.vals[i] = A22t_dn.vals[r * sc.ncols + c];
    //     }
    // }
}

template <typename T>
void Dirichlet<T>::_print(const step::Step &step)
{
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: feti/preconditioner/{Dirichlet}\n");
        for (size_t d = 0; d < feti.K.size(); ++d) {
            math::store(sc[d], utils::filename(utils::debugDirectory(step) + "/feti/precondition", (std::string("dirichlet") + std::to_string(d)).c_str()).c_str());
        }
    }
}

template struct Dirichlet<double>;
template struct Dirichlet<std::complex<double> >;

}
