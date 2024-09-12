
#include "math.h"
#include "math/primitives/matrix_info.h"

#include <algorithm>
#include <numeric>
#include <random>

namespace espreso {
namespace math {

template <typename T, typename I> void _orthonormalize(Matrix_Dense<T, I> &A)
{
    // Gram-Schmidt
//    Vector_Dense<T, I> mlt; mlt.resize(A.nrows);
//    Vector_Dense<T, I> row; row.size = A.ncols;
//    Matrix_Dense<T, I> mat; mat.ncols = A.ncols; mat.vals = A.vals;
//    for (esint r = 0; r < A.nrows; ++r) {
//        row.vals  = A.vals + r * A.ncols;
//        mat.nrows = r;
//        mlt.size  = r;
//        math::blas::multiply(T{ 1}, mat, row, T{0}, mlt);
//        math::blas::multiply(T{-1}, mat, mlt, T{1}, row, true);
//        math::blas::scale(row.size, T{1} / math::blas::norm(row.size, row.vals, 1), row.vals, 1);
//    }

    // MGS
    for (esint r = 0; r < A.nrows; ++r) {
        for (esint rr = 0; rr < r; ++rr) {
            T scale = math::blas::dot(A.ncols, A.vals + rr * A.ncols, 1, A.vals + r * A.ncols, 1);
            math::blas::add(A.ncols, A.vals + r * A.ncols, 1, -scale, A.vals + rr * A.ncols, 1);
        }
        math::blas::scale(A.ncols, T{1.} / math::blas::norm(A.ncols, A.vals + r * A.ncols, 1), A.vals + r * A.ncols, 1);
    }
}

template <typename T, typename I> void _permute(Matrix_CSR<T, I> &A, const std::vector<I> &perm)
{
    if (A.nrows != A.ncols) { eslog::error("cannot permute non-square matrix.\n"); }

    struct __ijv__ {
        T v; I r; I c;
        __ijv__(I r, I c, T v): v(v), r(r), c(c) {}
        bool operator<(const __ijv__ &other) const { if (r == other.r) return c < other.c; else return r < other.r; }
    };
    std::vector<__ijv__> ijv; ijv.reserve(A.ncols);

    for (I r = 0; r < A.nrows; ++r) {
        for (I c = A.rows[r]; c < A.rows[r + 1]; ++c) {
            I pr = std::min(perm[r], perm[A.cols[c - Indexing::CSR] - Indexing::CSR]);
            I pc = std::max(perm[r], perm[A.cols[c - Indexing::CSR] - Indexing::CSR]);
            ijv.push_back({pr, pc, A.vals[c - Indexing::CSR]});
        }
    }
    std::sort(ijv.begin(), ijv.end());

    A.rows[0] = Indexing::CSR;
    for (I i = 0, r = A.nrows; i < A.nnz; ++i) {
        if (r != ijv[i].r) A.rows[ijv[i].r + 1] = A.rows[ijv[i].r];
        r = ijv[i].r;
        A.rows[r + 1]++;
        A.cols[i] = ijv[i].c + Indexing::CSR;
        A.vals[i] = ijv[i].v;
    }
}

template <typename T, typename I> void _getNullPivots(Matrix_Dense<T, I> &R, std::vector<I> &pivots)
{
    std::vector<T> N(R.vals, R.vals + R.nnz);
    std::vector<I> piv(R.ncols); std::iota(piv.begin(), piv.end(), 0);

    for (I p = 0; p < R.nrows; ++p) {
        I index = std::max_element(N.begin(), N.end() - p * R.ncols, [] (const double &i, const double &j) { return std::fabs(i) <= std::fabs(j); }) - N.begin();
        I col = index % R.ncols;
        I row = index / R.ncols;
        // swap pivot to the last row / col
        std::swap(piv[col], piv[R.ncols - p - 1]);
        for (I k = 0; k < R.ncols    ; ++k) { std::swap(N[row * R.ncols + k], N[(R.nrows - p - 1) * R.ncols + k]); }
        for (I k = 0; k < R.nrows - p; ++k) { std::swap(N[k * R.ncols + col], N[k * R.ncols + R.ncols - p - 1]); }
        T pivot = N[(R.nrows - p - 1) * R.ncols + R.ncols - p - 1];
        for (I r = 0; r < R.nrows - p - 1; ++r) {
            for (I c = 0; c < R.ncols - p; ++c) {
                N[r * R.ncols + c] -= N[r * R.ncols + R.ncols - p - 1] * N[(R.nrows - p - 1) * R.ncols + c] / pivot;
            }
        }
    }
    for (I p = 0; p < R.nrows; ++p) {
        pivots.push_back(piv[piv.size() - p - 1]);
    }
    std::sort(pivots.begin(), pivots.end());
}


template <typename T, typename I> void _getKernel(Matrix_CSR<T, I> &A, Matrix_Dense<T, I> &R, Matrix_CSR<T, I> &regMat, I maxDefect, I scSize)
{
//
// Routine calculates kernel Kplus_R of K satisfied euqality K * Kplus_R = O,
// where O is zero matrix, and it makes the matrix K non-singular (A_reg)
// utilizing spectral conditions of Schur complement. Then ||K-K*inv(A_reg)*K||=0.0
//
//    1) diagonalScaling
//  reducing of big jump coefficient effect (TODO include diagonal scaling into whole ESPRESO)
    bool diagonalScaling                                = true;

//    2) permutVectorActive
//  random selection of singular DOFs
// 0 - no permut., 1 - std::vector shuffle
    esint permutVectorActive                            = 0;

//    3) use_null_pivots_or_s_set
// NtN_Mat from null pivots or fixing DOFs
//    bool use_null_pivots_or_s_set                       = true;

//    4) diagonalRegularization
//  regularization only on diagonal elements (big advantage: patern of K and A_regular is the same !!!)
//  size of set 's' = defect(K)
//  It's is active, only if and only if 'use_null_pivots_or_s_set = true'
    bool diagonalRegularization                         = true;

//    5) get_n_first_and_n_last_eigenvals_from_dense_A
// get and print 2*n K eigenvalues (K is temporarily converted to dense);
//    esint get_n_first_and_n_last_eigenvals_from_dense_A = 10;

//    6) get_n_first_and_n_last_eigenvals_from_dense_S
// get and print 2*n S eigenvalues
//    esint get_n_first_and_n_last_eigenvals_from_dense_S = 10;

//    7) plot_n_first_n_last_eigenvalues
// get of K eigenvalues (K is temporarily converted to dense matrix);
//    esint plot_n_first_n_last_eigenvalues               = 0;

//    8) fixing_nodes_or_dof
// non-singular part determined by fixing nodes (FN),
// min(fixing_nodes_or_dof)>=3; if variable is nonzero,
// parameter sc_size is set to fixing_nodes_or_dof*dofPerNode
//    esint fixing_nodes_or_dof                           = 0;
//    esint dofPerNode                                    = 3;
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    1) cond_numb_for_singular_matrix
//  If cond(K) > cond_numb_for_singular_matrix, K is considered as singular matrix.
//    double cond_numb_for_singular_matrix                = 1e13;

//    2) check_nonsing
// if check_nonsing>0, checking of A_rr non-singularity is activated and it is repeated
// (check_nonsing) times.
//    esint check_nonsing                                 = 0;

//    3) max_size_of_dense_matrix_to_get_eigs
// if size of K is less then CHECA_N..., K is converted to dense format to get eigenvalues.
//    esint max_size_of_dense_matrix_to_get_eigs          = 2500;

//    4) sc_size
// specification of size of Schur complement used for detection of zero eigenvalues.
//esint  sc_size >= expected defect 'd' (e.g. in elasticity d=6).
//    esint sc_size                                       = scSize;

//    5) twenty
// testing last twenty eigenvalues of S to distinguish, if d-last ones are zero or not.
//    esint twenty                                        = 20;
// twenty eigenvalues are ascstd::endly ordered in d = d[0],d[1], ..., d[n-2],d[n-1]

//    6) jump_in_eigenvalues_alerting_singularity
// if d[i]/d[i+1]< jump_in_eigenvalues_alerting_singularity, d[i] is last nonzero eigenvalue
    double jump_in_eigenvalues_alerting_singularity     = 1.0e-5;

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

    double rho = A.vals[0]; // max value on diagonal
    for (esint r = 0; r < A.nrows; ++r) {
        rho = std::max(rho, A.vals[A.rows[r] - Indexing::CSR]);
    }

    Matrix_CSR<T, I> _A(A);

    // scaling
    double di = 1, dj = 1;
    for (esint r = 0; r < _A.nrows; r++) {
        if (diagonalScaling) { di = _A.vals[_A.rows[r] - Indexing::CSR]; }
        for (esint c = _A.rows[r]; c < _A.rows[r + 1]; c++){
            if (diagonalScaling) {
                dj = _A.vals[_A.rows[_A.cols[c - Indexing::CSR] - Indexing::CSR] - Indexing::CSR];
            }
            _A.vals[c - Indexing::CSR] = A.vals[c - Indexing::CSR] / sqrt(di * dj);
        }
    }

    if (A.nrows < scSize){
        scSize = A.nrows;
    }

    esint nonsing_size = A.nrows - scSize;
    std::vector<I> permVec(A.nrows);

    // permutation
    std::iota(permVec.begin(), permVec.end(), 0);
    if (permutVectorActive) {
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(permVec.begin(), permVec.end(), g);
        std::sort(permVec.begin(), permVec.begin() + nonsing_size);
        std::sort(permVec.begin() + nonsing_size, permVec.end());
        std::vector<I> reversePerm(A.nrows); std::iota(reversePerm.begin(), reversePerm.end(), 0);
        std::sort(reversePerm.begin(), reversePerm.end(), [&] (int i, int j) { return permVec[i] < permVec[j]; });
        math::permute(_A, reversePerm);
    }

    // | A_rr  A_rs |
    // |       A_ss |
    Matrix_CSR<T, I> A_rr;
    Matrix_Dense<T, I> A_rs, A_ss;
    math::spblas::submatrix(_A, A_rr, 0, nonsing_size, 0, nonsing_size);
    math::spblas::submatrix(_A, A_rs, 0, nonsing_size, nonsing_size, nonsing_size + scSize, true); // transpose == get COL_MAJOR format required for the solver
    math::spblas::submatrix(_A, A_ss, nonsing_size, nonsing_size + scSize, nonsing_size, nonsing_size + scSize);

    DirectSparseSolver<T, I> A_rr_solver(A_rr);
    if (nonsing_size) {
        Matrix_Dense<T, I> invKrrKrs, KsrInvKrrKrs;
        A_rr_solver.symbolicFactorization();
        A_rr_solver.numericalFactorization();
        A_rr_solver.solve(A_rs, invKrrKrs);
        KsrInvKrrKrs.resize(invKrrKrs.nrows, invKrrKrs.nrows);
        math::blas::multiply(T{1}, A_rs, invKrrKrs, T{0}, KsrInvKrrKrs, false, true);
        // A_ss -= KsrInvKrrKrs (A_ss is UPPER)
        for (I r = 0, i = 0; r < A_ss.nrows; ++r) {
            for (I c = r; c < A_ss.ncols; ++c, ++i) {
                A_ss.vals[i] -= KsrInvKrrKrs.vals[r * KsrInvKrrKrs.ncols + c];
            }
        }
    }

    Vector_Dense<T, I> eigval; eigval.resize(A_ss.nrows);
    Matrix_Dense<T, I> eigvec; eigvec.resize(A_ss.nrows, A_ss.nrows);
//    double tt = eslog::time();
//    math::lapack::get_eig_sym(A_ss, eigval, eigvec);
    math::lapack::get_eig_sym(A_ss, eigval, eigvec, 1, maxDefect + 1);
//    printf("TIME %f\n", eslog::time() - tt);
//
//    for (esint i = 0; i < std::min(maxDefect + 1, eigval.size); ++i) {
//       printf("%+e\n", eigval.vals[i]);
//    }

    // identification of defect in K
    esint defect = std::min(maxDefect + 1, A_ss.nrows) - 1;
    while (defect && std::fabs(eigval.vals[defect - 1] / eigval.vals[defect]) > jump_in_eigenvalues_alerting_singularity) {
        --defect;
    }

    if (defect == 0) {
        return;
    }

    // CREATING KERNEL R_s FOR SINGULAR PART (SCHUR COMPLEMENT)
    Matrix_Dense<T, I> R_s;
    R_s.resize(defect, A_ss.ncols);
    for (I r = 0; r < R_s.nrows; ++r) {
        for (I c = 0; c < R_s.ncols; ++c) {
            R_s.vals[r * R_s.ncols + c] = eigvec.vals[c * eigvec.ncols + r];
        }
    }

    // CREATING KERNEL R_r FOR NON-SINGULAR PART
    Matrix_Dense<T, I> R_r;
    if (A_rr.ncols) {
        Matrix_Dense<T, I> X; X.resize(R_s.nrows, A_rs.ncols); // get X in COL_MAJOR
        math::blas::multiply(T{1}, R_s, A_rs, T{0}, X);
        A_rr_solver.solve(X, R_r); // inv(A_rr)*A_rs*R_s
    }

    // CREATING WHOLE KERNEL Kplus_R = [ (R_r)^T (R_s)^T ]^T
    R.resize(R_s.nrows, A.ncols);
    for (I r = 0; r < R.nrows; r++) {
        for (I c = 0; c < R_r.ncols; c++) {
            if (diagonalScaling) {
                di = A.vals[A.rows[permVec[c]] - Indexing::CSR];
            }
            R.vals[r * R.ncols + permVec[c]] = R_r.vals[r * R_r.ncols + c] / sqrt(di);
        }
        for (esint c = 0; c < R_s.ncols; c++) {
            if (diagonalScaling) {
                di = A.vals[A.rows[permVec[c + R_r.ncols]] - Indexing::CSR];
            }
            R.vals[r * R.ncols + permVec[c + R_r.ncols]] = -R_s.vals[r * R_s.ncols + c] / sqrt(di);
        }
    }

    math::orthonormalize(R);
    std::vector<I> pivots;
    _getNullPivots(R, pivots);

    if (diagonalRegularization) {
        regMat.resize(A.nrows, A.ncols, pivots.size());
        regMat.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
        regMat.shape = Matrix_Shape::UPPER;

        I row = 0;
        regMat.rows[row] = Indexing::CSR;
        for (size_t i = 0; i < pivots.size(); ++i, ++row) {
            while (row < pivots[i]) {
                regMat.rows[row + 1] = regMat.rows[row];
                ++row;
            }
            regMat.rows[row + 1] = regMat.rows[row] + 1;
            regMat.cols[regMat.rows[row] - Indexing::CSR] = pivots[i] + Indexing::CSR;
            regMat.vals[regMat.rows[row] - Indexing::CSR] = rho;
        }
        while (row < regMat.nrows) {
            regMat.rows[row + 1] = regMat.rows[row];
            ++row;
        }
    } else {
//        SparseMatrix N;
//        if (use_null_pivots_or_s_set){
//            N.CreateMatFromRowsFromMatrix( Kplus_R, null_pivots);
//        }
//        else
//        {
//            N.CreateMatFromRowsFromMatrix( Kplus_R, fix_dofs);
//        }
//        //null_pivots
//        SparseMatrix Nt;
//        N.MatTranspose( Nt );
//        SparseMatrix NtN_Mat;
//        NtN_Mat.MatMat( Nt,'N',N );
//        NtN_Mat.MatTranspose();
//        NtN_Mat.RemoveLower();
//        SparseSolverCPU NtN;
//        NtN.ImportMatrix(NtN_Mat);
//        NtN_Mat.Clear();
//        std::stringstream sss;
//        sss << "get kernel from K -> rank: " << info::mpi::rank;
//        NtN.Factorization(sss.str());
//        NtN.SolveMat_Sparse(Nt);
//        NtN.Clear();
//        NtN_Mat.MatMat(N,'N',Nt);
//        NtN_Mat.MatScale(rho);
//        NtN_Mat.RemoveLower();
//        K.MatAddInPlace (NtN_Mat,'N', 1);
//        // IF d_sub == -1, it is GGt0 of cluster and regMat is no need
//        if (d_sub!=-1)
//        {
//            regMat=NtN_Mat;
//            regMat.ConvertToCOO(1);
//        }
    }
}

template <> void orthonormalize(Matrix_Dense<double, int> &A) { _orthonormalize(A); }
template <> void permute(Matrix_CSR<double, int> &A, const std::vector<int> &perm) { _permute(A, perm); }
template <> void getKernel<double, int>(Matrix_CSR<double, int> &A, Matrix_Dense<double, int> &R, Matrix_CSR<double, int> &regMat, int maxDefect, int scSize) { _getKernel<double, int>(A, R, regMat, maxDefect, scSize); };

}
}



