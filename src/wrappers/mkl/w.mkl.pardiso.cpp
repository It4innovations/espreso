
#include "math/math.h"
#include "gpu/gpu_management.h"
#include "esinfo/eslog.h"

#include <numeric>

#ifdef HAVE_MKL
#ifdef ESPRESO_USE_WRAPPER_SPSOLVER_MKL

#include "wrappers/pardiso/w.pardiso.type.h"
#include "wrappers/pardiso/w.pardiso.h"
#include "mkl_pardiso.h"

namespace espreso {

template<typename T, typename I>
struct Solver_External_Representation
{
    PARDISOParameters<I> pp;
    const Matrix_CSR<T, I> *matrix = nullptr;
    I rows, nnzA, nnzL;
};

template<typename T>
void _pardiso(std::unique_ptr<Solver_External_Representation<T, int> > & ext, int nrhs, T *rhs, T *solution)
{
    pardiso(
            ext->pp.pt, &ext->pp.maxfct, &ext->pp.mnum,
            &ext->pp.mtype,
            &ext->pp.phase,
            &ext->matrix->nrows, ext->matrix->vals, ext->matrix->rows, ext->matrix->cols,
            ext->pp.perm, &nrhs, ext->pp.iparm, &ext->pp.msglvl,
            rhs, solution,
            &ext->pp.error);
}

template<typename T>
void _pardiso(std::unique_ptr<Solver_External_Representation<T, long> > & ext, long nrhs, T *rhs, T *solution)
{
    eslog::error("Implement MKL PARDISO 64bit interface.\n");
//    pardiso_64(
//            ext->pp.pt, &ext->pp.maxfct, &ext->pp.mnum,
//            &ext->pp.mtype,
//            &ext->pp.phase,
//            &ext->matrix->nrows, ext->matrix->vals, ext->matrix->rows, ext->matrix->cols,
//            ext->pp.perm, &nrhs, ext->pp.iparm, &ext->pp.msglvl,
//            rhs, solution,
//            &ext->pp.error);
}

template<typename T, typename I>
bool _callPardiso(I phase, std::unique_ptr<Solver_External_Representation<T,I> > & ext, I nrhs, T *rhs, T *solution)
{
    ext->pp.phase = phase;
    _pardiso(ext, nrhs, rhs, solution);

    switch (ext->pp.error) {
    case   0: break;
    case  -1: eslog::error("MKL PARDISO: input inconsistent.\n"); break;
    case  -2: eslog::error("MKL PARDISO: not enough memory.\n"); break;
    case  -3: eslog::error("MKL PARDISO: reordering problem.\n"); break;
    case  -4: eslog::error("MKL PARDISO: zero pivot, numerical factorization or iterative refinement problem.\n"); break;
    case  -5: eslog::error("MKL PARDISO: unclassified (internal) error.\n"); break;
    case  -6: eslog::error("MKL PARDISO: reordering failed.\n"); break;
    case  -7: eslog::error("MKL PARDISO: diagonal matrix is singular.\n"); break;
    case  -8: eslog::error("MKL PARDISO: 32-bit integer overflow problem.\n"); break;
    case  -9: eslog::error("MKL PARDISO: not enough memory for OOC.\n"); break;
    case -10: eslog::error("MKL PARDISO: error opening OOC files.\n"); break;
    case -11: eslog::error("MKL PARDISO: read/write error with OOC files.\n"); break;
    case -12: eslog::error("MKL PARDISO: (pardiso_64 only) pardiso_64 called from 32-bit library.\n"); break;
    case -13: eslog::error("MKL PARDISO: interrupted by the (user-defined) mkl_progress function.\n"); break;
    case -15: eslog::error("MKL PARDISO: internal error which can appear for iparm[23]=10 and iparm[12]=1. Try switch matching off (set iparm[12]=0 and rerun.\n"); break;
    }
    return ext->pp.error == 0;
}

template <typename T, typename I>
const char* DirectSparseSolver<T, I>::name()
{
    return "MKL PARDISO";
}

template <typename T, typename I>
bool DirectSparseSolver<T, I>::provideFactors()
{
    return false;
}

template <typename T, typename I>
bool DirectSparseSolver<T, I>::provideSC()
{
    return true;
}

template <typename T, typename I>
Solver_Factors DirectSparseSolver<T, I>::factorsSymmetry()
{
    return Solver_Factors::NONE;
}

template <typename T, typename I>
DirectSparseSolver<T, I>::DirectSparseSolver()
{

}

template <typename T, typename I>
DirectSparseSolver<T, I>::DirectSparseSolver(DirectSparseSolver<T, I> &&other) = default;

template <typename T, typename I>
DirectSparseSolver<T, I> & DirectSparseSolver<T, I>::operator=(DirectSparseSolver<T, I> &&other) = default;

template <typename T, typename I>
DirectSparseSolver<T, I>::~DirectSparseSolver()
{
    if (ext) {
        _callPardiso<T, I>(-1, ext, 0, nullptr, nullptr);
    }
}

template <typename T, typename I>
DirectSparseSolver<T, I>::DirectSparseSolver(const Matrix_CSR<T, I> &a)
{
    commit(a);
}

template<typename T>
void _init(std::unique_ptr<Solver_External_Representation<T, int> > & ext)
{
    pardisoinit(ext->pp.pt, &ext->pp.mtype, ext->pp.iparm);
}

template<typename T>
void _init(std::unique_ptr<Solver_External_Representation<T, long> > & ext)
{
    eslog::error("Implement MKL PARDISO 64bit interface.\n");
    // pardisoinit(ext->pp.pt, &ext->pp.mtype, ext->pp.iparm);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::commit(const Matrix_CSR<T, I> &a)
{
    ext = std::make_unique<Solver_External_Representation<T,I>>();
    ext->matrix = &a;
    ext->pp.mtype = _pardisoType(*ext->matrix);
    _init(ext);
    ext->pp.iparm[0] = 1;            /* No solver default */
    ext->pp.iparm[1] = 2;            /* Fill-in reordering from METIS */
    ext->pp.iparm[7] = 0;            /* Max number of refinement iterations (default changed between 2024.2 and 2025.0) */
    ext->pp.iparm[9] = 13;           /* Perturb the pivot elements with 1E-13 */

    ext->pp.iparm[4] = 2;            /* Return permutation vector */
    ext->pp.iparm[34] = 1 - Indexing::CSR; /* One- or Zero-based indexing */
    ext->pp.perm = new I[ext->matrix->nrows];
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::symbolicFactorization()
{
    _callPardiso<T, I>(11, ext, 0, nullptr, nullptr);
    ext->rows = ext->matrix->nrows;
    ext->nnzA = ext->matrix->nnz;
    // ext->memoryL = 1024 * (ext->pp.iparm[15] + ext->pp.iparm[16]);
    ext->nnzL = ext->pp.iparm[17];
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::numericalFactorization()
{
    _callPardiso<T, I>(22, ext, 0, nullptr, nullptr);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solve(const Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution)
{
    solution.resize(rhs.size);
    _callPardiso<T, I>(33, ext, 1, rhs.vals, solution.vals);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solve(const Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution)
{
    solution.resize(rhs.nrows, rhs.ncols);
    _callPardiso<T, I>(33, ext, rhs.nrows, rhs.vals, solution.vals);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveForward (const Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution)
{
    solution.resize(rhs.size);
    _callPardiso<T, I>(331, ext, 1, rhs.vals, solution.vals);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveDiagonal(const Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution)
{
    solution.resize(rhs.size);
    _callPardiso<T, I>(332, ext, 1, rhs.vals, solution.vals);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveBackward(const Vector_Dense<T, I> &rhs, Vector_Dense<T, I> &solution)
{
    solution.resize(rhs.size);
    _callPardiso<T, I>(333, ext, 1, rhs.vals, solution.vals);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveForward (const Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution)
{
    solution.resize(rhs.nrows, rhs.ncols);
    _callPardiso<T, I>(331, ext, rhs.nrows, rhs.vals, solution.vals);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveDiagonal(const Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution)
{
    solution.resize(rhs.nrows, rhs.ncols);
    _callPardiso<T, I>(332, ext, rhs.nrows, rhs.vals, solution.vals);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::solveBackward(const Matrix_Dense<T, I> &rhs, Matrix_Dense<T, I> &solution)
{
    solution.resize(rhs.nrows, rhs.ncols);
    _callPardiso<T, I>(333, ext, rhs.nrows, rhs.vals, solution.vals);
}

template <typename T, typename I>
I DirectSparseSolver<T, I>::getMatrixSize()
{
    return ext->rows;
}

template <typename T, typename I>
I DirectSparseSolver<T, I>::getMatrixNnz()
{
    return ext->nnzA;
}

template <typename T, typename I>
I DirectSparseSolver<T, I>::getFactorNnz()
{
    return ext->nnzL;
}

template <typename T, typename I>
template<typename A>
void DirectSparseSolver<T, I>::getFactorL(Matrix_CSR<T, I, A> &/*L*/, bool /*copyPattern*/, bool /*copyValues*/)
{
    eslog::error("MKL PARDISO does not provide factors.\n");
}

template <typename T, typename I>
template<typename A>
void DirectSparseSolver<T, I>::getFactorU(Matrix_CSR<T, I, A> &/*U*/, bool /*copyPattern*/, bool /*copyValues*/)
{
    eslog::error("MKL PARDISO does not provide factors.\n");
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::getPermutation(Permutation<I> &perm)
{
    if(ext->pp.phase < 11) eslog::error("getPermutation: invalid order of operations in spsolver\n");

    if(perm.size != ext->matrix->nrows) {
        perm.resize(ext->matrix->nrows);
    }

    std::copy_n(ext->pp.perm, ext->matrix->nrows, perm.dst_to_src);
    perm.invert(perm.dst_to_src, perm.src_to_dst);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::getPermutation(Vector_Dense<I> &perm)
{
    if(ext->pp.phase < 11) eslog::error("getPermutation: invalid order of operations in spsolver\n");

    if(perm.size != ext->matrix->nrows) {
        perm.resize(ext->matrix->nrows);
    }

    std::copy_n(ext->pp.perm, ext->matrix->nrows, perm.vals);
}

template <typename T, typename I>
void DirectSparseSolver<T, I>::getSC(Matrix_Dense<T, I> &sc, std::vector<I> &indices, bool symmetric_packed)
{
    Matrix_Dense<T, I> full;
    T * out_vals = sc.vals;
    if(symmetric_packed) {
        full.resize(sc.nrows, sc.nrows);
        out_vals = full.vals;
    }

    ext->pp.iparm[35] = 1;
    I *perm = ext->pp.perm;
    ext->pp.perm = indices.data();
    _callPardiso<T, I>(12, ext, 0, nullptr, out_vals);
    ext->pp.perm = perm;
    ext->pp.iparm[35] = 0;

    if(symmetric_packed) {
        for (I r = 0, i = 0; r < full.nrows; ++r) {
            for (I c = r; c < full.ncols; ++c, ++i) {
                sc.vals[i] = full.vals[r * full.ncols + c];
            }
        }
    }
}

}

#include "math/wrappers/math.spsolver.inst.hpp"

#endif
#endif
