
#ifdef HAVE_SUITESPARSE

#include "wrappers/suitesparse/operations/solver_csx.cholmod.h"

#include "wrappers/suitesparse/w.suitesparse.cholmod.h"
#include "math/primitives_new/allocator_new.h"
#include "math/operations/copy_dnx.h"
#include "math/operations/convert_csx_dny.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
struct solver_csx_cholmod_data
{
    cholmod_common cm_common;
    cholmod_factor * cm_factor_super = nullptr;
    cholmod_factor * cm_factor_simpl = nullptr;
    cholmod_sparse cm_A_view;
    VectorDenseData_new<I> map_simpl_super;
};



template<typename T, typename I>
solver_csx_cholmod<T,I>::solver_csx_cholmod()
{
    data = std::make_unique<solver_csx_cholmod_data<T,I>>();
}



template<typename T, typename I>
solver_csx_cholmod<T,I>::~solver_csx_cholmod()
{
    if (data->cm_factor_super != nullptr) _free<I>(data->cm_factor_super, data->cm_common);

    if (data->cm_factor_simpl != nullptr) _free<I>(data->cm_factor_simpl, data->cm_common);

    _finish<I>(data->cm_common);
}



template<typename T, typename I>
void solver_csx_cholmod<T,I>::internal_factorize_symbolic()
{
    if(!is_hermitian<T>(A->prop.symm)) eslog::error("matrix has to be hermitian\n");
    if(A->prop.dfnt != MatrixDefinitness_new::positive_definite) eslog::error("matrix has to be positive definite\n");

    _start<I>(data->cm_common);
    data->cm_common.final_ll = 1;
    data->cm_common.nthreads_max = 1;
    data->cm_common.nmethods = 1;
    data->cm_common.method[0].ordering = CHOLMOD_METIS;
    data->cm_common.itype = _getCholmodItype<I>();
    data->cm_common.supernodal = CHOLMOD_SUPERNODAL;

    MatrixCsxView_new<T,I> A_rt = A->get_transposed_reordered_view();
    MatrixCsxView_new<T,I> * A_to_use = (A->order == 'C') ? A : &A_rt;

    // cholmod assumes CSC
    data->cm_A_view.nrow = A_to_use->nrows;
    data->cm_A_view.ncol = A_to_use->ncols;
    data->cm_A_view.nzmax = A_to_use->nnz;
    data->cm_A_view.nz = nullptr;
    data->cm_A_view.z = nullptr;
    data->cm_A_view.stype = _getCholmodStype(A_to_use->prop.uplo);
    data->cm_A_view.itype = _getCholmodItype<I>();
    data->cm_A_view.xtype = _getCholmodXtype<T>();
    data->cm_A_view.dtype = _getCholmodDtype<T>();
    data->cm_A_view.sorted = 1;
    data->cm_A_view.packed = 1;

    data->cm_A_view.p = A->ptrs;
    data->cm_A_view.i = A->idxs;
    data->cm_A_view.x = A->vals;

    if (data->cm_A_view.nrow == 0 || data->cm_A_view.ncol == 0) return;

    data->cm_factor_super = _analyze<I>(&data->cm_A_view, data->cm_common);

    if(data->cm_common.lnz > (double)std::numeric_limits<I>::max()) {
        eslog::error("symbolicFactorization: factor nnz too large for the used integer type\n");
    }
    // https://github.com/DrTimothyAldenDavis/SuiteSparse/issues/523
    size_t factor_nnz = static_cast<I>(data->cm_common.lnz);
    nnz_L = factor_nnz;
    nnz_U = factor_nnz;

    if(need_factors) {
        // preprocess extraction of factors from supernodal factor directly to csx matrix

        if (data->cm_factor_super->xsize > utils::get_max_val_no_precision_loss_in_fp<T>()) {
            eslog::error("symbolicFactorization: factor nnz too large for my super->simpl map\n");
        }

        data->cm_factor_simpl = _copyFactor<I>(data->cm_factor_super, data->cm_common);
        _changeFactor<I>(_getCholmodXtype<T>(), true, true, true, true, data->cm_factor_simpl, data->cm_common);

        for(size_t i = 0; i < data->cm_factor_simpl->xsize; i++) {
            reinterpret_cast<T*>(data->cm_factor_simpl->x)[i] = static_cast<T>(i);
        }

        _changeFactor<I>(_getCholmodXtype<T>(), true, false, true, true, data->cm_factor_simpl, data->cm_common);

        _resymbol<I>(&data->cm_A_view, nullptr, 0, 1, data->cm_factor_simpl, data->cm_common);

        if (data->cm_factor_simpl->nzmax != factor_nnz) {
            eslog::error("symbolicFactorization: some weird error in cholmod. %zu %zu\n", data->cm_factor_simpl->nzmax, factor_nnz);
        }

        data->map_simpl_super.set(factor_nnz, AllocatorCPU_new::get_singleton());
        data->map_simpl_super.alloc();
        for(size_t i = 0; i < data->map_simpl_super.size; i++) {
            data->map_simpl_super.vals[i] = static_cast<I>(std::real(reinterpret_cast<T*>(data->cm_factor_simpl->x)[i]));
        }
    }
}



template<typename T, typename I>
void solver_csx_cholmod<T,I>::internal_factorize_numeric()
{
    data->cm_A_view.p = A->ptrs;
    data->cm_A_view.i = A->idxs;
    data->cm_A_view.x = A->vals;

    if (data->cm_A_view.nrow == 0 || data->cm_A_view.ncol == 0) return;

    _factorize<I>(data->cm_factor_super, &data->cm_A_view, data->cm_common);
}



template<typename T, typename I>
void solver_csx_cholmod<T,I>::internal_get_permutation(PermutationView_new<I> & perm)
{
    std::copy_n(static_cast<I*>(data->cm_factor_super->Perm), data->cm_factor_super->n, perm.dst_to_src);

    if (data->cm_factor_super->IPerm != nullptr) {
        std::copy_n(static_cast<I*>(data->cm_factor_super->IPerm), data->cm_factor_super->n, perm.src_to_dst);
    }
    else {
        PermutationView_new<I>::invert(perm.dst_to_src, perm.src_to_dst, perm.size);
    }
}



template<typename T, typename I>
void solver_csx_cholmod<T,I>::get_factor_impl(MatrixCsxView_new<T,I> & factor, bool pattern, bool values)
{
    // correct order/uplo/size/nnz is handled in base class, here I will just copy the data

    if (pattern) {
        std::copy_n(static_cast<I*>(data->cm_factor_simpl->p), data->cm_factor_simpl->n+1, factor.ptrs);
        std::copy_n(static_cast<I*>(data->cm_factor_simpl->i), data->cm_factor_simpl->nzmax, factor.idxs);
    }
    if (values) {
        for(size_t i = 0; i < data->map_simpl_super.size; i++) {
            factor.vals[i] = reinterpret_cast<T*>(data->cm_factor_super->x)[data->map_simpl_super.vals[i]];
        }
    }
}



template<typename T, typename I>
void solver_csx_cholmod<T,I>::internal_get_factor_L(MatrixCsxView_new<T,I> & L, bool pattern, bool values)
{
    get_factor_impl(L, pattern, values);
}



template<typename T, typename I>
void solver_csx_cholmod<T,I>::internal_get_factor_U(MatrixCsxView_new<T,I> & U, bool pattern, bool values)
{
    get_factor_impl(U, pattern, values);
}



template<typename T, typename I>
void solver_csx_cholmod<T,I>::internal_solve(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
{
    MatrixDenseView_new<T> rhs_mat;
    rhs_mat.set_view(rhs.size, 1, rhs.size, 'C', rhs.vals, rhs.ator);

    MatrixDenseView_new<T> sol_mat;
    sol_mat.set_view(sol.size, 1, sol.size, 'C', sol.vals, sol.ator);

    this->internal_solve(rhs_mat, sol_mat);
}



template<typename T, typename I>
void solver_csx_cholmod<T,I>::internal_solve(MatrixDenseView_new<T> & rhs, MatrixDenseView_new<T> & sol)
{
    if(rhs.order != 'C') eslog::error("only support colmajor rhs/sol for now\n");
    if (data->cm_A_view.nrow == 0 || data->cm_A_view.ncol == 0) return;

    cholmod_dense cm_rhs;
    cm_rhs.nrow = rhs.nrows;
    cm_rhs.ncol = rhs.ncols;
    cm_rhs.d = rhs.ld;
    cm_rhs.nzmax = cm_rhs.d * rhs.nrows;
    cm_rhs.x = rhs.vals;
    cm_rhs.xtype = _getCholmodXtype<T>();
    cm_rhs.dtype = _getCholmodDtype<T>();

    cholmod_dense * cm_sol = _solve<I>(CHOLMOD_A, data->cm_factor_super, &cm_rhs, data->cm_common);

    {
        MatrixDenseView_new<T> sol_from_cm;
        sol_from_cm.set_view(cm_sol->nrow, cm_sol->ncol, cm_sol->d, 'C', reinterpret_cast<T*>(cm_sol->x), AllocatorDummy_new::get_singleton(true, false));

        copy_dnx<T>::do_all(&sol_from_cm, &sol);
    }

    _free<I>(cm_sol, data->cm_common);
}



template<typename T, typename I>
void solver_csx_cholmod<T,I>::internal_solve(MatrixCsxView_new<T,I> & rhs, MatrixDenseView_new<T> & sol)
{
    if(sol.order != 'C') eslog::error("only support colmajor sol for now\n");

    MatrixDenseData_new<T> rhs_dn;
    rhs_dn.set(rhs.nrows, rhs.ncols, 'C', AllocatorCPU_new::get_singleton());
    rhs_dn.alloc();

    convert_csx_dny<T,I>::do_all(&rhs, &rhs_dn);

    this->internal_solve(rhs_dn, sol);
}



#define INSTANTIATE_T_I(T,I) \
template class solver_csx_cholmod<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T, int32_t) \
    /* INSTANTIATE_T_I(T, int64_t) */

        #define INSTANTIATE \
        /* INSTANTIATE_T(float) */ \
        INSTANTIATE_T(double) \
        /* INSTANTIATE_T(std::complex<float>) */ \
        INSTANTIATE_T(std::complex<double>)

            INSTANTIATE

        #undef INSTANTIATE
    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
}
}

#endif
