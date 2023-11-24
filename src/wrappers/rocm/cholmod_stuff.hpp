

#pragma once

extern "C" {
#include <cholmod.h>
}

#include "matrices.hpp"

namespace espreso {

template<typename I>
void my_cholmod_start(cholmod_common * common)
{
    if constexpr(sizeof(I) == 4) cholmod_start(common);
    if constexpr(sizeof(I) == 8) cholmod_l_start(common);
}
template<typename I>
cholmod_factor * my_cholmod_analyze(cholmod_sparse * A, cholmod_common * common)
{
    if constexpr(sizeof(I) == 4) return cholmod_analyze(A, common);
    if constexpr(sizeof(I) == 8) return cholmod_l_analyze(A, common);
}
template<typename I>
void my_cholmod_factorize(cholmod_sparse * A, cholmod_factor * F, cholmod_common * common)
{
    if constexpr(sizeof(I) == 4) cholmod_factorize(A, F, common);
    if constexpr(sizeof(I) == 8) cholmod_l_factorize(A, F, common);
}
template<typename I>
cholmod_sparse * my_cholmod_factor_to_sparse(cholmod_factor * F, cholmod_common * common)
{
    if constexpr(sizeof(I) == 4) return cholmod_factor_to_sparse(F, common);
    if constexpr(sizeof(I) == 8) return cholmod_l_factor_to_sparse(F, common);
}
template<typename I>
cholmod_dense * my_cholmod_solve(int sys, cholmod_factor * L, cholmod_dense * B, cholmod_common * common)
{
    if constexpr(sizeof(I) == 4) return cholmod_solve(sys, L, B, common);
    if constexpr(sizeof(I) == 8) return cholmod_l_solve(sys, L, B, common);
}
template<typename I>
void my_cholmod_free_sparse(cholmod_sparse ** A, cholmod_common * common)
{
    if constexpr(sizeof(I) == 4) cholmod_free_sparse(A, common);
    if constexpr(sizeof(I) == 8) cholmod_l_free_sparse(A, common);
}
template<typename I>
void my_cholmod_free_factor(cholmod_factor ** F, cholmod_common * common)
{
    if constexpr(sizeof(I) == 4) cholmod_free_factor(F, common);
    if constexpr(sizeof(I) == 8) cholmod_l_free_factor(F, common);
}
template<typename I>
int my_cholmod_sdmult(cholmod_sparse * A, int transpose, double alpha[2], double beta[2], cholmod_dense * X, cholmod_dense * Y, cholmod_common * Common)
{
    if constexpr(sizeof(I) == 4) return cholmod_sdmult(A, transpose, alpha, beta, X, Y, Common);
    if constexpr(sizeof(I) == 8) return cholmod_l_sdmult(A, transpose, alpha, beta, X, Y, Common);
}
template<typename I>
cholmod_dense * my_cholmod_sparse_to_dense(cholmod_sparse * A, cholmod_common * Common)
{
    if constexpr(sizeof(I) == 4) return cholmod_sparse_to_dense(A, Common);
    if constexpr(sizeof(I) == 8) return cholmod_l_sparse_to_dense(A, Common);
}
template<typename I>
cholmod_dense * my_cholmod_allocate_dense(size_t nrow, size_t ncol, size_t d, int xtype, cholmod_common *Common)
{
    if constexpr(sizeof(I) == 4) return cholmod_allocate_dense(nrow, ncol, d, xtype, Common);
    if constexpr(sizeof(I) == 8) return cholmod_l_allocate_dense(nrow, ncol, d, xtype, Common);
}
template<typename I>
int my_cholmod_free_dense(cholmod_dense **X, cholmod_common *Common)
{
    if constexpr(sizeof(I) == 4) return cholmod_free_dense(X, Common);
    if constexpr(sizeof(I) == 8) return cholmod_l_free_dense(X, Common);
}
template<typename I>
cholmod_factor * my_cholmod_copy_factor(cholmod_factor * F, cholmod_common * Common)
{
    if constexpr(sizeof(I) == 4) return cholmod_copy_factor(F, Common);
    if constexpr(sizeof(I) == 8) return cholmod_l_copy_factor(F, Common);
}

template<typename T, typename I>
void my_cholmod_start_and_init(cholmod_common ** cm_common, char ordering = 'M')
{
    *cm_common = new cholmod_common();
    cholmod_common * common = *cm_common;
    my_cholmod_start<I>(common);
    common->final_ll = 1;
    common->nthreads_max = 1;
    common->nmethods = 1;
    common->method[0].ordering = (ordering == 'A' ? CHOLMOD_AMD : (ordering == 'M' ? CHOLMOD_METIS : CHOLMOD_NATURAL));
//    common->dtype = (std::is_same_v<T,double>) ? (CHOLMOD_DOUBLE) : (CHOLMOD_SINGLE);
    common->itype = (std::is_same_v<I,int64_t>) ? (CHOLMOD_LONG) : (CHOLMOD_INT);
    common->supernodal = CHOLMOD_SUPERNODAL;
}

void my_cholmod_finish_and_delete(cholmod_common ** cm_common)
{
    cholmod_finish(*cm_common);
    delete *cm_common;
    *cm_common = nullptr;
}

template<typename T, typename I, template<typename> typename A>
cholmod_sparse my_cholmod_sparse_view(const MatrixCSR<T,I,A> & matrix, int triangle)
{
    // triangle: -1 lower, 1 upper, 0 general matrix
    cholmod_sparse cm_sparse;
    cm_sparse.nrow = matrix.ncols;
    cm_sparse.ncol = matrix.nrows;
    cm_sparse.nzmax = matrix.nvals;
    cm_sparse.p = matrix.rowptrs;
    cm_sparse.i = matrix.colidxs;
    cm_sparse.nz = nullptr;
    cm_sparse.x = matrix.vals;
    cm_sparse.z = nullptr;
    cm_sparse.stype = -1 * triangle;
    cm_sparse.itype = (sizeof(I) == 4) ? (CHOLMOD_INT) : (CHOLMOD_LONG);
    cm_sparse.xtype = CHOLMOD_REAL;
    cm_sparse.dtype = (std::is_same_v<T,double>) ? (CHOLMOD_DOUBLE) : (CHOLMOD_SINGLE);
    cm_sparse.sorted = 1;
    cm_sparse.packed = 1;
    return cm_sparse;
}

template<typename I>
size_t my_cholmod_factor_get_nnz(cholmod_factor * cm_factor)
{
    size_t lnz = 0;

    if(cm_factor->is_super)
    {
        I * Lpi = (I*)cm_factor->pi ;
        I * Super = (I*)cm_factor->super ;

        for (I s = 0 ; s < cm_factor->nsuper ; s++)
        {
            I k1 = Super [s] ;
            I k2 = Super [s+1] ;
            I psi = Lpi [s] ;
            I psend = Lpi [s+1] ;
            I nsrow = psend - psi ;
            I nscol = k2 - k1 ;

            lnz += nscol * nsrow - (nscol*nscol - nscol)/2 ;
        }
    }
    else
    {
        MY_ABORT("Only supernodal factors are supported");
    }

    return lnz;
}

template<typename T, typename I, template<typename> typename A>
void my_cholmod_extract_factor_matrix(cholmod_factor * cm_factor, MatrixCSR<T,I,A> & matrix, cholmod_common * cm_common, bool copy_pattern = true, bool copy_vals = true)
{
    if(matrix.nrows != cm_factor->n || matrix.ncols != cm_factor->n) MY_ABORT("Output matrix has wrong dimensions");

    cholmod_factor * factor_copy = my_cholmod_copy_factor<I>(cm_factor, cm_common);
    cholmod_sparse * cm_L = my_cholmod_factor_to_sparse<I>(factor_copy, cm_common);
    if(matrix.nvals != cm_L->nzmax) MY_ABORT("Output matrix has wrong dimensions");

    if(copy_pattern) std::copy_n(static_cast<I*>(cm_L->p), cm_L->nrow+1, matrix.rowptrs);
    if(copy_pattern) std::copy_n(static_cast<I*>(cm_L->i), cm_L->nzmax,  matrix.colidxs);
    if(copy_vals)    std::copy_n(static_cast<T*>(cm_L->x), cm_L->nzmax,  matrix.vals);

    my_cholmod_free_sparse<I>(&cm_L, cm_common);
    my_cholmod_free_factor<I>(&factor_copy, cm_common);
}

template<typename I, template<typename> typename A>
void my_cholmod_extract_factor_perm(cholmod_factor * cm_factor, Permutation<I,A> & perm, cholmod_common * cm_common)
{
    if(perm.size != cm_factor->n) MY_ABORT("Output permutation has wrong size");

    std::copy_n(static_cast<I*>(cm_factor->Perm), cm_factor->n, perm.forward);
    inverse_permutation(perm, PermutationDirection::Backward);
}



struct tm_cholmodfact
{
    my_timer symbolic, numeric, extract;
};

template<typename T, typename I, template<typename> typename Au, template<typename> typename Ap, template<typename> typename Ak>
void factorize_with_cholmod(MatrixCSR<T,I,Au> * out_Ufactor, Permutation<I,Ap> & out_perm, const MatrixCSR<T,I,Ak> & Kreg, char ordering, tm_cholmodfact & tm, cholmod_common * cm_common = nullptr)
{
    if(out_perm.size != Kreg.nrows)
    {
        MY_ABORT("factorize_with_cholmod: output permutation has wrong dimension");
    }
    bool delete_common = false;
    if(cm_common == nullptr)
    {
        my_cholmod_start_and_init<T,I>(&cm_common, ordering);
        delete_common = true;
    }

    cholmod_sparse cm_Kreg = my_cholmod_sparse_view(Kreg, 1);

    tm.symbolic.start();
    cholmod_factor * cm_fac = my_cholmod_analyze<I>(&cm_Kreg, cm_common);
    tm.symbolic.stop();
    tm.numeric.start();
    my_cholmod_factorize<I>(&cm_Kreg, cm_fac, cm_common);
    tm.numeric.stop();

    tm.extract.start();
    out_Ufactor->resize(cm_fac->n, cm_fac->n, my_cholmod_factor_get_nnz<I>(cm_fac), true);
    my_cholmod_extract_factor_matrix(cm_fac, *out_Ufactor, cm_common, true, true);
    my_cholmod_extract_factor_perm(cm_fac, out_perm, cm_common);
    tm.extract.stop();

    my_cholmod_free_factor<I>(&cm_fac, cm_common);

    if(delete_common)
    {
        my_cholmod_finish_and_delete(&cm_common);
    }
}


template<typename T, typename I>
void print_matrix(const cholmod_dense & M, const char * name)
{
    printf("Matrix %s\n", name);
    printf("Size %lldx%lld, ld %lld\n", (long long)M.nrow, (long long)M.ncol, (long long)M.d);
    printf("Printed transposed because of colmajor\n");
    for(I c = 0; c < M.ncol; c++)
    {
        for(I r = 0; r < M.nrow; r++)
        {
            if constexpr(std::is_floating_point_v<T>)
            {
                double v = reinterpret_cast<double*>(M.x)[c * M.d + r];
                char str[100];
                sprintf(str, "%+10.3e", v);
                if(strstr(str, "nan") != nullptr) printf("  nan      ");
                else if(strstr(str, "inf") != nullptr) printf(" %cinf      ", v > 0 ? '+' : '-');
                else if(std::abs(v) >= 1e100) printf("  INF      ");
                else if(std::abs(v) <= 1e-100) printf("  0        ");
                else printf(" %+10.3e", v);
            }
            if constexpr(std::is_integral_v<T>) printf(" %+9lld", (long long)reinterpret_cast<I*>(M.x)[c * M.d + r]);
        }
        printf("\n");
    }
    fflush(stdout);
}

}
