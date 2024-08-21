
#ifdef HAVE_ONEAPI

#include "gpu/gpu_spblas.h"
#include "w.oneapi.gpu_management.h"
#include "basis/utilities/utils.h"
#include <oneapi/mkl.hpp>

namespace espreso {
namespace gpu {
namespace spblas {

    namespace onemkl = oneapi::mkl;
    namespace onesparse = onemkl::sparse;

    namespace
    {
        static inline onemkl::uplo _char_to_uplofill(char c)
        {
            switch(c) {
                case 'U': return onemkl::uplo::upper;
                case 'L': return onemkl::uplo::lower;
                default: eslog::error("invalid fill '%c'\n", c);
            }
        }

        static inline onemkl::transpose _char_to_operation(char c)
        {
            switch(c) {
                case 'N': return onemkl::transpose::nontrans;
                case 'T': return onemkl::transpose::trans;
                case 'H': return onemkl::transpose::conjtrans;
                default: eslog::error("invalid operation '%c'\n", c);
            }
        }

        static inline onemkl::layout _char_to_layout(char c)
        {
            switch(c) {
                case 'R': return onemkl::layout::row_major;
                case 'C': return onemkl::layout::col_major;
                default: eslog::error("invalid layout '%c'\n", c);
            }
        }
        
        template<typename T, typename I, char order>
        struct kernel_functor_sp2dn
        {
            I * sparse_rowptrs;
            I * sparse_colidxs;
            T * sparse_vals;
            T * dense_vals;
            I dense_ld;
            kernel_functor_sp2dn(I*sr, I*sc, T*sv, T*dv, I ld)
                : sparse_rowptrs(sr)
                , sparse_colidxs(sc)
                , sparse_vals(sv)
                , dense_vals(dv)
                , dense_ld(ld)
            {}
            void operator()(sycl::nd_item<1> item) const {
                sycl::group g = item.get_group();
                I row = g.get_group_linear_id();
                I start = sparse_rowptrs[row];
                I end = sparse_rowptrs[row+1];
                for(I i = start + g.get_local_linear_id(); i < end; i += g.get_local_linear_range()) {
                    I col = sparse_colidxs[i];
                    I val = sparse_vals[i];
                    if constexpr(order == 'R') dense_vals[row * dense_ld + col] = val;
                    if constexpr(order == 'C') dense_vals[row + dense_ld * col] = val;
                }
            }
        };
    }

    spblas_wrapper_impl get_implementation()
    {
        return spblas_wrapper_impl::ONEMKL_SPARSE;
    }
    
    struct _handle
    {
        mgm::queue q;
        sycl::queue qq;
    };

    struct _descr_matrix_csr
    {
        onesparse::matrix_handle_t d;
        void * rowptrs = nullptr;
        void * colidxs = nullptr;
        void * vals = nullptr;
        int64_t nrows = -1;
        int64_t nnz = -1;
        char fill = '_';
    };

    struct _descr_matrix_dense
    {
        void * vals = nullptr;
        int64_t nrows = -1;
        int64_t ncols = -1;
        int64_t ld = -1;
        char order = '_';
        _descr_matrix_dense get_complementary()
        {
            _descr_matrix_dense ret = *this;
            std::swap(ret.nrows, ret.ncols);
            ret.order = mgm::order_change(ret.order);
            return ret;
        }
    };

    struct _descr_vector_dense
    {
        void * vals = nullptr;
    };

    struct _descr_sparse_trsv {};

    struct _descr_sparse_trsm {};

    struct _descr_sparse_mv {};

    void handle_create(handle & h, mgm::queue & q)
    {
        h = std::make_shared<_handle>();
        h->q = q;
        h->qq = q->q;
    }

    void handle_destroy(handle & h)
    {
        h.reset();
    }

    template<typename T, typename I>
    void descr_matrix_csr_create(handle & h, descr_matrix_csr & descr, I nrows, I ncols, I nnz, char fill)
    {
        descr = std::make_shared<_descr_matrix_csr>();
        onesparse::init_matrix_handle(&descr->d);
        descr->nrows = nrows;
        descr->nnz = nnz;
        descr->fill = fill;
    }

    template<typename T, typename I, typename A>
    void descr_matrix_csr_link_data(handle & h, descr_matrix_csr & descr, Matrix_CSR<T,I,A> & matrix)
    {
        onesparse::set_csr_data(h->qq, descr->d, matrix.nrows, matrix.ncols, onemkl::index_base::zero, matrix.rows, matrix.cols, matrix.vals);
    }

    void descr_matrix_csr_destroy(handle & h, descr_matrix_csr & descr)
    {
        onesparse::release_matrix_handle(h->qq, &descr->d);

        descr.reset();
    }

    template<typename T, typename I>
    void descr_matrix_dense_create(handle & h, descr_matrix_dense & descr, I nrows, I ncols, I ld, char order)
    {
        descr = std::make_shared<_descr_matrix_dense>();
        descr->order = order;
    }

    template<typename T, typename I, typename A>
    void descr_matrix_dense_link_data(handle & h, descr_matrix_dense & descr, Matrix_Dense<T,I,A> & matrix)
    {
        descr->vals = matrix.vals;
        descr->nrows = matrix.nrows;
        descr->ncols = matrix.ncols;
        descr->ld = matrix.get_ld();
    }

    void descr_matrix_dense_destroy(handle & h, descr_matrix_dense & descr)
    {
        descr.reset();
    }

    template<typename T, typename I>
    void descr_vector_dense_create(handle & h, descr_vector_dense & descr, I size)
    {
        descr = std::make_shared<_descr_vector_dense>();
    }

    template<typename T, typename I, typename A>
    void descr_vector_dense_link_data(handle & h, descr_vector_dense & descr, Vector_Dense<T,I,A> & vector)
    {
        descr->vals = vector.vals;
    }

    template<typename T, typename I, typename A>
    void descr_vector_dense_link_data(handle & h, descr_vector_dense & descr, Matrix_Dense<T,I,A> & matrix, I colidx)
    {
        descr->vals = matrix.vals + colidx * matrix.get_ld();
    }

    void descr_vector_dense_destroy(handle & h, descr_vector_dense & descr)
    {
        descr.reset();
    }

    void descr_sparse_trsv_create(handle & h, descr_sparse_trsv & descr) {}

    void descr_sparse_trsv_destroy(handle & h, descr_sparse_trsv & descr) {}

    void descr_sparse_trsm_create(handle & h, descr_sparse_trsm & descr) {}

    void descr_sparse_trsm_destroy(handle & h, descr_sparse_trsm & descr) {}

    void descr_sparse_mv_create(handle & h, descr_sparse_mv & descr) {}

    void descr_sparse_mv_destroy(handle & h, descr_sparse_mv & descr) {}

    template<typename T, typename I>
    void transpose(handle & h, descr_matrix_csr & output, descr_matrix_csr & input, bool conjugate, size_t & buffersize, void * /*buffer*/, char stage)
    {
        if constexpr(utils::is_real<T>()) conjugate = false;
        char op = (conjugate ? 'H' : 'T');
        if(stage == 'B') buffersize = 0;
        // if(stage == 'P') ;
        if(stage == 'C') onesparse::omatcopy(h->qq, _char_to_operation(op), input->d, output->d);
    }

    template<typename T, typename I>
    void sparse_to_dense(handle & h, char op, descr_matrix_csr & sparse, descr_matrix_dense & dense, size_t & /*buffersize*/, void * /*buffer*/, char stage)
    {
        h->qq.fill(dense->vals, T{0}, dense->nrows * dense->ld);

        I avg_nnz_per_row = sparse->nnz / sparse->nrows;
        int workitems_in_workgroup = std::clamp(avg_nnz_per_row, 64, 1024);
        sycl::nd_range<1> range(sycl::range<1>(sparse->nrows * workitems_in_workgroup), sycl::range<1>(workitems_in_workgroup));
        if(dense->order == 'R') {
            h->qq.parallel_for(
                range,
                kernel_functor_sp2dn<T,I,'R'>((I*)sparse->rowptrs, (I*)sparse->colidxs, (T*)sparse->vals, (T*)dense->vals, dense->ld)
            );
        }
        else if(dense->order == 'C') {
            h->qq.parallel_for(
                range,
                kernel_functor_sp2dn<T,I,'C'>((I*)sparse->rowptrs, (I*)sparse->colidxs, (T*)sparse->vals, (T*)dense->vals, dense->ld)
            );
        }
        else {
            eslog::error("invalid order '%c'\n", dense->order);
        }
    }

    template<typename T, typename I>
    void trsv(handle & h, char op, descr_matrix_csr & matrix, descr_vector_dense & rhs, descr_vector_dense & sol, descr_sparse_trsv & /*descr_trsv*/, size_t & buffersize, void * /*buffer*/, char stage)
    {
        T one = 1.0;
        if(stage == 'B') buffersize = 0;
        // if(stage == 'P') ;
        // if(stage == 'U') ;
        if(stage == 'C') onesparse::trsv(h->qq, _char_to_uplofill(matrix->fill), _char_to_operation(op), onemkl::diag::nonunit, one, matrix->d, (T*)rhs->vals, (T*)sol->vals);
    }

    template<typename T, typename I>
    void trsm(handle & h, char op_mat, char op_rhs, char op_sol, descr_matrix_csr & matrix, descr_matrix_dense & rhs, descr_matrix_dense & sol, descr_sparse_trsm & descr_trsm, size_t & buffersize, void * buffer, char stage)
    {
        if(op_sol == 'N') {
            if(rhs->order == sol->order) {
                T one = 1.0;
                if(stage == 'B') buffersize = 0;
                // if(stage == 'P') ;
                // if(stage == 'U') ;
                if(stage == 'C') onesparse::trsm(h->qq, _char_to_layout(sol->order), _char_to_operation(op_mat), _char_to_operation(op_rhs), _char_to_uplofill(matrix->fill), onemkl::diag::nonunit, one, matrix->d, (T*)rhs->vals, sol->ncols, rhs->ld, (T*)sol->vals, sol->ld);
            }
            else {
                descr_matrix_dense rhs_compl = std::make_shared<_descr_matrix_dense>(rhs->get_complementary());
                char op_rhs_compl = mgm::operation_combine(op_rhs, 'T');
                trsm<T,I>(h, op_mat, op_rhs_compl, op_sol, matrix, rhs_compl, sol, descr_trsm, buffersize, buffer, stage);
            }
        }
        else if(op_sol == 'T') {
            descr_matrix_dense sol_compl = std::make_shared<_descr_matrix_dense>(sol->get_complementary());
            char op_sol_compl = mgm::operation_combine(op_sol, 'T'); // 'N'
            trsm<T,I>(h, op_mat, op_rhs, op_sol_compl, matrix, rhs, sol_compl, descr_trsm, buffersize, buffer, stage);
        }
        else {
            eslog::error("op_sol '%c' not supported\n", op_sol);
            // char op_mat_compl = mgm::operation_combine(op_mat, 'C');
            // char op_rhs_compl = mgm::operation_combine(op_rhs, 'C');
            // char op_sol_compl = mgm::operation_combine(op_sol, 'C');
            // trsm<T,I>(h, op_mat_compl, op_rhs_compl, op_sol_compl, matrix, rhs, sol, descr_trsm, buffersize, buffer, stage);
        }
    }

    template<typename T, typename I>
    void mv(handle & h, char op, descr_matrix_csr & A, descr_vector_dense & x, descr_vector_dense & y, descr_sparse_mv & descr_mv, size_t & buffersize, void * /*buffer*/, char stage)
    {
        T zero = 0.0;
        T one = 1.0;
        if(stage == 'B') buffersize = 0;
        // if(stage == 'P') ;
        if(stage == 'C') onesparse::gemv(h->qq, _char_to_operation(op), one, A->d, (T*)x->vals, zero, (T*)y->vals);
    }

    template<typename T, typename I>
    void mm(handle & h, char op_A, char op_B, descr_matrix_csr & A, descr_matrix_dense & B, descr_matrix_dense & C, size_t & buffersize, void * buffer, char stage)
    {
        if(B->order == C->order) {
            T zero = 0.0;
            T one = 1.0;
            if(stage == 'B') buffersize = 0;
            // if(stage == 'P') ;
            if(stage == 'C') onesparse::gemm(h->qq, _char_to_layout(C->order), _char_to_operation(op_A), _char_to_operation(op_B), one, A->d, (T*)B->vals, C->ncols, B->ld, zero, (T*)C->vals, C->ld);
        }
        else {
            descr_matrix_dense B_compl = std::make_shared<_descr_matrix_dense>(B->get_complementary());
            char op_B_compl = mgm::operation_combine(op_B, 'T');
            mm<T,I>(h, op_A, op_B_compl, A, B_compl, C, buffersize, buffer, stage);
        }
    }

}
}
}

#include "gpu/gpu_spblas.inst.hpp"

#endif
