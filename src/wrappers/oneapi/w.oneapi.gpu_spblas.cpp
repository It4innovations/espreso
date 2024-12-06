
#ifdef HAVE_ONEAPI

#include "gpu/gpu_spblas.h"
#include "w.oneapi.gpu_management.h"
#include "basis/utilities/utils.h"
#include <oneapi/mkl.hpp>
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-local-typedef"
#include <oneapi/dpl/execution>
#pragma clang diagnostic pop
#include <oneapi/dpl/async>



namespace espreso {
namespace gpu {
namespace spblas {

    namespace onemkl = oneapi::mkl;
    namespace onesparse = onemkl::sparse;
    namespace onedpl = oneapi::dpl;

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
                    T val = sparse_vals[i];
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
        int64_t ncols = -1;
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

    struct _descr_sparse_trsv
    {
        descr_matrix_csr d_matrix;
    };

    struct _descr_sparse_trsm
    {
        descr_matrix_csr d_matrix;
    };

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
    void descr_matrix_csr_create(handle & /*h*/, descr_matrix_csr & descr, I nrows, I ncols, I nnz, char fill)
    {
        descr = std::make_shared<_descr_matrix_csr>();
        onesparse::init_matrix_handle(&descr->d);
        descr->nrows = nrows;
        descr->ncols = ncols;
        descr->nnz = nnz;
        descr->fill = fill;
    }

    template<typename T, typename I>
    void descr_matrix_csr_link_data_internal(handle & h, descr_matrix_csr & descr, I * rowptrs, I * colidxs, T * vals)
    {
        onesparse::set_csr_data(h->qq, descr->d, descr->nrows, descr->ncols, onemkl::index_base::zero, rowptrs, colidxs, vals);
        descr->rowptrs = rowptrs;
        descr->colidxs = colidxs;
        descr->vals = vals;
        onesparse::set_matrix_property(descr->d, onesparse::property::sorted);
    }

    template<typename T, typename I, typename A>
    void descr_matrix_csr_link_data(handle & h, descr_matrix_csr & descr, Matrix_CSR<T,I,A> & matrix)
    {
        descr_matrix_csr_link_data_internal(h, descr, matrix.rows, matrix.cols, matrix.vals);
    }

    void descr_matrix_csr_destroy(handle & h, descr_matrix_csr & descr)
    {
        if(descr.get() == nullptr) return;

        onesparse::release_matrix_handle(h->qq, &descr->d);
        descr.reset();
    }

    template<typename T, typename I>
    void descr_matrix_dense_create(handle & /*h*/, descr_matrix_dense & descr, I nrows, I ncols, I ld, char order)
    {
        descr = std::make_shared<_descr_matrix_dense>();
        descr->nrows = nrows;
        descr->ncols = ncols;
        descr->ld = ld;
        descr->order = order;
    }

    template<typename T, typename I, typename A>
    void descr_matrix_dense_link_data(handle & /*h*/, descr_matrix_dense & descr, Matrix_Dense<T,I,A> & matrix)
    {
        descr->vals = matrix.vals;
    }

    void descr_matrix_dense_destroy(handle & /*h*/, descr_matrix_dense & descr)
    {
        if(descr.get() == nullptr) return;

        descr.reset();
    }

    template<typename T, typename I>
    void descr_vector_dense_create(handle & /*h*/, descr_vector_dense & descr, I /*size*/)
    {
        descr = std::make_shared<_descr_vector_dense>();
    }

    template<typename T, typename I, typename A>
    void descr_vector_dense_link_data(handle & /*h*/, descr_vector_dense & descr, Vector_Dense<T,I,A> & vector)
    {
        descr->vals = vector.vals;
    }

    template<typename T, typename I, typename A>
    void descr_vector_dense_link_data(handle & /*h*/, descr_vector_dense & descr, Matrix_Dense<T,I,A> & matrix, I colidx)
    {
        descr->vals = matrix.vals + colidx * matrix.get_ld();
    }

    void descr_vector_dense_destroy(handle & /*h*/, descr_vector_dense & descr)
    {
        if(descr.get() == nullptr) return;

        descr.reset();
    }

    void descr_sparse_trsv_create(handle & /*h*/, descr_sparse_trsv & descr)
    {
        descr = std::make_shared<_descr_sparse_trsv>();
    }

    void descr_sparse_trsv_destroy(handle & h, descr_sparse_trsv & descr)
    {
        if(descr.get() == nullptr) return;

        descr_matrix_csr_destroy(h, descr->d_matrix);
        descr.reset();
    }

    void descr_sparse_trsm_create(handle & /*h*/, descr_sparse_trsm & descr)
    {
        descr = std::make_shared<_descr_sparse_trsm>();
    }

    void descr_sparse_trsm_destroy(handle & h, descr_sparse_trsm & descr)
    {
        if(descr.get() == nullptr) return;

        descr_matrix_csr_destroy(h, descr->d_matrix);
        descr.reset();
    }

    void descr_sparse_mv_create(handle & /*h*/, descr_sparse_mv & /*descr*/) {}

    void descr_sparse_mv_destroy(handle & /*h*/, descr_sparse_mv & /*descr*/) {}

    place get_place_trsm()
    {
        return place::OUT_OF_PLACE;
    }

    template<typename T, typename I>
    void transpose(handle & h, descr_matrix_csr & output, descr_matrix_csr & input, bool conjugate, size_t & buffersize, void * buffer, char stage)
    {
        auto exec_pol = onedpl::execution::make_device_policy(h->qq);

        if(stage == 'A') {
            buffersize = 0;
        }
        if(stage == 'B') {
            buffersize = 0;
            buffersize += input->nnz * sizeof(I); // map
            buffersize += input->nnz * sizeof(I); // output_ijv_rowidxs
        }
        if(stage == 'P') {
            I * map = (I*)buffer;
            buffer = (char*)buffer + input->nnz * sizeof(I);
            I * output_ijv_rowidxs = (I*)buffer;
            I * output_ijv_colidxs = (I*)output->colidxs;
            I * input_ijv_rowidxs = output_ijv_colidxs;
            const I * input_ijv_colidxs = (I*)input->colidxs;

            auto counting_begin = onedpl::counting_iterator<I>(0);
            onedpl::experimental::copy_async(exec_pol, counting_begin, counting_begin + input->nnz, map);

            {
                I * csr_rowptrs = (I*)input->rowptrs;
                I * ijv_rowidxs = input_ijv_rowidxs;
                int wpw = 64;
                int wpk = input->nrows;
                h->qq.parallel_for(
                    sycl::nd_range<1>(sycl::range<1>(wpk*wpw), sycl::range<1>(wpw)),
                    [=](sycl::nd_item<1> item) {
                        I row = item.get_group_linear_id();
                        I start = csr_rowptrs[row];
                        I end = csr_rowptrs[row+1];
                        for(I i = start + item.get_local_linear_id(); i < end; i += item.get_group().get_local_linear_range()) {
                            ijv_rowidxs[i] = row;
                        }
                    }
                );
            }

            h->qq.copy<I>(input_ijv_colidxs, output_ijv_rowidxs, input->nnz);

            // cant use stable_sort_by_key or sort_by_key, because asynchronicity is not specified
            auto zip_iter = onedpl::make_zip_iterator(output_ijv_rowidxs, output_ijv_colidxs, map);
            onedpl::experimental::sort_async(exec_pol, zip_iter, zip_iter + input->nnz, [](auto l, auto r){
                // if(onedpl::get<0>(l) == onedpl::get<0>(r)) {
                //     return onedpl::get<1>(l) < onedpl::get<1>(r);
                // }
                // return onedpl::get<0>(l) < onedpl::get<0>(r);

                auto [rowl, coll, mapl] = l;
                auto [rowr, colr, mapr] = r;
                if(rowl == rowr) {
                    return coll < colr;
                }
                return rowl < rowr;
            });

            {
                I * ijv_rowidxs = output_ijv_rowidxs;
                I * csr_rowptrs = (I*)output->rowptrs;
                I nrows = output->nrows;
                I nnz = input->nnz;
                I wpw = 256;
                I wpk = 64;
                h->qq.parallel_for(
                    sycl::nd_range<1>(sycl::range<1>(wpk*wpw), sycl::range<1>(wpw)),
                    [=](sycl::nd_item<1> item) {
                        for(I i = item.get_global_linear_id(); i < nnz-1; i += item.get_global_range(0)) {
                            I curr_row = ijv_rowidxs[i];
                            I next_row = ijv_rowidxs[i+1];
                            for(I r = curr_row; r < next_row; r++) csr_rowptrs[r+1] = i+1;
                        }
                        if(item.get_global_linear_id() == item.get_global_range(0)-1) {
                            I lastrow = ijv_rowidxs[nnz-1];
                            for(I r = lastrow; r < nrows; r++) csr_rowptrs[r+1] = nnz;
                            I firstrow = ijv_rowidxs[0];
                            for(I r = 0; r <= firstrow; r++) csr_rowptrs[r] = 0;
                        }
                    }
                );
            }
        }
        if(stage == 'C') {
            I * map = (I*)buffer;
            auto perm_iter = onedpl::make_permutation_iterator((T*)input->vals, map);
            if(utils::is_complex<T>() && conjugate) {
                auto perm_conj_iter = onedpl::make_transform_iterator(perm_iter, [](T x){ if constexpr(utils::is_complex<T>()){ utils::imag_ref(x) *= -1; } return x; });
                onedpl::experimental::copy_async(exec_pol, perm_conj_iter, perm_conj_iter + output->nnz, (T*)output->vals);
            }
            else {
                onedpl::experimental::copy_async(exec_pol, perm_iter, perm_iter + output->nnz, (T*)output->vals);
            }
        }
    }

    template<typename T, typename I>
    void sparse_to_dense(handle & h, char op, descr_matrix_csr & sparse, descr_matrix_dense & dense, size_t & buffersize, void * buffer, char stage)
    {
        if(stage == 'A') {
            buffersize = 0;
        }
        // if(stage == 'B') ;
        // if(stage == 'P') ;
        if(stage == 'C') {
            if(op == 'N') {
                I avg_nnz_per_row = sparse->nnz / sparse->nrows;
                int workitems_in_workgroup = std::clamp(avg_nnz_per_row, 64, 1024);
                sycl::nd_range<1> range(sycl::range<1>(sparse->nrows * workitems_in_workgroup), sycl::range<1>(workitems_in_workgroup));
                if(dense->order == 'R') {
                    h->qq.fill(dense->vals, T{0}, dense->nrows * dense->ld);
                    h->qq.parallel_for(
                        range,
                        kernel_functor_sp2dn<T,I,'R'>((I*)sparse->rowptrs, (I*)sparse->colidxs, (T*)sparse->vals, (T*)dense->vals, dense->ld)
                    );
                }
                else if(dense->order == 'C') {
                    h->qq.fill(dense->vals, T{0}, dense->ncols * dense->ld);
                    h->qq.parallel_for(
                        range,
                        kernel_functor_sp2dn<T,I,'C'>((I*)sparse->rowptrs, (I*)sparse->colidxs, (T*)sparse->vals, (T*)dense->vals, dense->ld)
                    );
                }
                else {
                    eslog::error("invalid order '%c'\n", dense->order);
                }
            }
            else if(op == 'C') {
                sparse_to_dense<T,I>(h, op, sparse, dense, buffersize, buffer, stage);
                if constexpr(utils::is_complex<T>()) {
                    utils::remove_complex_t<T> mone = -1;
                    utils::remove_complex_t<T> * ptr = (utils::remove_complex_t<T>*)sparse->vals + 1;
                    if(stage == 'C') onemkl::blas::row_major::scal(h->qq, sparse->nnz, mone, ptr, 2);
                }
            }
            else {
                char op_compl = mgm::operation_combine(op, 'T');
                descr_matrix_dense descr_dense_compl = std::make_shared<_descr_matrix_dense>(dense->get_complementary());
                sparse_to_dense<T,I>(h, op_compl, sparse, descr_dense_compl, buffersize, buffer, stage);
            }
        }
    }


    template<typename T, typename I>
    void trsv(handle & h, char op, descr_matrix_csr & matrix, descr_vector_dense & rhs, descr_vector_dense & sol, descr_sparse_trsv & descr_trsv, buffer_sizes & buffersizes, void * /*buffer_persistent*/, void * /*buffer_tmp*/, char stage)
    {
        if(op != 'N') {
            eslog::error("gpu::spblas::trsv in oneapi: only 'N' is supported for parameter op\n");
        }

        T one = 1.0;
        if(stage == 'A') {
            // no idea and hard to figure out, just guess and hope
            buffersizes.allocsize_internal = ((gpu::mgm::operation_remove_conj(transpose) == 'N') ? (0) : (8 * matrix->nnz));
        }
        if(stage == 'B') {
            buffersizes.persistent = 0;
            buffersizes.tmp_preprocess = 0;
            buffersizes.tmp_update = 0;
            buffersizes.tmp_compute = 0;
        }
        // if(stage == 'P') ;
        if(stage == 'U') {
            // data in sparse matrix cannot change, according to the matrix handle contract. So I have to create a new matrix handle, and possibly preprocess again
            // the operation descriptor contains the matrix handle
            descr_matrix_csr_destroy(h, descr_trsv->d_matrix);
            descr_matrix_csr_create<T,I>(h, descr_trsv->d_matrix, matrix->nrows, matrix->ncols, matrix->nnz, matrix->fill);
            descr_matrix_csr_link_data_internal<T,I>(h, descr_trsv->d_matrix, (I*)matrix->rowptrs, (I*)matrix->colidxs, (T*)matrix->vals);
            onesparse::optimize_trsv(h->qq, _char_to_uplofill(matrix->fill), _char_to_operation(op), onemkl::diag::nonunit, descr_trsv->d_matrix->d);
        }
        if(stage == 'C') onesparse::trsv(h->qq, _char_to_uplofill(matrix->fill), _char_to_operation(op), onemkl::diag::nonunit, one, descr_trsv->d_matrix->d, (T*)rhs->vals, (T*)sol->vals);
    }

    template<typename T, typename I>
    void trsm(handle & h, char op_mat, char op_rhs, char op_sol, descr_matrix_csr & matrix, descr_matrix_dense & rhs, descr_matrix_dense & sol, descr_sparse_trsm & descr_trsm, buffer_sizes & buffersizes, void * buffer_persistent, void * buffer_tmp, char stage)
    {
        if(rhs.get() == sol.get()) eslog::error("wrong rhs and sol parameters: must not be the same, because trsm in oneapi::sparse is out-of-place\n");

        if(op_sol == 'N') {
            if(rhs->order == sol->order) {
                if(op_mat != 'N') {
                    eslog::error("gpu::spblas::trsm in oneapi: only 'N' is supported for parameter op_mat\n");
                }
                if(op_rhs != 'N') {
                    eslog::error("gpu::spblas::trsm in oneapi: only 'N' is supported for parameter op_rhs\n");
                }

                T one = 1.0;
                if(stage == 'A') {
                    // no idea and hard to figure out, just guess and hope
                    buffersizes.allocsize_internal = 8 * matrix->nnz;
                }
                if(stage == 'B') {
                    buffersizes.persistent = 0;
                    buffersizes.tmp_preprocess = 0;
                    buffersizes.tmp_update = 0;
                    buffersizes.tmp_compute = 0;
                }
                // if(stage == 'P') ;
                if(stage == 'U') {
                    // data in sparse matrix cannot change, according to the matrix handle contract. So I have to create a new matrix handle, and possibly preprocess again
                    // the operation descriptor contains the matrix handle
                    descr_matrix_csr_destroy(h, descr_trsm->d_matrix);
                    descr_matrix_csr_create<T,I>(h, descr_trsm->d_matrix, matrix->nrows, matrix->ncols, matrix->nnz, matrix->fill);
                    descr_matrix_csr_link_data_internal<T,I>(h, descr_trsm->d_matrix, (I*)matrix->rowptrs, (I*)matrix->colidxs, (T*)matrix->vals);
                    // onesparse::optimize_trsm(h->qq, _char_to_layout(sol->order), _char_to_uplofill(matrix->fill), _char_to_operation(op_mat), onemkl::diag::nonunit, descr_trsm->d_matrix->d, sol->ncols);
                }
                if(stage == 'C') onesparse::trsm(h->qq, _char_to_layout(sol->order), _char_to_operation(op_mat), _char_to_operation(op_rhs), _char_to_uplofill(matrix->fill), onemkl::diag::nonunit, one, descr_trsm->d_matrix->d, (T*)rhs->vals, sol->ncols, rhs->ld, (T*)sol->vals, sol->ld);
            }
            else {
                descr_matrix_dense rhs_compl = std::make_shared<_descr_matrix_dense>(rhs->get_complementary());
                char op_rhs_compl = mgm::operation_combine(op_rhs, 'T');
                trsm<T,I>(h, op_mat, op_rhs_compl, op_sol, matrix, rhs_compl, sol, descr_trsm, buffersizes, buffer_persistent, buffer_tmp, stage);
            }
        }
        else if(op_sol == 'T') {
            descr_matrix_dense sol_compl = std::make_shared<_descr_matrix_dense>(sol->get_complementary());
            char op_sol_compl = mgm::operation_combine(op_sol, 'T'); // 'N'
            trsm<T,I>(h, op_mat, op_rhs, op_sol_compl, matrix, rhs, sol_compl, descr_trsm, buffersizes, buffer_persistent, buffer_tmp, stage);
        }
        else {
            eslog::error("op_sol '%c' not supported\n", op_sol);
            // char op_mat_compl = mgm::operation_combine(op_mat, 'C');
            // char op_rhs_compl = mgm::operation_combine(op_rhs, 'C');
            // char op_sol_compl = mgm::operation_combine(op_sol, 'C');
            // trsm<T,I>(h, op_mat_compl, op_rhs_compl, op_sol_compl, matrix, rhs, sol, descr_trsm, buffersizes, buffer_persistend, buffer_tmp, stage);
        }
    }

    template<typename T, typename I>
    void mv(handle & h, char op, descr_matrix_csr & A, descr_vector_dense & x, descr_vector_dense & y, descr_sparse_mv & /*descr_mv*/, size_t & buffersize, void * /*buffer*/, char stage)
    {
        T zero = 0.0;
        T one = 1.0;
        if(stage == 'A') buffersize = 0;
        if(stage == 'B') buffersize = 0;
        if(stage == 'P') onesparse::optimize_gemv(h->qq, _char_to_operation(op), A->d);
        if(stage == 'C') onesparse::gemv(h->qq, _char_to_operation(op), one, A->d, (T*)x->vals, zero, (T*)y->vals);
    }

    template<typename T, typename I>
    void mm(handle & h, char op_A, char op_B, descr_matrix_csr & A, descr_matrix_dense & B, descr_matrix_dense & C, size_t & buffersize, void * buffer, char stage)
    {
        if(B->order == C->order) {
            if(op_B != 'N') {
                eslog::error("gpu::spblas::mm in oneapi: only 'N' is supported for parameter op_B\n");
            }

            T zero = 0.0;
            T one = 1.0;
            if(stage == 'A') buffersize = 0;
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
