
#ifdef HAVE_ROCM

#include "gpu/gpu_spblas.h"
#include "w.rocm.gpu_management.h"

#include <rocsparse.h>



inline void _check(rocsparse_status status, const char *file, int line)
{
    if (status != rocsparse_status_success)
    {
        espreso::eslog::error("ROCSPARSE Error %d. In file '%s' on line %d\n", status, file, line);
    }
}



namespace espreso {
namespace gpu {
namespace spblas {

    namespace
    {
        template<typename I>
        static rocsparse_indextype _sparse_index_type()
        {
            if constexpr(std::is_same_v<I, int32_t>) return rocsparse_indextype_i32;
            if constexpr(std::is_same_v<I, int64_t>) return rocsparse_indextype_i64;
        }

        template<typename T>
        static rocsparse_datatype _sparse_data_type()
        {
            if constexpr(std::is_same_v<T, float>)  return rocsparse_datatype_f32_r;
            if constexpr(std::is_same_v<T, double>) return rocsparse_datatype_f64_r;
            if constexpr(std::is_same_v<T, std::complex<float>>)  return rocsparse_datatype_f32_c;
            if constexpr(std::is_same_v<T, std::complex<double>>) return rocsparse_datatype_f64_c;
        }

        static rocsparse_operation _char_to_operation(char c)
        {
            switch(c)
            {
                case 'N': return rocsparse_operation_none;
                case 'T': return rocsparse_operation_transpose;
                case 'H': return rocsparse_operation_conjugate_transpose;
                default: eslog::error("invalid operation '%c'\n", c);
            }
        }
    }

    struct _handle
    {
        rocsparse_handle h;
    };

    struct _descr_matrix_csr
    {
        rocsparse_spmat_descr d;
    };

    struct _descr_matrix_dense
    {
        rocsparse_dnmat_descr d;
        rocsparse_dnmat_descr d_complementary;
        char order;
        _descr_matrix_dense get_complementary()
        {
            _descr_matrix_dense ret;
            ret.d = d_complementary;
            ret.d_complementary = d;
            ret.order = mgm::change_order(order);
            return ret;
        }
    };

    struct _descr_vector_dense
    {
        rocsparse_dnvec_descr d;
    };

    struct _descr_sparse_trsv { };

    struct _descr_sparse_trsm { };

    void handle_create(handle & h, mgm::queue & q)
    {
        h = std::make_shared<_handle>();
        CHECK(rocsparse_create_handle(&h->h));
        CHECK(rocsparse_set_stream(h->h, q->stream));
    }

    void handle_destroy(handle & h)
    {
        if(h.get() == nullptr) return;

        CHECK(rocsparse_destroy_handle(h->h));
        h.reset();
    }

    template<typename T, typename I>
    void descr_matrix_csr_create(descr_matrix_csr & descr, I nrows, I ncols, I nnz, char fill)
    {
        descr = std::make_shared<_descr_matrix_csr>();
        void * dummyptr = reinterpret_cast<void*>(1);
        CHECK(rocsparse_create_csr_descr(&descr->d, nrows, ncols, nnz, dummyptr, dummyptr, dummyptr, _sparse_index_type<I>(), _sparse_index_type<I>(), rocsparse_index_base_zero, _sparse_data_type<T>()));
        
        rocsparse_fill_mode upper = rocsparse_fill_mode_upper;
        rocsparse_fill_mode lower = rocsparse_fill_mode_lower;
        rocsparse_diag_type nonunit = rocsparse_diag_type_non_unit;
        // rocsparse_matrix_type triangular = rocsparse_matrix_type_triangular;
        rocsparse_storage_mode sorted = rocsparse_storage_mode_sorted;
        if(fill == 'L') CHECK(rocsparse_spmat_set_attribute(descr->d, rocsparse_spmat_fill_mode, &lower, sizeof(lower)));
        if(fill == 'U') CHECK(rocsparse_spmat_set_attribute(descr->d, rocsparse_spmat_fill_mode, &upper, sizeof(upper)));
        if(fill != 'N') CHECK(rocsparse_spmat_set_attribute(descr->d, rocsparse_spmat_diag_type, &nonunit, sizeof(nonunit)));
        // if(fill != 'N') CHECK(rocsparse_spmat_set_attribute(descr->d, rocsparse_spmat_matrix_type, &triangular, sizeof(triangular)));
        CHECK(rocsparse_spmat_set_attribute(descr->d, rocsparse_spmat_storage_mode, &sorted, sizeof(sorted)));
    }

    template<typename T, typename I, typename A>
    void descr_matrix_csr_link_data(descr_matrix_csr & descr, Matrix_CSR<T,I,A> & matrix)
    {
        static_assert(A::is_data_device_accessible, "matrix data must be device accessible");
        CHECK(rocsparse_csr_set_pointers(descr->d, matrix.rows, matrix.cols, matrix.vals));
    }

    void descr_matrix_csr_destroy(descr_matrix_csr & descr)
    {
        if(descr.get() == nullptr) return;

        CHECK(rocsparse_destroy_spmat_descr(descr->d));
        descr.reset();
    }

    template<typename T, typename I>
    void descr_matrix_dense_create(descr_matrix_dense & descr, I nrows, I ncols, I ld, char order)
    {
        descr = std::make_shared<_descr_matrix_dense>();
        void * dummyptr = reinterpret_cast<void*>(1);
        if(order == 'R') CHECK(rocsparse_create_dnmat_descr(&descr->d, nrows, ncols, ld, dummyptr, _sparse_data_type<T>(), rocsparse_order_row));
        if(order == 'C') CHECK(rocsparse_create_dnmat_descr(&descr->d, nrows, ncols, ld, dummyptr, _sparse_data_type<T>(), rocsparse_order_column));
        if(order == 'R') CHECK(rocsparse_create_dnmat_descr(&descr->d_complementary, ncols, nrows, ld, dummyptr, _sparse_data_type<T>(), rocsparse_order_column));
        if(order == 'C') CHECK(rocsparse_create_dnmat_descr(&descr->d_complementary, ncols, nrows, ld, dummyptr, _sparse_data_type<T>(), rocsparse_order_row));
        descr->order = order;
    }

    template<typename T, typename I, typename A>
    void descr_matrix_dense_link_data(descr_matrix_dense & descr, Matrix_Dense<T,I,A> & matrix)
    {
        static_assert(A::is_data_device_accessible, "matrix data must be device accessible");
        CHECK(rocsparse_dnmat_set_values(descr->d, matrix.vals));
        CHECK(rocsparse_dnmat_set_values(descr->d_complementary, matrix.vals));
    }

    void descr_matrix_dense_destroy(descr_matrix_dense & descr)
    {
        if(descr.get() == nullptr) return;

        CHECK(rocsparse_destroy_dnmat_descr(descr->d));
        CHECK(rocsparse_destroy_dnmat_descr(descr->d_complementary));
        descr.reset();
    }

    template<typename T, typename I>
    void descr_vector_dense_create(descr_vector_dense & descr, I size)
    {
        descr = std::make_shared<_descr_vector_dense>();
        void * dummyptr = reinterpret_cast<void*>(1);
        CHECK(rocsparse_create_dnvec_descr(&descr->d, size, dummyptr, _sparse_data_type<T>()));
    }

    template<typename T, typename I, typename A>
    void descr_vector_dense_link_data(descr_vector_dense & descr, Vector_Dense<T,I,A> & vector)
    {
        static_assert(A::is_data_device_accessible, "matrix data must be device accessible");
        CHECK(rocsparse_dnvec_set_values(descr->d, vector.vals));
    }

    template<typename T, typename I, typename A>
    void descr_vector_dense_link_data(descr_vector_dense & descr, Matrix_Dense<T,I,A> & matrix, I colidx)
    {
        static_assert(A::is_data_device_accessible, "vector data must be device accessible");
        CHECK(rocsparse_dnvec_set_values(descr->d, matrix.vals + colidx * matrix.get_ld()));
    }

    void descr_vector_dense_destroy(descr_vector_dense & descr)
    {
        if(descr.get() == nullptr) return;

        CHECK(rocsparse_destroy_dnvec_descr(descr->d));
        descr.reset();
    }

    void descr_sparse_trsv_create(descr_sparse_trsv & /*descr*/) { }

    void descr_sparse_trsv_destroy(descr_sparse_trsv & /*descr*/) { }

    void descr_sparse_trsm_create(descr_sparse_trsm & /*descr*/) { }

    void descr_sparse_trsm_destroy(descr_sparse_trsm & /*descr*/) { }

    template<typename T, typename I>
    void sparse_to_dense(handle & h, char transpose, descr_matrix_csr & sparse, descr_matrix_dense & dense, size_t & buffersize, void * buffer, char stage)
    {
        if(transpose == 'N')
        {
            if(stage == 'B') CHECK(rocsparse_sparse_to_dense(h->h, sparse->d, dense->d, rocsparse_sparse_to_dense_alg_default, &buffersize, nullptr));
            if(stage == 'C') CHECK(rocsparse_sparse_to_dense(h->h, sparse->d, dense->d, rocsparse_sparse_to_dense_alg_default, &buffersize, buffer));
        }
        else if(transpose == 'T')
        {
            descr_matrix_dense descr_dense_complementary = std::make_shared<_descr_matrix_dense>(dense->get_complementary());
            sparse_to_dense<T,I>(h, 'N', sparse, descr_dense_complementary, buffersize, buffer, stage);
        }
        else
        {
            eslog::error("transpose '%c' not supported\n", transpose);
        }
    }

    template<typename T, typename I>
    void trsv(handle & h, char transpose, descr_matrix_csr & matrix, descr_vector_dense & rhs, descr_vector_dense & sol, descr_sparse_trsv & /*descr_trsv*/, size_t & buffersize, void * buffer, char stage)
    {
        T one = 1.0;
        if(stage == 'B') CHECK(rocsparse_spsv(h->h, _char_to_operation(transpose), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), rocsparse_spsv_alg_default, rocsparse_spsv_stage_buffer_size, &buffersize, buffer));
        if(stage == 'P') CHECK(rocsparse_spsv(h->h, _char_to_operation(transpose), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), rocsparse_spsv_alg_default, rocsparse_spsv_stage_preprocess,  &buffersize, buffer));
        // if(stage == 'U') ; // no update matrix function, hopefully dont need to to anything. otherwise redo preprocessing
        if(stage == 'C') CHECK(rocsparse_spsv(h->h, _char_to_operation(transpose), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), rocsparse_spsv_alg_default, rocsparse_spsv_stage_compute,     &buffersize, buffer));
    }

    template<typename T, typename I>
    void trsm(handle & h, char transpose_mat, char transpose_rhs, char transpose_sol, descr_matrix_csr & matrix, descr_matrix_dense & rhs, descr_matrix_dense & sol, descr_sparse_trsm & descr_trsm, size_t & buffersize, void * buffer, char stage)
    {
#if HIP_VERSION_MAJOR >= 6
        static_assert(false, "not sure how (not) buggy it is in newer versions, check and redo");
#else
        // older rocsparse assumes colmajor for both dense matrices
        // older rocsparse mistakenly thinks that trans_B is applied to both B and C
        if(rhs->order == sol->order)
        {
            if(transpose_rhs == transpose_sol)
            {
                if(rhs->order == 'C')
                {
                    T one = 1.0;
                    if(stage == 'B') CHECK(rocsparse_spsm(h->h, _char_to_operation(transpose_mat), _char_to_operation(transpose_rhs), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), rocsparse_spsm_alg_default, rocsparse_spsm_stage_buffer_size, &buffersize, buffer));
                    if(stage == 'P') CHECK(rocsparse_spsm(h->h, _char_to_operation(transpose_mat), _char_to_operation(transpose_rhs), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), rocsparse_spsm_alg_default, rocsparse_spsm_stage_preprocess,  &buffersize, buffer));
                    // if(stage == 'U') ; // no update matrix function, hopefully dont need to to anything. otherwise redo preprocessing
                    if(stage == 'C') CHECK(rocsparse_spsm(h->h, _char_to_operation(transpose_mat), _char_to_operation(transpose_rhs), &one, matrix->d, rhs->d, sol->d, _sparse_data_type<T>(), rocsparse_spsm_alg_default, rocsparse_spsm_stage_compute,     &buffersize, buffer));
                }
                else if(rhs->order == 'R')
                {
                    descr_matrix_dense descr_rhs_compl = std::make_shared<_descr_matrix_dense>(rhs->get_complementary());
                    descr_matrix_dense descr_sol_compl = std::make_shared<_descr_matrix_dense>(sol->get_complementary());
                    char transpose_compl = mgm::change_operation_array_transpose(transpose_rhs);
                    trsm<T,I>(h, transpose_mat, transpose_compl, transpose_compl, matrix, descr_rhs_compl, descr_sol_compl, descr_trsm, buffersize, buffer, stage);
                }
                else
                {
                    eslog::error("wrong dense matrix order '%c'\n", rhs->order);
                }
            }
            else
            {
                eslog::error("unsupported combimation of matrix orders and transpositions '%c'\n", rhs->order);
            }
        }
        else
        {
            descr_matrix_dense descr_rhs_compl = std::make_shared<_descr_matrix_dense>(rhs->get_complementary());
            char transpose_rhs_compl = mgm::change_operation_array_transpose(transpose_rhs);
            trsm<T,I>(h, transpose_mat, transpose_rhs_compl, transpose_sol, matrix, descr_rhs_compl, sol, descr_trsm, buffersize, buffer, stage);
        }
#endif
    }

    template<typename T, typename I>
    void mv(handle & h, char transpose, descr_matrix_csr & A, descr_vector_dense & x, descr_vector_dense & y, size_t & buffersize, void * buffer, char stage)
    {
        T one = 1.0;
        T zero = 0.0;
// should check for ROCSPARSE_VERSION_*, but I don't understand their versioning
#if HIP_VERSION_MAJOR >= 6 // rocparse that comes with hip/rocm 6.0.x, replaced old spmv with new spmv with stage (aka renamed the previous spmv_ex)
        if(stage == 'B') CHECK(rocsparse_spmv(h->h, _char_to_operation(transpose), &one, A->d, x->d, &zero, y->d, _sparse_data_type<T>(), rocsparse_spmv_alg_default, rocsparse_spmv_stage_buffer_size, &buffersize, buffer));
        if(stage == 'P') CHECK(rocsparse_spmv(h->h, _char_to_operation(transpose), &one, A->d, x->d, &zero, y->d, _sparse_data_type<T>(), rocsparse_spmv_alg_default, rocsparse_spmv_stage_preprocess,  &buffersize, buffer));
        if(stage == 'C') CHECK(rocsparse_spmv(h->h, _char_to_operation(transpose), &one, A->d, x->d, &zero, y->d, _sparse_data_type<T>(), rocsparse_spmv_alg_default, rocsparse_spmv_stage_compute,     &buffersize, buffer));
#elif HIP_VERSION_MAJOR == 5 && HIP_VERSION_MINOR >= 4 // rocsparse that comes with hip/rocm 5.4.x, deprecated original spmv and added spmv_ex with stage
        if(stage == 'B') CHECK(rocsparse_spmv_ex(h->h, _char_to_operation(transpose), &one, A->d, x->d, &zero, y->d, _sparse_data_type<T>(), rocsparse_spmv_alg_default, rocsparse_spmv_stage_buffer_size, &buffersize, buffer));
        if(stage == 'P') CHECK(rocsparse_spmv_ex(h->h, _char_to_operation(transpose), &one, A->d, x->d, &zero, y->d, _sparse_data_type<T>(), rocsparse_spmv_alg_default, rocsparse_spmv_stage_preprocess,  &buffersize, buffer));
        if(stage == 'C') CHECK(rocsparse_spmv_ex(h->h, _char_to_operation(transpose), &one, A->d, x->d, &zero, y->d, _sparse_data_type<T>(), rocsparse_spmv_alg_default, rocsparse_spmv_stage_compute,     &buffersize, buffer));
#else // older
        if(stage == 'B') CHECK(rocsparse_spmv(h->h, _char_to_operation(transpose), &one, A->d, x->d, &zero, y->d, _sparse_data_type<T>(), rocsparse_spmv_alg_default, &buffersize, nullptr));
        if(stage == 'C') CHECK(rocsparse_spmv(h->h, _char_to_operation(transpose), &one, A->d, x->d, &zero, y->d, _sparse_data_type<T>(), rocsparse_spmv_alg_default, &buffersize, buffer));
#endif
    }

    template<typename T, typename I>
    void mm(handle & h, char transpose_A, char transpose_B, descr_matrix_csr & A, descr_matrix_dense & B, descr_matrix_dense & C, size_t & buffersize, void * buffer, char stage)
    {
        T zero = 0.0;
        T one = 1.0;
        if(stage == 'B') CHECK(rocsparse_spmm(h->h, _char_to_operation(transpose_A), _char_to_operation(transpose_B), &one, A->d, B->d, &zero, C->d, _sparse_data_type<T>(), rocsparse_spmm_alg_default, rocsparse_spmm_stage_buffer_size, &buffersize, buffer));
        if(stage == 'P') CHECK(rocsparse_spmm(h->h, _char_to_operation(transpose_A), _char_to_operation(transpose_B), &one, A->d, B->d, &zero, C->d, _sparse_data_type<T>(), rocsparse_spmm_alg_default, rocsparse_spmm_stage_preprocess,  &buffersize, buffer));
        if(stage == 'C') CHECK(rocsparse_spmm(h->h, _char_to_operation(transpose_A), _char_to_operation(transpose_B), &one, A->d, B->d, &zero, C->d, _sparse_data_type<T>(), rocsparse_spmm_alg_default, rocsparse_spmm_stage_compute,     &buffersize, buffer));
    }

}
}
}

#include "gpu/gpu_spblas_inst.hpp"

#endif
