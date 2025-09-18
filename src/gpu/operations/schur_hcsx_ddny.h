
#ifndef SRC_GPU_OPERATIONS_SC_HCSX_DDNY_H
#define SRC_GPU_OPERATIONS_SC_HCSX_DDNY_H

#include "gpu/gpu_management.h"
#include "gpu/gpu_spblas.h"
#include "gpu/gpu_dnblas.h"
#include "math/primitives_new/allocator_new.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
class schur_hcsx_ddny
{
public:
    using Treal = utils::remove_complex_t<T>;
    enum struct implementation_selector
    {
        autoselect,
        manual_simple,
        triangular
    };
protected:
    schur_hcsx_ddny() = default;
public:
    schur_hcsx_ddny(const schur_hcsx_ddny &) = delete;
    schur_hcsx_ddny(schur_hcsx_ddny &&) = default;
    schur_hcsx_ddny & operator=(const schur_hcsx_ddny &) = delete;
    schur_hcsx_ddny & operator=(schur_hcsx_ddny &&) = default;
    virtual ~schur_hcsx_ddny() = default;
public:
    static std::unique_ptr<schur_hcsx_ddny<T,I>> make(implementation_selector is = implementation_selector::autoselect);
    virtual const char * get_name() { return "UNDEFINED"; }
public:
    void set_handles(gpu::mgm::queue q_, gpu::spblas::handle spblas_handle_, gpu::dnblas::handle dnblas_handle_);
    void set_coefficients(Treal alpha_);
    void set_matrix(MatrixCsxView_new<T,I> * h_A11_, MatrixCsxView_new<T,I> * h_A12_, MatrixCsxView_new<T,I> * h_A21_, MatrixCsxView_new<T,I> * h_A22_);
    void set_matrix(MatrixCsxView_new<T,I> * h_A_, size_t size_sc_);
    void set_sc(MatrixDenseView_new<T> * d_sc_);
    void set_need_solve_A11(bool need_solve_A11_);
    void setup();
    size_t get_wss_internal();
    size_t get_wss_persistent();
    size_t get_wss_tmp_preprocess();
    size_t get_wss_tmp_perform();
    void set_ws_persistent(void * ws_persistent_);
    void preprocess_submit(void * ws_tmp);
    void perform_1_submit();
    void perform_2_submit(void * ws_tmp);
    void solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol);
    void solve_A11(MatrixDenseView_new<T> & rhs, MatrixDenseView_new<T> & sol);
protected:
    gpu::mgm::queue q;
    gpu::spblas::handle handle_spblas;
    gpu::dnblas::handle handle_dnblas;
    MatrixCsxView_new<T,I> * h_A11 = nullptr;
    MatrixCsxView_new<T,I> * h_A12 = nullptr;
    MatrixCsxView_new<T,I> * h_A21 = nullptr;
    MatrixCsxView_new<T,I> * h_A22 = nullptr;
    MatrixCsxView_new<T,I> * h_A = nullptr;
    MatrixDenseView_new<T> * d_sc = nullptr;
    Treal alpha = Treal{1};
    bool need_solve_A11 = false;
    bool called_set_handles = false;
    char called_set_matrix = '_';
    bool called_setup = false;
    bool called_preprocess = false;
    bool called_perform = false;
    void * ws_persistent = nullptr;
    size_t wss_internal = 0;
    size_t wss_persistent = 0;
    size_t wss_tmp_preprocess = 0;
    size_t wss_tmp_perform = 0;
protected:
    size_t size_matrix = 0;
    size_t size_sc = 0;
    size_t size_A11 = 0;
    bool is_matrix_hermitian = false;
protected:
    std::unique_ptr<AllocatorArena_new> ator_ws_persistent;
    std::unique_ptr<AllocatorArena_new> ator_ws_tmp_linear;
    std::unique_ptr<AllocatorSinglePointer_new> ator_ws_tmp_overlap;
    size_t wss_tmp_preprocess_linear = 0;
    size_t wss_tmp_preprocess_overlap = 0;
    size_t wss_tmp_perform_linear = 0;
    size_t wss_tmp_perform_overlap = 0;
protected:
    virtual void internal_setup() {};
    virtual void internal_preprocess_submit() {};
    virtual void internal_perform_1_submit() {};
    virtual void internal_perform_2_submit() {};
    virtual void internal_solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol) {};
    virtual void internal_solve_A11(MatrixDenseView_new<T> & rhs, MatrixDenseView_new<T> & sol) {};
};



}
}
}

#endif /* SRC_GPU_OPERATIONS_SC_HCSX_DDNY_H */
