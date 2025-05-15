
#include "gpu/operations/sc_hcsx_ddny.h"

#include "gpu/operations/sc_hcsx_ddny.tria.h"
#include "basis/utilities/stacktimer.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
std::unique_ptr<sc_hcsx_ddny<T,I>> sc_hcsx_ddny<T,I>::make(implementation_selector is)
{
    auto autoselect_implementation = [](){
        return implementation_selector::triangular;
    };

    switch(is) {
        case implementation_selector::autoselect:
            return sc_hcsx_ddny<T,I>::make(autoselect_implementation());
        case implementation_selector::triangular:
            return std::make_unique<sc_hcsx_ddny_tria<T,I>>();
        default:
            eslog::error("invalid implementation selector\n");
    }
}



template<typename T, typename I>
void sc_hcsx_ddny<T,I>::set_handles(gpu::mgm::queue q_, gpu::spblas::handle handle_spblas_, gpu::dnblas::handle handle_dnblas_)
{
    if(called_set_handles) eslog::error("handles are already set\n");

    q = q_;
    handle_spblas = handle_spblas_;
    handle_dnblas = handle_dnblas_;

    called_set_handles = true;
}



template<typename T, typename I>
void sc_hcsx_ddny<T,I>::set_coefficients(Treal alpha_)
{
    alpha = alpha_;
}



template<typename T, typename I>
void sc_hcsx_ddny<T,I>::set_matrix(MatrixCsxView_new<T,I> * h_A11_, MatrixCsxView_new<T,I> * h_A12_, MatrixCsxView_new<T,I> * h_A21_, MatrixCsxView_new<T,I> * h_A22_)
{
    if(called_set_matrix != '_') eslog::error("matrix is already set\n");

    h_A11 = h_A11_;
    h_A12 = h_A12_;
    h_A21 = h_A21_;
    h_A22 = h_A22_;

    called_set_matrix = '4';
}



template<typename T, typename I>
void sc_hcsx_ddny<T,I>::set_matrix(MatrixCsxView_new<T,I> * h_A_, size_t size_sc_)
{
    if(called_set_matrix != '_') eslog::error("matrix is already set\n");

    h_A = h_A_;
    size_sc = size_sc_;

    called_set_matrix = '1';
}



template<typename T, typename I>
void sc_hcsx_ddny<T,I>::set_sc(MatrixDenseView_new<T> * d_sc_)
{
    if(d_sc != nullptr) eslog::error("matrix d_sc is already set\n");

    d_sc = d_sc_;
}



template<typename T, typename I>
void sc_hcsx_ddny<T,I>::set_need_solve_A11(bool need_solve_A11_)
{
    if(called_setup) eslog::error("cannot re-set need_solve_A11 after setup was called\n");

    need_solve_A11 = need_solve_A11_;
}



template<typename T, typename I>
void sc_hcsx_ddny<T,I>::setup()
{
    stacktimer::push("sc_csx_dny::setup");

    if(called_setup) eslog::error("setup was already called\n");
    if(called_set_matrix == '_') eslog::error("matrix is not set\n");
    if(d_sc == nullptr) eslog::error("d_sc is not set\n");
    if(!d_sc->ator->is_data_accessible_gpu()) eslog::error("matrix sc must be gpu-accessible\n");
    if(d_sc->nrows != d_sc->ncols) eslog::error("d_sc has to be square\n");

    if(called_set_matrix == '4') {
        if(h_A11 != nullptr && !h_A11->ator->is_data_accessible_cpu()) eslog::error("matrix h_A11 must be cpu-accessible\n");
        if(h_A12 != nullptr && !h_A12->ator->is_data_accessible_cpu()) eslog::error("matrix h_A12 must be cpu-accessible\n");
        if(h_A21 != nullptr && !h_A21->ator->is_data_accessible_cpu()) eslog::error("matrix h_A21 must be cpu-accessible\n");
        if(h_A22 != nullptr && !h_A22->ator->is_data_accessible_cpu()) eslog::error("matrix h_A22 must be cpu-accessible\n");
        if(h_A11 == nullptr) eslog::error("h_A11 cannot be nullptr\n");
        int num_sides_set = (int)(h_A12 != nullptr) + (int)(h_A21 != nullptr);
        if(is_hermitian<T>(h_A11->prop.symm) && num_sides_set == 0) eslog::error("at least one of h_A12 and h_A21 has to be set\n");
        if(!is_hermitian<T>(h_A11->prop.symm) && num_sides_set <= 1) eslog::error("both h_A12 and h_A21 have to be set\n");
        if(h_A11->nrows != h_A11->ncols) eslog::error("h_A11 has to be square\n");
        if(h_A22 != nullptr && h_A22->nrows != h_A22->ncols) eslog::error("h_A22 has to be square\n");
        if(h_A22 != nullptr && h_A22->prop.symm != h_A11->prop.symm) eslog::error("h_A11 and h_A22 must have equal symmetry\n");
        if(!is_symmetry_equal<T>(h_A11->prop.symm, d_sc->prop.symm)) eslog::error("h_A11 and d_sc must have equal symmetry\n");
        if(h_A12 != nullptr && h_A12->nrows != h_A11->nrows) eslog::error("wrong matrix h_A12 size (does not match h_A11)\n");
        if(h_A21 != nullptr && h_A21->ncols != h_A11->ncols) eslog::error("wrong matrix h_A21 size (does not match h_A11)\n");
        if(h_A12 != nullptr && h_A22 != nullptr && h_A12->ncols != h_A22->ncols) eslog::error("wrong matrix h_A12 size (does not natch h_A22)\n");
        if(h_A21 != nullptr && h_A22 != nullptr && h_A21->nrows != h_A22->nrows) eslog::error("wrong matrix h_A21 size (does not natch h_A22)\n");
        // if only one of h_A12 or h_A21 is set and the system is hermitian, it is assumed the other is its conjugate transpose
        // if h_A22 is nullptr, it is considered to be a null matrix

        if(h_A12 != nullptr) size_sc = h_A12->ncols;
        if(h_A21 != nullptr) size_sc = h_A21->nrows;
        size_A11 = h_A11->nrows;
        size_matrix = size_A11 + size_sc;

        is_matrix_hermitian = is_hermitian<T>(h_A11->prop.symm);

        if(is_matrix_hermitian && h_A11->prop.uplo != 'U' && h_A11->prop.uplo != 'L') eslog::error("wrong h_A11 uplo\n");
        if(is_matrix_hermitian && d_sc->prop.uplo != 'U' && d_sc->prop.uplo != 'L') eslog::error("wrong sc uplo\n");
        if(is_matrix_hermitian && h_A22 != nullptr && h_A22->prop.uplo != 'U' && h_A22->prop.uplo != 'L') eslog::error("wrong h_A22 uplo\n");
    }

    if(called_set_matrix == '1') {
        if(!h_A->ator->is_data_accessible_cpu()) eslog::error("matrix h_A must be cpu-accessible\n");

        size_matrix = h_A->nrows;
        size_A11 = size_matrix - size_sc;

        if(h_A->prop.symm != d_sc->prop.symm) eslog::error("h_A and d_sc must have equal symmetry\n");

        is_matrix_hermitian = is_hermitian<T>(h_A->prop.symm);
        
        if(is_matrix_hermitian && h_A->prop.uplo != 'U' && h_A->prop.uplo != 'L') eslog::error("wrong h_A uplo\n");
        if(is_matrix_hermitian && d_sc->prop.uplo != 'U' && d_sc->prop.uplo != 'L') eslog::error("wrong d_sc uplo\n");
    }

    if(d_sc->nrows != size_sc) eslog::error("mismatch between provided size_sc and SC matrix size\n");

    ator_ws_persistent = std::make_unique<AllocatorArena_new>(AllocatorGPU_new::get_singleton());
    ator_ws_tmp_linear = std::make_unique<AllocatorArena_new>(AllocatorGPU_new::get_singleton());
    ator_ws_tmp_overlap = std::make_unique<AllocatorSinglePointer_new>(AllocatorGPU_new::get_singleton());

    this->internal_setup();

    called_setup = true;

    stacktimer::pop();
}



template<typename T, typename I>
size_t sc_hcsx_ddny<T,I>::get_wss_internal()
{
    return wss_internal;
}



template<typename T, typename I>
size_t sc_hcsx_ddny<T,I>::get_wss_persistent()
{
    return wss_persistent;
}



template<typename T, typename I>
size_t sc_hcsx_ddny<T,I>::get_wss_tmp_preprocess()
{
    return wss_tmp_preprocess;
}



template<typename T, typename I>
size_t sc_hcsx_ddny<T,I>::get_wss_tmp_perform()
{
    return wss_tmp_perform;
}



template<typename T, typename I>
void sc_hcsx_ddny<T,I>::set_ws_persistent(void * ws_persistent_)
{
    if(ws_persistent_ == nullptr && wss_persistent > 0) eslog::error("persistent workspace is null\n");
    if(ws_persistent != nullptr) eslog::error("cannot re-set persistent workspace\n");

    ws_persistent = ws_persistent_;
}



template<typename T, typename I>
void sc_hcsx_ddny<T,I>::preprocess_submit(void * ws_tmp)
{
    stacktimer::push("sc_hcsx_ddny::preprocess_submit");

    if(!called_setup) eslog::error("setup has not been called\n");
    if(called_preprocess) eslog::error("preprocess has already been called\n");
    if(ws_tmp == nullptr && wss_tmp_preprocess > 0) eslog::error("temporary workspace is null\n");

    ator_ws_persistent->set(ws_persistent, wss_persistent);

    ator_ws_tmp_linear->set(ws_tmp, wss_tmp_preprocess_linear);
    ws_tmp = utils::pointer_advance(ws_tmp, wss_tmp_preprocess_linear);
    ator_ws_tmp_overlap->set(ws_tmp, wss_tmp_preprocess_overlap);

    this->internal_preprocess_submit();

    ator_ws_tmp_linear->unset();
    ator_ws_tmp_overlap->unset();

    called_preprocess = true;

    stacktimer::pop();
}



template<typename T, typename I>
void sc_hcsx_ddny<T,I>::perform_1_submit()
{
    stacktimer::push("sc_hcsx_ddny::perform_1_submit");

    if(!called_preprocess) eslog::error("preprocess has not been called\n");

    this->internal_perform_1_submit();

    stacktimer::pop();
}



template<typename T, typename I>
void sc_hcsx_ddny<T,I>::perform_2_submit(void * ws_tmp)
{
    stacktimer::push("sc_hcsx_ddny::perform_2_submit");

    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(ws_tmp == nullptr && wss_tmp_perform > 0) eslog::error("temporary workspace is null\n");

    ator_ws_tmp_linear->set(ws_tmp, wss_tmp_perform_linear);
    ws_tmp = utils::pointer_advance(ws_tmp, wss_tmp_perform_linear);
    ator_ws_tmp_overlap->set(ws_tmp, wss_tmp_perform_overlap);

    this->internal_perform_2_submit();

    ator_ws_tmp_linear->unset();
    ator_ws_tmp_overlap->unset();

    called_perform = true;

    stacktimer::pop();
}



template<typename T, typename I>
void sc_hcsx_ddny<T,I>::solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol)
{
    stacktimer::push("sc_hcsx_ddny::solve_A11");

    if(!called_preprocess) eslog::error("preprocess has not been called\n");
    if(!need_solve_A11) eslog::error("need_solve_A11 is not set, so cannot solve A11\n");
    if(rhs.size != size_A11) eslog::error("wrong rhs size\n");
    if(sol.size != size_A11) eslog::error("wrong sol size\n");

    this->internal_solve_A11(rhs, sol);

    stacktimer::pop();
}



#define INSTANTIATE_T_I(T,I) \
template class sc_hcsx_ddny<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T,int32_t) \
    /* INSTANTIATE_T_I(T,int64_t) */

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
