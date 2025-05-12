
#ifndef SRC_GPU_OPERATIONS_SC_HCSX_DDNY_TRIA_H
#define SRC_GPU_OPERATIONS_SC_HCSX_DDNY_TRIA_H

#include "gpu/operations/sc_hcsx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
struct sc_hcsx_ddny_tria_data;

template<typename T, typename I>
class sc_hcsx_ddny_tria : public sc_hcsx_ddny<T,I>
{
public:
    sc_hcsx_ddny_tria();
    sc_hcsx_ddny_tria(const sc_hcsx_ddny_tria &) = delete;
    sc_hcsx_ddny_tria(sc_hcsx_ddny_tria &&) = default;
    sc_hcsx_ddny_tria & operator=(const sc_hcsx_ddny_tria &) = delete;
    sc_hcsx_ddny_tria & operator=(sc_hcsx_ddny_tria &&) = default;
    virtual ~sc_hcsx_ddny_tria();
protected:
    virtual void internal_setup();
    virtual void internal_preprocess_submit();
    virtual void internal_perform_1_submit();
    virtual void internal_perform_2_submit();
    virtual void internal_solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol);
private:
    std::unique_ptr<sc_hcsx_ddny_tria_data<T,I>> data;
protected:
    using sc_hcsx_ddny<T,I>::q;
    using sc_hcsx_ddny<T,I>::handle_spblas;
    using sc_hcsx_ddny<T,I>::handle_dnblas;
    using sc_hcsx_ddny<T,I>::h_A11;
    using sc_hcsx_ddny<T,I>::h_A12;
    using sc_hcsx_ddny<T,I>::h_A21;
    using sc_hcsx_ddny<T,I>::h_A22;
    using sc_hcsx_ddny<T,I>::h_A;
    using sc_hcsx_ddny<T,I>::d_sc;
    using sc_hcsx_ddny<T,I>::alpha;
    using sc_hcsx_ddny<T,I>::need_solve_A11;
    using sc_hcsx_ddny<T,I>::called_set_handles;
    using sc_hcsx_ddny<T,I>::called_set_matrix;
    using sc_hcsx_ddny<T,I>::called_setup;
    using sc_hcsx_ddny<T,I>::called_preprocess;
    using sc_hcsx_ddny<T,I>::called_perform;
    using sc_hcsx_ddny<T,I>::ws_persistent;
    using sc_hcsx_ddny<T,I>::wss_internal;
    using sc_hcsx_ddny<T,I>::wss_persistent;
    using sc_hcsx_ddny<T,I>::wss_tmp_preprocess;
    using sc_hcsx_ddny<T,I>::wss_tmp_perform;
    using sc_hcsx_ddny<T,I>::size_matrix;
    using sc_hcsx_ddny<T,I>::size_sc;
    using sc_hcsx_ddny<T,I>::size_A11;
    using sc_hcsx_ddny<T,I>::is_matrix_hermitian;
    using sc_hcsx_ddny<T,I>::ator_ws_persistent;
    using sc_hcsx_ddny<T,I>::ator_ws_tmp_linear;
    using sc_hcsx_ddny<T,I>::ator_ws_tmp_overlap;
    using sc_hcsx_ddny<T,I>::wss_tmp_preprocess_linear;
    using sc_hcsx_ddny<T,I>::wss_tmp_preprocess_overlap;
    using sc_hcsx_ddny<T,I>::wss_tmp_perform_linear;
    using sc_hcsx_ddny<T,I>::wss_tmp_perform_overlap;
};



}
}
}




#endif /* SRC_GPU_OPERATIONS_SC_HCSX_DDNY_TRIA_H */
