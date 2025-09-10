
#ifndef SRC_GPU_OPERATIONS_SCHUR_HCSX_DDNY_MANUAL_SIMPLE_H
#define SRC_GPU_OPERATIONS_SCHUR_HCSX_DDNY_MANUAL_SIMPLE_H

#include "gpu/operations/schur_hcsx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



template<typename T, typename I>
struct schur_hcsx_ddny_manual_simple_data;

template<typename T, typename I>
class schur_hcsx_ddny_manual_simple : public schur_hcsx_ddny<T,I>
{
public:
    schur_hcsx_ddny_manual_simple();
    schur_hcsx_ddny_manual_simple(const schur_hcsx_ddny_manual_simple &) = delete;
    schur_hcsx_ddny_manual_simple(schur_hcsx_ddny_manual_simple &&) = default;
    schur_hcsx_ddny_manual_simple & operator=(const schur_hcsx_ddny_manual_simple &) = delete;
    schur_hcsx_ddny_manual_simple & operator=(schur_hcsx_ddny_manual_simple &&) = default;
    virtual ~schur_hcsx_ddny_manual_simple();
public:
    const char * get_name() override { return "schur_hcsx_ddny_manual_simple"; }
protected:
    void internal_setup() override;
    void internal_preprocess_submit() override;
    void internal_perform_1_submit() override;
    void internal_perform_2_submit() override;
    void internal_solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol) override;
private:
    std::unique_ptr<schur_hcsx_ddny_manual_simple_data<T,I>> data;
protected:
    using schur_hcsx_ddny<T,I>::q;
    using schur_hcsx_ddny<T,I>::handle_spblas;
    using schur_hcsx_ddny<T,I>::handle_dnblas;
    using schur_hcsx_ddny<T,I>::h_A11;
    using schur_hcsx_ddny<T,I>::h_A12;
    using schur_hcsx_ddny<T,I>::h_A21;
    using schur_hcsx_ddny<T,I>::h_A22;
    using schur_hcsx_ddny<T,I>::h_A;
    using schur_hcsx_ddny<T,I>::d_sc;
    using schur_hcsx_ddny<T,I>::alpha;
    using schur_hcsx_ddny<T,I>::need_solve_A11;
    using schur_hcsx_ddny<T,I>::called_set_handles;
    using schur_hcsx_ddny<T,I>::called_set_matrix;
    using schur_hcsx_ddny<T,I>::called_setup;
    using schur_hcsx_ddny<T,I>::called_preprocess;
    using schur_hcsx_ddny<T,I>::called_perform;
    using schur_hcsx_ddny<T,I>::ws_persistent;
    using schur_hcsx_ddny<T,I>::wss_internal;
    using schur_hcsx_ddny<T,I>::wss_persistent;
    using schur_hcsx_ddny<T,I>::wss_tmp_preprocess;
    using schur_hcsx_ddny<T,I>::wss_tmp_perform;
    using schur_hcsx_ddny<T,I>::size_matrix;
    using schur_hcsx_ddny<T,I>::size_sc;
    using schur_hcsx_ddny<T,I>::size_A11;
    using schur_hcsx_ddny<T,I>::is_matrix_hermitian;
    using schur_hcsx_ddny<T,I>::ator_ws_persistent;
    using schur_hcsx_ddny<T,I>::ator_ws_tmp_linear;
    using schur_hcsx_ddny<T,I>::ator_ws_tmp_overlap;
    using schur_hcsx_ddny<T,I>::wss_tmp_preprocess_linear;
    using schur_hcsx_ddny<T,I>::wss_tmp_preprocess_overlap;
    using schur_hcsx_ddny<T,I>::wss_tmp_perform_linear;
    using schur_hcsx_ddny<T,I>::wss_tmp_perform_overlap;
};



}
}
}




#endif /* SRC_GPU_OPERATIONS_SCHUR_HCSX_DDNY_MANUAL_SIMPLE_H */
