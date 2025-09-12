
#ifndef SRC_WRAPPERS_ROCM_OPERATIONS_W_ROCM_SUBMATRIX_DCSX_DCSX_H
#define SRC_WRAPPERS_ROCM_OPERATIONS_W_ROCM_SUBMATRIX_DCSX_DCSX_H

#include "gpu/operations/submatrix_dcsx_dcsx.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_ROCM



template<typename T, typename I>
class w_rocm_submatrix_dcsx_dcsx : public submatrix_dcsx_dcsx<T,I>
{
public:
    w_rocm_submatrix_dcsx_dcsx() = default;
    virtual ~w_rocm_submatrix_dcsx_dcsx() = default;
protected:
    using submatrix_dcsx_dcsx<T,I>::q;
    using submatrix_dcsx_dcsx<T,I>::primary_start;
    using submatrix_dcsx_dcsx<T,I>::primary_end;
    using submatrix_dcsx_dcsx<T,I>::secdary_start;
    using submatrix_dcsx_dcsx<T,I>::secdary_end;
    using submatrix_dcsx_dcsx<T,I>::M_src;
    using submatrix_dcsx_dcsx<T,I>::M_dst;
    using submatrix_dcsx_dcsx<T,I>::ws_persistent;
    using submatrix_dcsx_dcsx<T,I>::wss_internal;
    using submatrix_dcsx_dcsx<T,I>::wss_persistent;
    using submatrix_dcsx_dcsx<T,I>::wss_tmp_preprocess;
    using submatrix_dcsx_dcsx<T,I>::wss_tmp_perform;
protected:
    void internal_setup() override;
    void internal_preprocess(void * ws_tmp) override;
    void internal_perform(void * ws_tmp) override;
private:
    size_t wss_pers_startptrs = 0;
    size_t wss_pers_endptrs = 0;
    size_t wss_pers_outptrs = 0;
    size_t wss_scan = 0;
    I * src_start_ptrs = nullptr;
    I * src_end_ptrs = nullptr;
    I * dst_ptrs = nullptr;
};



#else



template<typename T, typename I>
class w_rocm_submatrix_dcsx_dcsx : public submatrix_dcsx_dcsx<T,I>
{
public:
    w_rocm_submatrix_dcsx_dcsx() { eslog::error("rocm wrapper is not available\n"); }
    virtual ~w_rocm_submatrix_dcsx_dcsx() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_ROCM_OPERATIONS_W_ROCM_SUBMATRIX_DCSX_DCSX_H */
