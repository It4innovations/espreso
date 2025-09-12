
#ifndef SRC_WRAPPERS_ROCM_OPERATIONS_W_ROCM_SUBMATRIX_DDNX_DDNX_H
#define SRC_WRAPPERS_ROCM_OPERATIONS_W_ROCM_SUBMATRIX_DDNX_DDNX_H

#include "gpu/operations/submatrix_ddnx_ddnx_noncontig.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_ROCM



template<typename T, typename I>
class w_rocm_submatrix_ddnx_ddnx_noncontig : public submatrix_ddnx_ddnx_noncontig<T,I>
{
public:
    w_rocm_submatrix_ddnx_ddnx_noncontig() = default;
    virtual ~w_rocm_submatrix_ddnx_ddnx_noncontig() = default;
protected:
    using submatrix_ddnx_ddnx_noncontig<T,I>::q;
    using submatrix_ddnx_ddnx_noncontig<T,I>::d_M_src;
    using submatrix_ddnx_ddnx_noncontig<T,I>::d_M_dst;
    using submatrix_ddnx_ddnx_noncontig<T,I>::d_row_map;
    using submatrix_ddnx_ddnx_noncontig<T,I>::d_col_map;
protected:
    void internal_perform() override;
};



#else



template<typename T, typename I>
class w_rocm_submatrix_ddnx_ddnx_noncontig : public submatrix_ddnx_ddnx_noncontig<T,I>
{
public:
    w_rocm_submatrix_ddnx_ddnx_noncontig() { eslog::error("rocm wrapper is not available\n"); }
    virtual ~w_rocm_submatrix_ddnx_ddnx_noncontig() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_ROCM_OPERATIONS_W_ROCM_SUBMATRIX_DDNX_DDNX_H */
