
#ifndef SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEAPI_SUPERMATRIX_DDNX_DDNX_H
#define SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEAPI_SUPERMATRIX_DDNX_DDNX_H

#include "gpu/operations/supermatrix_ddnx_ddnx_noncontig.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_ONEAPI



template<typename T, typename I>
class w_oneapi_supermatrix_ddnx_ddnx_noncontig : public supermatrix_ddnx_ddnx_noncontig<T,I>
{
public:
    w_oneapi_supermatrix_ddnx_ddnx_noncontig() = default;
    virtual ~w_oneapi_supermatrix_ddnx_ddnx_noncontig() = default;
protected:
    using mode = typename supermatrix_ddnx_ddnx_noncontig<T,I>::mode;
    using supermatrix_ddnx_ddnx_noncontig<T,I>::q;
    using supermatrix_ddnx_ddnx_noncontig<T,I>::d_M_src;
    using supermatrix_ddnx_ddnx_noncontig<T,I>::d_M_dst;
    using supermatrix_ddnx_ddnx_noncontig<T,I>::d_row_map;
    using supermatrix_ddnx_ddnx_noncontig<T,I>::d_col_map;
    using supermatrix_ddnx_ddnx_noncontig<T,I>::mode_val;
protected:
    void internal_perform() override;
};



#else



template<typename T, typename I>
class w_oneapi_supermatrix_ddnx_ddnx_noncontig : public supermatrix_ddnx_ddnx_noncontig<T,I>
{
public:
    w_oneapi_supermatrix_ddnx_ddnx_noncontig() { eslog::error("oneapi wrapper is not available\n"); }
    virtual ~w_oneapi_supermatrix_ddnx_ddnx_noncontig() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEAPI_SUPERMATRIX_DDNX_DDNX_H */
