
#ifndef SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEAPI_SUBMATRIX_DCSX_DDNY_H
#define SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEAPI_SUBMATRIX_DCSX_DDNY_H

#include "gpu/operations/submatrix_dcsx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_ONEAPI



template<typename T, typename I>
class w_oneapi_submatrix_dcsx_ddny : public submatrix_dcsx_ddny<T,I>
{
public:
    w_oneapi_submatrix_dcsx_ddny() = default;
    virtual ~w_oneapi_submatrix_dcsx_ddny() = default;
protected:
    using submatrix_dcsx_ddny<T,I>::q;
    using submatrix_dcsx_ddny<T,I>::primary_start;
    using submatrix_dcsx_ddny<T,I>::primary_end;
    using submatrix_dcsx_ddny<T,I>::secdary_start;
    using submatrix_dcsx_ddny<T,I>::secdary_end;
    using submatrix_dcsx_ddny<T,I>::M_src;
    using submatrix_dcsx_ddny<T,I>::M_dst;
    using submatrix_dcsx_ddny<T,I>::ws_persistent;
    using submatrix_dcsx_ddny<T,I>::wss_internal;
    using submatrix_dcsx_ddny<T,I>::wss_persistent;
    using submatrix_dcsx_ddny<T,I>::wss_tmp_preprocess;
    using submatrix_dcsx_ddny<T,I>::wss_tmp_perform;
protected:
    void internal_setup() override;
    void internal_preprocess(void * ws_tmp) override;
    void internal_perform(void * ws_tmp) override;
};



#else



template<typename T, typename I>
class w_oneapi_submatrix_dcsx_ddny : public submatrix_dcsx_ddny<T,I>
{
public:
    w_oneapi_submatrix_dcsx_ddny() { eslog::error("oneapi wrapper is not available\n"); }
    virtual ~w_oneapi_submatrix_dcsx_ddny() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_ONEAPI_OPERATIONS_W_ONEAPI_SUBMATRIX_DCSX_DDNY_H */
