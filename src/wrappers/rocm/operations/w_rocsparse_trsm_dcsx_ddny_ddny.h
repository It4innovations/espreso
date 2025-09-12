
#ifndef SRC_WRAPPERS_ROCM_OPERATIONS_W_ROCSPARSE_TRSM_DCSX_DDNY_DDNY_H
#define SRC_WRAPPERS_ROCM_OPERATIONS_W_ROCSPARSE_TRSM_DCSX_DDNY_DDNY_H

#include "gpu/operations/trsm_dcsx_ddny_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_ROCM



struct w_rocsparse_trsm_dcsx_ddny_ddny_data;

template<typename T, typename I>
class w_rocsparse_trsm_dcsx_ddny_ddny : public trsm_dcsx_ddny_ddny<T,I>
{
public:
    w_rocsparse_trsm_dcsx_ddny_ddny();
    virtual ~w_rocsparse_trsm_dcsx_ddny_ddny();
protected:
    using trsm_dcsx_ddny_ddny<T,I>::q;
    using trsm_dcsx_ddny_ddny<T,I>::handle_spblas;
    using trsm_dcsx_ddny_ddny<T,I>::A;
    using trsm_dcsx_ddny_ddny<T,I>::X;
    using trsm_dcsx_ddny_ddny<T,I>::B;
    using trsm_dcsx_ddny_ddny<T,I>::ws_persistent;
    using trsm_dcsx_ddny_ddny<T,I>::wss_internal;
    using trsm_dcsx_ddny_ddny<T,I>::wss_persistent;
    using trsm_dcsx_ddny_ddny<T,I>::wss_tmp_preprocess;
    using trsm_dcsx_ddny_ddny<T,I>::wss_tmp_perform;
    using trsm_dcsx_ddny_ddny<T,I>::place;
protected:
    char internal_get_native_place() override;
    void internal_setup() override;
    void internal_preprocess(void * ws_tmp) override;
    void internal_perform(void * ws_tmp) override;
private:
    std::unique_ptr<w_rocsparse_trsm_dcsx_ddny_ddny_data> data;
};



#else



template<typename T, typename I>
class w_rocsparse_trsm_dcsx_ddny_ddny : public trsm_dcsx_ddny_ddny<T,I>
{
public:
    w_rocsparse_trsm_dcsx_ddny_ddny() { eslog::error("rocm wrapper is not available\n"); }
    virtual ~w_rocsparse_trsm_dcsx_ddny_ddny() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_ROCM_OPERATIONS_W_ROCSPARSE_TRSM_DCSX_DDNY_DDNY_H */
