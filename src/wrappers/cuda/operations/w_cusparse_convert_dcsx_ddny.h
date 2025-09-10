
#ifndef SRC_WRAPPERS_CUDA_OPERATIONS_W_CUSPARSE_CONVERT_DCSX_DDNY_H
#define SRC_WRAPPERS_CUDA_OPERATIONS_W_CUSPARSE_CONVERT_DCSX_DDNY_H

#include "gpu/operations/convert_dcsx_ddny.h"



namespace espreso {
namespace gpu {
namespace operations {



#ifdef HAVE_CUDA



struct w_cusparse_convert_dcsx_ddny_data;

template<typename T, typename I>
class w_cusparse_convert_dcsx_ddny : public convert_dcsx_ddny<T,I>
{
public:
    w_cusparse_convert_dcsx_ddny();
    virtual ~w_cusparse_convert_dcsx_ddny();
protected:
    using convert_dcsx_ddny<T,I>::q;
    using convert_dcsx_ddny<T,I>::handle_spblas;
    using convert_dcsx_ddny<T,I>::M_src;
    using convert_dcsx_ddny<T,I>::M_dst;
    using convert_dcsx_ddny<T,I>::ws_persistent;
    using convert_dcsx_ddny<T,I>::wss_internal;
    using convert_dcsx_ddny<T,I>::wss_persistent;
    using convert_dcsx_ddny<T,I>::wss_tmp_preprocess;
    using convert_dcsx_ddny<T,I>::wss_tmp_perform;
protected:
    void internal_setup() override;
    void internal_preprocess(void * ws_tmp) override;
    void internal_perform(void * ws_tmp) override;
private:
    std::unique_ptr<w_cusparse_convert_dcsx_ddny_data> data;
};



#else



template<typename T, typename I>
class w_cusparse_convert_dcsx_ddny : public convert_dcsx_ddny<T,I>
{
public:
    w_cusparse_convert_dcsx_ddny() { eslog::error("cuda wrapper is not available\n"); }
    virtual ~w_cusparse_convert_dcsx_ddny() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_CUDA_OPERATIONS_W_CUSPARSE_CONVERT_DCSX_DDNY_H */
