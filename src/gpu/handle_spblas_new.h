
#ifndef SRC_GPU_HANDLE_SPBLAS_NEW_H
#define SRC_GPU_HANDLE_SPBLAS_NEW_H



namespace espreso {
namespace gpu {



class handle_spblas_new
{
public:
    handle_spblas_new() = default;
    handle_spblas_new(const handle_spblas_new &) = delete;
    handle_spblas_new(handle_spblas_new &&) = delete;
    handle_spblas_new & operator=(const handle_spblas_new &) = delete;
    handle_spblas_new & operator=(handle_spblas_new &&) = delete;
    virtual ~handle_spblas_new() = default;
    static std::unique_ptr<handle_spblas_new> make();
};



}
}

#endif /* SRC_GPU_HANDLE_SPBLAS_NEW_H */
