
#ifndef SRC_FETI_DUALOPERATOR_DUALOP_EXPLICIT_APPLICATOR_H
#define SRC_FETI_DUALOPERATOR_DUALOP_EXPLICIT_APPLICATOR_H

#include "math/primitives_new/matrix_dense_data_new.h"
#include "math/primitives_new/vector_dense_data_new.h"
#include "math/primitives_new/multi_vector_dense_data_new.h"
#include "gpu/gpu_management.h"
#include "gpu/gpu_dnblas.h"
#include "feti/feti.h"
#include "math/primitives_new/allocator_new.h"


namespace espreso {



template<typename T, typename I>
class dualop_explicit_applicator
{
    // apply is without inter-process synchronization
public:
    dualop_explicit_applicator() = default;
    dualop_explicit_applicator(const dualop_explicit_applicator &) = delete;
    dualop_explicit_applicator(dualop_explicit_applicator &&) = delete;
    dualop_explicit_applicator & operator=(const dualop_explicit_applicator &) = delete;
    dualop_explicit_applicator & operator=(dualop_explicit_applicator &&) = delete;
    ~dualop_explicit_applicator() = default;
public:
    void set_config(bool timers_inner_);
    void set_handles(gpu::mgm::queue * main_q_, std::vector<gpu::mgm::queue> * queues_, std::vector<gpu::dnblas::handle> * handles_dense_);
    void set_feti(FETI<T> * feti);
    void set_memory(char vector_mem_, char Fs_mem_);
    void set_D2C_map(std::vector<std::vector<I>> * D2C_old_);
    void set_Fs(std::vector<MatrixDenseView_new<T>*> Fs_);
    void set_apply_target(char target_);
    void setup();
    size_t get_wss_gpu_persistent();
    void set_ws_gpu_persistent(void * ws_gpu_persistent_);
    void preprocess();
    void update_F(size_t di);
    void apply(VectorDenseView_new<T> & cluster_x, VectorDenseView_new<T> & cluster_y, void * ws_gpu_tmp, size_t wss_gpu_tmp, const std::function<void(void)> & func_while_waiting = [](){});
    void apply(MatrixDenseView_new<T> & cluster_X, MatrixDenseView_new<T> & cluster_Y, void * ws_gpu_tmp, size_t wss_gpu_tmp, const std::function<void(void)> & func_while_waiting = [](){});
private:
    size_t n_domains = 0;
    size_t n_queues = 0;
    size_t n_dofs_cluster_interface = 0; // size of cluster-wide dual vector
    size_t total_dofs_interface = 0; // sum of subdomain inteface sizes
    gpu::mgm::queue * main_q = nullptr;
    std::vector<gpu::mgm::queue> * queues = nullptr;
    std::vector<gpu::dnblas::handle> * handles_dnblas = nullptr;
    size_t wss_gpu_persistent = 0;
    void * ws_gpu_persistent = nullptr;
    std::unique_ptr<AllocatorArena_new> ator_ws_gpu_persistent;
    std::unique_ptr<AllocatorArena_new> ator_ws_gpu_tmp;
    std::vector<std::vector<I>> * D2C_old = nullptr;
    MultiVectorDenseData_new<I,I> D2C;
    MultiVectorDenseData_new<I,I> d_D2C;
    std::vector<size_t> n_dofs_interfaces;
    I * offsets = nullptr;
    std::vector<MatrixDenseView_new<T>*> Fs;
    std::vector<MatrixDenseData_new<T>> Fs_2;
    MultiVectorDenseData_new<T,I> xs_vec;
    MultiVectorDenseData_new<T,I> ys_vec;
    bool use_gpu = false;
    bool timers_inner = false;
    char apply_target = '_';
    char Fs_mem = '_';
    char vector_mem = '_';
    bool need_copy_vectors = false;
    bool need_copy_Fs = false;
    bool called_set_handles = false;
    bool called_setup = false;
};



}



#endif /* SRC_FETI_DUALOPERATOR_DUALOP_EXPLICIT_APPLICATOR_H */
