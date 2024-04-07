
#ifndef SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_ACC_H_
#define SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_ACC_H_

#include <vector>
#include <stdexcept>
#include <memory>

#include "totalfeti.explicit.h"
#include "math/wrappers/math.acc.feti.dual.h"
#include "gpu/gpu_management.h"
#include "gpu/gpu_dnblas.h"
#include "gpu/gpu_spblas.h"
#include "gpu/gpu_kernels.h"
#include "basis/utilities/cbmb_allocator.h"

namespace espreso {

/*
 * K+: KxK : block diagonal
 * B : LxK : from primal to dual
 *
 * y = F * x = (B * K+ * Bt) * x
 *
 * Btx = Bt * x          :: x        -> Btx     : L -> K (per domain)
 * KplusBtx = K+ * Btx   :: Btx      -> KplusBtx: K -> K (per domain)
 * y = B * KplusBtx      :: KplusBtx -> y       : K -> L (per domain except RBM and mortars)
 *
 */

template <typename T, typename I>
class TotalFETIExplicitAcc: public DualOperator<T> {
public:
    TotalFETIExplicitAcc(FETI<T> &feti);
    virtual ~TotalFETIExplicitAcc();

    virtual void info() override;
    virtual void set(const step::Step &step) override;
    virtual void update(const step::Step &step) override;

    virtual void clear_gpu_cache() override;

    // y = F * x
    virtual void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) override;
    // y = K+(f - Bt * x)
    virtual void toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y) override;

protected:
    void print(const step::Step &step);

    void config_replace_defaults();

    using DualOperator<T>::feti;
    using DualOperator<T>::d;

    using Ad = gpu::mgm::Ad;
    using Ah = gpu::mgm::Ah;

    struct per_domain_stuff {
        DirectSparseSolver<T> solver_Kreg;
        Matrix_Dense<T,I,Ad> d_F;
        std::unique_ptr<Matrix_Dense<T,I,cbmba_d>> d_F_tmp;
        Matrix_CSR<T,I> Kreg;
        Matrix_CSR<T,I,Ah> h_L_sp;
        Matrix_CSR<T,I,Ah> h_LH_sp;
        Matrix_CSR<T,I,Ah> h_U_sp;
        Matrix_CSR<T,I,Ah> h_UH_sp;
        Matrix_CSR<T,I,Ah> h_Bperm_sp;    
        Matrix_CSR<T,I,Ad> d_L_sp;
        Matrix_CSR<T,I,Ad> d_LH_sp;
        Matrix_CSR<T,I,Ad> d_U_sp;
        Matrix_CSR<T,I,Ad> d_UH_sp;
        Matrix_CSR<T,I,Ad> d_Bperm_sp;
        std::unique_ptr<Matrix_Dense<T,I,cbmba_d>> d_L_dn;
        std::unique_ptr<Matrix_Dense<T,I,cbmba_d>> d_LH_dn;
        std::unique_ptr<Matrix_Dense<T,I,cbmba_d>> d_U_dn;
        std::unique_ptr<Matrix_Dense<T,I,cbmba_d>> d_UH_dn;
        std::unique_ptr<Matrix_Dense<T,I,cbmba_d>> d_X;
        std::unique_ptr<Matrix_Dense<T,I,cbmba_d>> d_Y;
        gpu::spblas::descr_matrix_dense descr_F;
        gpu::spblas::descr_matrix_dense descr_F_tmp;
        gpu::spblas::descr_matrix_csr descr_L_sp;
        gpu::spblas::descr_matrix_csr descr_LH_sp;
        gpu::spblas::descr_matrix_csr descr_U_sp;
        gpu::spblas::descr_matrix_csr descr_UH_sp;
        gpu::spblas::descr_matrix_csr descr_Bperm_sp;
        gpu::spblas::descr_matrix_dense descr_L_dn;
        gpu::spblas::descr_matrix_dense descr_LH_dn;
        gpu::spblas::descr_matrix_dense descr_U_dn;
        gpu::spblas::descr_matrix_dense descr_UH_dn;
        gpu::spblas::descr_matrix_dense descr_X;
        gpu::spblas::descr_matrix_dense descr_Y;
        gpu::spblas::descr_sparse_trsm descr_sparse_trsm1;
        gpu::spblas::descr_sparse_trsm descr_sparse_trsm2;
        void * buffer_transL2LH = nullptr;
        void * buffer_transU2UH = nullptr;
        void * buffer_sptrs1 = nullptr;
        void * buffer_sptrs2 = nullptr;
        void * buffer_spmm = nullptr;
        size_t buffersize_transL2LH = 0;
        size_t buffersize_transU2UH = 0;
        size_t buffersize_sptrs1 = 0;
        size_t buffersize_sptrs2 = 0;
        size_t buffersize_spmm = 0;
        size_t buffersize_other = 0;
        Vector_Dense<I,I> transmap_L2LH;
        Vector_Dense<I,I> transmap_U2UH;
        Vector_Dense<I,I,Ad> d_applyg_D2C;
        Vector_Dense<T,I,Ad> d_apply_x;
        Vector_Dense<T,I,Ad> d_apply_y;
        Vector_Dense<T,I,Ah> h_applyc_x;
        Vector_Dense<T,I,Ah> h_applyc_y;
        size_t allocated_F_index;
        char hermitian_F_fill;
        I n_dofs_domain;
        I n_dofs_interface;
        I n_nz_factor;
        I ld_domain;
        I ld_interface;
        I ld_X;
        I ld_F;
        bool should_allocate_d_F;
    };

    DualOperatorExplicitGpuConfig * config = nullptr;
    char order_X, order_F;
    bool is_system_hermitian;
    bool is_factor1_dense, is_factor2_dense, is_factor1_sparse, is_factor2_sparse;
    bool is_path_trsm, is_path_herk;
    bool need_Y;
    bool solver_get_L, solver_get_U;
    bool can_use_LH_is_U_h_sp, can_use_UH_is_L_h_sp;
    bool can_use_LH_is_U_d_sp, can_use_UH_is_L_d_sp;
    bool can_use_LH_is_U_d_dn, can_use_UH_is_L_d_dn;
    bool is_f_triangles_shared, need_f_tmp;
    bool timers_basic, timers_detailed;
    bool memory_info_basic, memory_info_detailed;
    bool do_conjtrans_L2LH_h, do_conjtrans_U2UH_h, do_conjtrans_L2LH_d, do_conjtrans_U2UH_d;
    bool do_sp2dn_L, do_sp2dn_LH, do_sp2dn_U, do_sp2dn_UH;
    bool do_trsm1_sp_L, do_trsm1_sp_LH, do_trsm2_sp_U, do_trsm2_sp_UH, do_trsm1_dn_L, do_trsm1_dn_LH, do_trsm2_dn_U, do_trsm2_dn_UH;
    bool do_descr_sp_L, do_descr_sp_LH, do_descr_sp_U, do_descr_sp_UH, do_descr_dn_L, do_descr_dn_LH, do_descr_dn_U, do_descr_dn_UH;
    bool do_alloc_h_sp_L, do_alloc_h_sp_U, do_alloc_h_sp_LH, do_alloc_h_sp_UH, do_link_h_sp_LH_U, do_link_h_sp_UH_L;
    bool do_alloc_d_sp_L, do_alloc_d_sp_U, do_alloc_d_sp_LH, do_alloc_d_sp_UH, do_link_d_sp_LH_U, do_link_d_sp_UH_L;
    bool do_alloc_d_dn_L, do_alloc_d_dn_U, do_alloc_d_dn_LH, do_alloc_d_dn_UH, do_link_d_dn_LH_U, do_link_d_dn_UH_L;
    bool do_copyin_L, do_copyin_U, do_copyin_LH, do_copyin_UH;
    static constexpr size_t align_B = 512;
    static constexpr size_t align_elem = align_B / sizeof(T);
    int stage = 0;
    gpu::mgm::device device;
    size_t n_domains;
    size_t n_queues;
    void * mem_pool_device;
    std::unique_ptr<cbmba_resource> cbmba_res_device;
    gpu::mgm::queue main_q;
    std::vector<gpu::mgm::queue> queues;
    std::vector<gpu::dnblas::handle> handles_dense;
    std::vector<gpu::spblas::handle> handles_sparse;
    std::vector<per_domain_stuff> domain_data;
    std::vector<Matrix_Dense<T,I,Ad>> d_Fs_allocated;
    Vector_Dense<T,I,Ad> d_applyg_x_cluster;
    Vector_Dense<T,I,Ad> d_applyg_y_cluster;
    Vector_Dense<T*,I,Ad> d_applyg_xs_pointers;
    Vector_Dense<T*,I,Ad> d_applyg_ys_pointers;
    Vector_Dense<I,I,Ad> d_applyg_n_dofs_interfaces;
    Vector_Dense<I*,I,Ad> d_applyg_D2Cs_pointers;

};

}

#endif /* SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_ACC_H_ */
