
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
    ~TotalFETIExplicitAcc();

    void info();
    void set(const step::Step &step);
    void update(const step::Step &step);

    // y = F * x
    void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);
    // y = K+(f - Bt * x)
    void toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y);

protected:
    void print(const step::Step &step);

    using DualOperator<T>::feti;
    using DualOperator<T>::d;

    using Ad = gpu::mgm::Ad;
    using Ah = gpu::mgm::Ah;

    struct per_domain_stuff {
        DirectSparseSolver<T> solver_Kreg;
        Matrix_Dense<T,I,Ad> d_F;
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
        void * buffer_sptrs1 = nullptr;
        void * buffer_sptrs2 = nullptr;
        void * buffer_spmm = nullptr;
        size_t buffersize_sptrs1;
        size_t buffersize_sptrs2;
        size_t buffersize_spmm;
        size_t buffersize_other;
        Vector_Dense<I,I> transmap_L2LH;
        Vector_Dense<I,I> transmap_U2UH;
        Vector_Dense<I,I,Ad> d_applyg_D2C;
        Vector_Dense<T,I,Ad> d_apply_x;
        Vector_Dense<T,I,Ad> d_apply_y;
        Vector_Dense<T,I,Ah> h_applyc_x;
        Vector_Dense<T,I,Ah> h_applyc_y;
    };

    DualOperatorExplicitGpuConfig * config = nullptr;
    char order_X;
    char order_F;
    char hermitian_F_fill;
    bool is_system_hermitian;
    bool is_factor1_dense, is_factor2_dense, is_factor1_sparse, is_factor2_sparse;
    bool is_path_trsm, is_path_herk;
    bool trsm1_use_L, trsm1_use_LH, trsm2_use_U, trsm2_use_UH;
    bool need_Y;
    bool solver_get_L, solver_get_U;
    bool can_use_LH_is_U_h_sp, can_use_UH_is_L_h_sp;
    bool can_use_LH_is_U_d_sp, can_use_UH_is_L_d_sp;
    bool can_use_LH_is_U_d_dn, can_use_UH_is_L_d_dn;
    bool need_conjtrans_L2LH, need_conjtrans_U2UH;
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
    Vector_Dense<T,I,Ad> d_applyg_x_cluster;
    Vector_Dense<T,I,Ad> d_applyg_y_cluster;
    Vector_Dense<T*,I,Ad> d_applyg_xs_pointers;
    Vector_Dense<T*,I,Ad> d_applyg_ys_pointers;
    Vector_Dense<I,I,Ad> d_applyg_n_dofs_interfaces;
    Vector_Dense<I*,I,Ad> d_applyg_D2Cs_pointers;

};

}

#endif /* SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_ACC_H_ */
