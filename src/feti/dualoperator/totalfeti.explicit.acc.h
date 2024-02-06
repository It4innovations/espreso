
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

	static constexpr size_t align_B = 512;
	static constexpr size_t align_elem = align_B / sizeof(T);
	int stage = 0;
	char wcpset, wcpupdate, wcpapply, spdnfactor, trsvtrsm, factrs1, factrs2syrk, handlecount, applyalg;
	gpu::mgm::device device;
	size_t n_domains;
	size_t n_queues;
	void * mem_pool_device;
	std::unique_ptr<cbmba_resource> cbmba_res_device;
	gpu::mgm::queue main_q;
	std::vector<gpu::mgm::queue> queues;
	std::vector<gpu::dnblas::handle> handles_dense;
	std::vector<gpu::spblas::handle> handles_sparse;
	std::vector<DirectSparseSolver<T>> solvers_Kreg;
	std::vector<Matrix_Dense<T,I,Ad>> d_Fs;
	std::vector<Matrix_CSR<T,I>> Kregs;
	std::vector<Matrix_CSR<T,I,Ah>> h_Us_sp;
	std::vector<Matrix_CSR<T,I,Ah>> h_Ls_sp;
	std::vector<Matrix_CSR<T,I,Ah>> h_Bperms_sp;
	std::vector<Matrix_CSR<T,I,Ad>> d_Us_sp;
	std::vector<Matrix_CSR<T,I,Ad>> d_Ls_sp;
	std::vector<Matrix_CSR<T,I,Ad>> d_Bperms_sp;
	std::vector<std::unique_ptr<Matrix_Dense<T,I,cbmba_d>>> d_Us_dn;
	std::vector<std::unique_ptr<Matrix_Dense<T,I,cbmba_d>>> d_Ls_dn;
	std::vector<std::unique_ptr<Matrix_Dense<T,I,cbmba_d>>> d_Xs_r;
	std::vector<std::unique_ptr<Matrix_Dense<T,I,cbmba_d>>> d_Xs_c;
	std::vector<std::unique_ptr<Matrix_Dense<T,I,cbmba_d>>> d_Ys_r;
	std::vector<std::unique_ptr<Matrix_Dense<T,I,cbmba_d>>> d_Ys_c;
	std::vector<gpu::spblas::descr_matrix_dense> descr_Fs_r;
	std::vector<gpu::spblas::descr_matrix_dense> descr_Fs_c;
	std::vector<gpu::spblas::descr_matrix_csr> descr_Us_sp1;
	std::vector<gpu::spblas::descr_matrix_csr> descr_Us_sp2;
	std::vector<gpu::spblas::descr_matrix_csr> descr_Ls_sp1;
	std::vector<gpu::spblas::descr_matrix_csr> descr_Ls_sp2;
	std::vector<gpu::spblas::descr_matrix_csr> descr_Bperms_sp;
	std::vector<gpu::spblas::descr_matrix_dense> descr_Us_dn;
	std::vector<gpu::spblas::descr_matrix_dense> descr_Ls_dn;
	std::vector<gpu::spblas::descr_matrix_dense> descr_Xs_r;
	std::vector<gpu::spblas::descr_matrix_dense> descr_Xs_c;
	std::vector<gpu::spblas::descr_matrix_dense> descr_Ys_r;
	std::vector<gpu::spblas::descr_matrix_dense> descr_Ys_c;
	std::vector<std::vector<gpu::spblas::descr_vector_dense>> descr_Xs_vecs;
	std::vector<std::vector<gpu::spblas::descr_vector_dense>> descr_Ys_vecs;
	std::vector<gpu::spblas::descr_sparse_trsv> descrs_sparse_trsv1;
	std::vector<gpu::spblas::descr_sparse_trsv> descrs_sparse_trsv2;
	std::vector<gpu::spblas::descr_sparse_trsm> descrs_sparse_trsm1;
	std::vector<gpu::spblas::descr_sparse_trsm> descrs_sparse_trsm2;
	std::vector<size_t> buffersizes_sptrs1;
	std::vector<size_t> buffersizes_sptrs2;
	std::vector<size_t> buffersizes_spmm;
	std::vector<size_t> buffersizes_other;
	std::vector<void*> buffers_sptrs1;
	std::vector<void*> buffers_sptrs2;
	std::vector<void*> buffers_spmm;
	std::vector<Vector_Dense<I,I>> transmaps_L2U;
	std::vector<Vector_Dense<I,I>> transmaps_U2L;
	std::vector<Vector_Dense<I,I,Ad>> d_applyg_D2Cs;
	std::vector<Vector_Dense<T,I,Ad>> d_apply_xs;
	std::vector<Vector_Dense<T,I,Ad>> d_apply_ys;
	std::vector<Vector_Dense<T,I,Ah>> h_applyc_xs;
	std::vector<Vector_Dense<T,I,Ah>> h_applyc_ys;
	Vector_Dense<T,I,Ad> d_applyg_x_cluster;
	Vector_Dense<T,I,Ad> d_applyg_y_cluster;
	Vector_Dense<T*,I,Ad> d_applyg_xs_pointers;
	Vector_Dense<T*,I,Ad> d_applyg_ys_pointers;
	Vector_Dense<I,I,Ad> d_applyg_n_dofs_interfaces;
	Vector_Dense<I*,I,Ad> d_applyg_D2Cs_pointers;

};

}

#endif /* SRC_FETI_DUALOPERATOR_TOTALFETI_EXPLICIT_ACC_H_ */
