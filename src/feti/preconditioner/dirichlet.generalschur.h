
#ifndef SRC_FETI_PREDONCITIONER_DIRICHLET_GENERALSCHUR_CPU_H_
#define SRC_FETI_PREDONCITIONER_DIRICHLET_GENERALSCHUR_CPU_H_

#include "preconditioner.h"
#include "math/primitives_new.h"
#include "math/wrappers/math.spsolver.h"
#include "math/operations/schur_csx_dny.h"
#include "math/operations/permute_csx_csx_map.h"
#include "math/operations/submatrix_dnx_dnx_view.h"
#include "gpu/operations/schur_hcsx_ddny.h"

namespace espreso {

template <typename T, typename I>
struct DirichletGeneralSchur: public Preconditioner<T> {
    DirichletGeneralSchur(FETI<T> &feti, char assemble_apply_where_);
    ~DirichletGeneralSchur();

    void setup() override;
    size_t get_wss_gpu_persistent() override { return wss_gpu_persistent; }
    size_t get_wss_gpu_internal() override { return wss_gpu_internal; }
    void set_ws_gpu_persistent(void * ws_gpu_persistent_) override { ws_gpu_persistent = ws_gpu_persistent_; }

    virtual void info() override;
    virtual void set(const step::Step &step) override;
    virtual void update(const step::Step &step) override;
    virtual void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) override;

protected:
    using Preconditioner<T>::feti;
private:
    using schur_impl_cpu_t = typename math::operations::schur_csx_dny<T,I>::implementation_selector;
    using schur_impl_gpu_t = typename gpu::operations::schur_hcsx_ddny<T,I>::implementation_selector;
    struct config {
        bool parallel_set;
        bool parallel_update;
        bool parallel_apply;
        bool outer_timers;
        bool inner_timers;
        bool print_config;
        char order_sc;
        char assemble_where;
        char apply_where;
        schur_impl_cpu_t schur_impl_cpu;
        schur_impl_gpu_t schur_impl_gpu;
    };
    void setup_config();
private:
    struct per_domain_data {
        size_t n_dofs_domain;
        size_t n_dofs_surface;
        size_t n_dofs_internal;
        MatrixCsxView_new<T,I> K_new;
        MatrixCsxView_new<T,I> B_new;
        MatrixCsxData_new<T,I> Kperm;
        MatrixDenseView_new<T> sc;
        math::operations::submatrix_dnx_dnx_view<T> op_sc_sub_from_allocd;
        PermutationData_new<I> perm_surface_to_botright;
        VectorDenseData_new<I> map_domain_to_surface;
        math::operations::permute_csx_csx_map<T,I> op_perm_K;
        std::unique_ptr<math::operations::schur_csx_dny<T,I>> op_sc;
        std::unique_ptr<gpu::operations::schur_hcsx_ddny<T,I>> op_sc_gpu;
        VectorDenseData_new<T> apply_w;
        VectorDenseData_new<T> apply_z;
        VectorDenseData_new<T> apply_w_gpu;
        VectorDenseData_new<T> apply_z_gpu;
    };
private:
    char assemble_apply_where;
    config cfg;
    size_t n_domains;
    std::unique_ptr<AllocatorArena_new> ator_ws_gpu_persistent;
    std::vector<per_domain_data> domain_data;
    std::vector<MatrixDenseData_new<T>> sc_allocated;
    void * ws_gpu_persistent = nullptr;
    size_t wss_gpu_persistent = 0;
    size_t wss_gpu_internal = 0;
    bool are_all_hermitian = false;
    bool called_setup = false;
};

}

#endif /* SRC_FETI_PREDONCITIONER_DIRICHLET_GENERALSCHUR_CPU_H_ */
