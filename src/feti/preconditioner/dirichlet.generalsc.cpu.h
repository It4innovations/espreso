
#ifndef SRC_FETI_PREDONCITIONER_DIRICHLET_GENERALSC_CPU_H_
#define SRC_FETI_PREDONCITIONER_DIRICHLET_GENERALSC_CPU_H_

#include "preconditioner.h"
#include "math/primitives_new.h"
#include "math/wrappers/math.spsolver.h"
#include "math/operations/sc_csx_dny.h"
#include "math/operations/permute_csx_csx_map.h"

namespace espreso {

template <typename T, typename I>
struct DirichletGeneralScCpu: public Preconditioner<T> {
    DirichletGeneralScCpu(FETI<T> &feti);
    ~DirichletGeneralScCpu();

    virtual void info() override;
    virtual void set(const step::Step &step) override;
    virtual void update(const step::Step &step) override;
    virtual void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) override;

protected:
    using Preconditioner<T>::feti;
private:
    using sc_is_t = typename math::operations::sc_csx_dny<T,I>::implementation_selector;
    struct config {
        bool parallel_set = true;
        bool parallel_update = true;
        bool parallel_apply = true;
        bool outer_timers = false;
        bool inner_timers = false;
        char order_sc = 'C';
        sc_is_t sc_is = sc_is_t::autoselect;
    } cfg;
    struct per_domain_data {
        size_t n_dofs_domain;
        size_t n_dofs_surface;
        size_t n_dofs_internal;
        MatrixCsxView_new<T,I> K_new;
        MatrixCsxView_new<T,I> B_new;
        MatrixCsxData_new<T,I> Kperm;
        MatrixDenseView_new<T> sc;
        PermutationData_new<I> perm_surface_to_botright;
        VectorDenseData_new<I> map_domain_to_surface;
        math::operations::permute_csx_csx_map<T,I> op_perm_K;
        std::unique_ptr<math::operations::sc_csx_dny<T,I>> op_sc;
        VectorDenseData_new<T> apply_w;
        VectorDenseData_new<T> apply_z;
    };
    size_t n_domains;
    std::vector<per_domain_data> domain_data;
    std::vector<MatrixDenseData_new<T>> sc_allocated;
    bool are_all_hermitian = false;
};

}

#endif /* SRC_FETI_PREDONCITIONER_DIRICHLET_GENERALSC_CPU_H_ */
