
#include "feti/preconditioner/dirichlet.generalsc.cpu.h"

#include "math/operations/mv_dnx.h"
#include "math/operations/sorting_permutation.h"
#include "math/primitives_new/allocator_new.h"
#include "basis/utilities/stacktimer.h"


namespace espreso {



template<typename T, typename I>
DirichletGeneralScCpu<T,I>::DirichletGeneralScCpu(FETI<T> &feti)
    : Preconditioner<T>(feti)
{
}



template<typename T, typename I>
DirichletGeneralScCpu<T,I>::~DirichletGeneralScCpu()
{
}



template<typename T, typename I>
void DirichletGeneralScCpu<T,I>::info()
{
}



template<typename T, typename I>
void DirichletGeneralScCpu<T,I>::set(const step::Step &/*step*/)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("DirichletGeneralScCpu::set");

    n_domains = feti.K.size();

    domain_data.resize(n_domains);

    stacktimer::push("DirichletGeneralScCpu::set setup");
    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_set)
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_data & my = domain_data[di];

        stacktimer::info("subdomain %zu", di);

        my.K_new = MatrixCsxView_new<T,I>::from_old(feti.K[di]);
        my.B_new = MatrixCsxView_new<T,I>::from_old(feti.B1[di]);
        if(my.B_new.order != 'R') eslog::error("B_new is assumed to have row major order\n");

        my.n_dofs_domain = my.K_new.nrows;

        my.perm_surface_to_botright.set(my.n_dofs_domain, AllocatorCPU_new::get_singleton());
        my.perm_surface_to_botright.alloc();

        {
            VectorDenseData_new<I> surface_indices;
            surface_indices.set(my.n_dofs_domain, AllocatorCPU_new::get_singleton());
            surface_indices.alloc();
            std::fill_n(surface_indices.vals, surface_indices.size, I{0});
            for(size_t i = 0; i < my.B_new.nnz; i++) {
                I col = my.B_new.idxs[i];
                surface_indices.vals[col] = I{1};
            }

            my.n_dofs_surface = std::count_if(surface_indices.vals, surface_indices.vals + surface_indices.size, [](I x){ return x != 0; });
            my.n_dofs_internal = my.n_dofs_domain - my.n_dofs_surface;
            
            math::operations::sorting_permutation<I,I>::do_all(&surface_indices, &my.perm_surface_to_botright);
        }

        my.map_domain_to_surface.set(my.n_dofs_domain, AllocatorCPU_new::get_singleton());
        my.map_domain_to_surface.alloc();
        std::fill_n(my.map_domain_to_surface.vals, my.map_domain_to_surface.size, std::numeric_limits<I>::max());
        for(size_t i = 0; i < my.n_dofs_domain; i++) {
            my.map_domain_to_surface.vals[i] = my.perm_surface_to_botright.src_to_dst[i] - my.n_dofs_internal;
        }

        // TODO: what if K has only upper triangle?
        my.Kperm.set(my.n_dofs_domain, my.n_dofs_domain, my.K_new.nnz, my.K_new.order, AllocatorCPU_new::get_singleton());
        my.Kperm.prop = my.K_new.prop;
        my.Kperm.alloc();

        my.op_perm_K.set_matrix_src(&my.K_new);
        my.op_perm_K.set_matrix_dst(&my.Kperm);
        my.op_perm_K.set_perm_rows(&my.perm_surface_to_botright);
        my.op_perm_K.set_perm_cols(&my.perm_surface_to_botright);

        my.op_sc = math::operations::schur_csx_dny<T,I>::make(cfg.sc_is);
        my.op_sc->set_coefficients(utils::remove_complex_t<T>{1});
        my.op_sc->set_matrix(&my.Kperm, my.n_dofs_surface);
        my.op_sc->set_sc(&my.sc);
        my.op_sc->set_need_solve_A11(false);

        my.apply_w.set(my.n_dofs_surface, AllocatorCPU_new::get_singleton());
        my.apply_z.set(my.n_dofs_surface, AllocatorCPU_new::get_singleton());

        my.apply_w.alloc();
        my.apply_z.alloc();
    }
    if(!cfg.inner_timers) stacktimer::enable();
    stacktimer::pop();

    are_all_hermitian = std::all_of(domain_data.begin(), domain_data.end(), [](const per_domain_data & my) { return is_hermitian<T>(my.K_new.prop.symm); });

    if(are_all_hermitian) {
        sc_allocated.resize((n_domains - 1) / 2 + 1);
        std::vector<size_t> domain_idxs_sorted_by_sc_size_desc(n_domains);
        for(size_t di = 0; di < n_domains; di++) {
            domain_idxs_sorted_by_sc_size_desc[di] = di;
        }
        std::sort(domain_idxs_sorted_by_sc_size_desc.rbegin(), domain_idxs_sorted_by_sc_size_desc.rend(), [&](size_t dl, size_t dr){ return domain_data[dl].n_dofs_surface < domain_data[dr].n_dofs_surface; });
        for(size_t i = 0; i < n_domains; i++) {
            size_t di = domain_idxs_sorted_by_sc_size_desc[i];
            size_t di_bigger = domain_idxs_sorted_by_sc_size_desc[(i / 2) * 2];
            size_t allocated_sc_index = i / 2;
            MatrixDenseData_new<T> & sc_allocd = sc_allocated[allocated_sc_index];
            per_domain_data & data_bigger = domain_data[di_bigger];
            per_domain_data & data_di = domain_data[di];
            if(di == di_bigger) {
                sc_allocd.set(data_bigger.n_dofs_surface + (cfg.order_sc == 'R'), data_bigger.n_dofs_surface + (cfg.order_sc == 'C'), cfg.order_sc, AllocatorCPU_new::get_singleton());
                sc_allocd.alloc();
                data_di.sc.set_view(data_di.n_dofs_surface, data_di.n_dofs_surface, sc_allocd.ld, sc_allocd.order, sc_allocd.vals + sc_allocd.ld, sc_allocd.ator);
                data_di.sc.prop.uplo = (sc_allocd.order == 'R') ? 'L' : 'U';
            }
            else {
                data_di.sc = sc_allocd.get_submatrix_view(0, data_di.n_dofs_surface, 0, data_di.n_dofs_surface);
                data_di.sc.prop.uplo = (sc_allocd.order == 'R') ? 'U' : 'L';
            }
            data_di.sc.prop.symm = MatrixSymmetry_new::hermitian;
        }
    }
    else {
        sc_allocated.resize(n_domains);
        for(size_t di = 0; di < n_domains; di++) {
            per_domain_data & my = domain_data[di];
            sc_allocated[di].set(my.n_dofs_surface, my.n_dofs_surface, cfg.order_sc, AllocatorCPU_new::get_singleton());
            my.sc = sc_allocated[di];
        }
    }

    stacktimer::push("DirichletGeneralScCpu::set preprocess");
    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_set)
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_data & my = domain_data[di];

        stacktimer::info("subdomain %zu", di);

        my.op_perm_K.perform_pattern();

        my.op_sc->preprocess();
    }
    if(!cfg.inner_timers) stacktimer::enable();
    stacktimer::pop();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void DirichletGeneralScCpu<T,I>::update(const step::Step &/*step*/)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("DirichletGeneralScCpu::update");

    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_set)
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_data & my = domain_data[di];

        stacktimer::info("subdomain %zu", di);

        my.op_perm_K.perform_values();

        my.op_sc->perform();
    }
    if(!cfg.inner_timers) stacktimer::enable();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void DirichletGeneralScCpu<T,I>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    stacktimer::push("DirichletGeneralScCpu::apply");

    std::fill_n(y.vals, y.size, T{0});
    
    stacktimer::push("apply_all_subdomains\n");
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_apply)
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_data & my = domain_data[di];
        std::vector<int> & myD2C = feti.D2C[di];

        std::fill_n(my.apply_w.vals, my.apply_w.size, T{0});

        // x -> select domain dofs from dual vector, apply Bt, select surface dofs -> apply_w
        for(size_t r = 0; r < my.B_new.nrows; r++) {
            I start = my.B_new.ptrs[r];
            I end = my.B_new.ptrs[r+1];
            for(I i = start; i < end; i++) {
                I c = my.B_new.idxs[i];
                T v = my.B_new.vals[i];
                I idx_x = myD2C[r];
                I idx_w = my.map_domain_to_surface.vals[c];

                my.apply_w.vals[idx_w] += v * x.vals[idx_x];
            }
        }

        // apply_w -> apply SC -> apply_z
        math::operations::mv_dnx<T>::do_all(&domain_data[di].sc, &my.apply_w, &my.apply_z, T{1}, T{0});

        // apply_z -> extend surface to domain, apply B, distribute to dual vector -> y
        for(size_t r = 0; r < my.B_new.nrows; r++) {
            I start = my.B_new.ptrs[r];
            I end = my.B_new.ptrs[r+1];
            for(I i = start; i < end; i++) {
                I c = my.B_new.idxs[i];
                T v = my.B_new.vals[i];
                I idx_z = my.map_domain_to_surface.vals[c];
                I idx_y = myD2C[r];

                utils::atomic_add<T>(y.vals[idx_y], v * my.apply_z.vals[idx_z]);
            }
        }
    }
    stacktimer::pop();
    
    y.synchronize();

    stacktimer::pop();
}



#define INSTANTIATE_T_I(T,I) \
template class DirichletGeneralScCpu<T,I>;

    #define INSTANTIATE_T(T) \
    INSTANTIATE_T_I(T, int32_t) \
    /* INSTANTIATE_T_I(T, int64_t) */

        /* INSTANTIATE_T(float) */ \
        INSTANTIATE_T(double) \
        /* INSTANTIATE_T(std::complex<float>) */ \
        INSTANTIATE_T(std::complex<double>)

    #undef INSTANTIATE_T
#undef INSTANTIATE_T_I



}
