
#include "feti/preconditioner/dirichlet.generalschur.h"

#include "esinfo/ecfinfo.h"
#include "basis/utilities/minmaxavg.h"
#include "math/operations/mv_dnx.h"
#include "math/operations/sorting_permutation.h"
#include "math/primitives_new/allocator_new.h"
#include "math/wrappers/math.blas.h"
#include "basis/utilities/stacktimer.h"


namespace espreso {



template<typename T, typename I>
DirichletGeneralSchur<T,I>::DirichletGeneralSchur(FETI<T> &feti, char assemble_apply_where_)
    : Preconditioner<T>(feti), assemble_apply_where(assemble_apply_where_)
{
    setup_config();
}



template<typename T, typename I>
DirichletGeneralSchur<T,I>::~DirichletGeneralSchur()
{
}



template<typename T, typename I>
void DirichletGeneralSchur<T,I>::info()
{
    eslog::info(" = DIRICHLET PRECONDITIONER USING GENERAL SCHUR OPERATION ON %cPU                             = \n", assemble_apply_where);

    if(cfg.print_config) {
        auto order_to_string = [](char order){ switch(order){ case 'R': return "ROW_MAJOR"; case 'C': return "COL_MAJOR"; default: return "UNDEFINED"; }};
        auto bool_to_string = [](bool val){ return val ? "TRUE" : "FALSE";};
        auto schur_impl_cpu_to_string = [](schur_impl_cpu_t schur_impl){ switch(schur_impl) { case schur_impl_cpu_t::autoselect: return "autoselect"; case schur_impl_cpu_t::manual_simple: return "manual_simple"; case schur_impl_cpu_t::triangular: return "triangular"; case schur_impl_cpu_t::mklpardiso: return "mklpardiso"; case schur_impl_cpu_t::sparse_solver: return "sparse_solver"; case schur_impl_cpu_t::mumps: return "mumps"; case schur_impl_cpu_t::pastix: return "pastix"; default: return "UNDEFINED"; }};
        auto schur_impl_cpu_to_string_actual = [](schur_impl_cpu_t schur_impl){ return math::operations::schur_csx_dny<T,I>::make(schur_impl)->get_name(); };
        auto schur_impl_gpu_to_string = [](schur_impl_gpu_t schur_impl){ switch(schur_impl) { case schur_impl_gpu_t::autoselect: return "autoselect"; case schur_impl_gpu_t::manual_simple: return "manual_simple"; case schur_impl_gpu_t::triangular: return "triangular"; default: return "UNDEFINED"; } };
        auto schur_impl_gpu_to_string_actual = [](schur_impl_gpu_t schur_impl){ return gpu::operations::schur_hcsx_ddny<T,I>::make(schur_impl)->get_name(); };

        eslog::info(" =   %-50s       %+30s = \n", "parallel_set", bool_to_string(cfg.parallel_set));
        eslog::info(" =   %-50s       %+30s = \n", "parallel_update", bool_to_string(cfg.parallel_update));
        eslog::info(" =   %-50s       %+30s = \n", "parallel_apply", bool_to_string(cfg.parallel_apply));
        eslog::info(" =   %-50s       %+30s = \n", "inner_timers", bool_to_string(cfg.inner_timers));
        eslog::info(" =   %-50s       %+30s = \n", "outer_timers", bool_to_string(cfg.outer_timers));
        eslog::info(" =   %-50s       %+30s = \n", "order_sc", order_to_string(cfg.order_sc));
        eslog::info(" =   %-50s       %+30s = \n", "schur implementation cpu", schur_impl_cpu_to_string(cfg.schur_impl_cpu));
        eslog::info(" =   %-50s       %+30s = \n", "schur implementation cpu actual", schur_impl_cpu_to_string_actual(cfg.schur_impl_cpu));
        eslog::info(" =   %-50s       %+30s = \n", "schur implementation gpu", schur_impl_gpu_to_string(cfg.schur_impl_gpu));
        eslog::info(" =   %-50s       %+30s = \n", "schur implementation gpu actual", schur_impl_gpu_to_string_actual(cfg.schur_impl_gpu));
    }
    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_data & data){ return data.n_dofs_surface; }).to_string("  Domain surface [dofs]").c_str());
    eslog::info(minmaxavg<size_t>::compute_from_allranks(domain_data.begin(), domain_data.end(), [](const per_domain_data & data){ return data.n_dofs_internal; }).to_string("  Domain interior [dofs]").c_str());
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}



template<typename T, typename I>
void DirichletGeneralSchur<T,I>::setup()
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("DirichletGeneralSchur::setup");

    if(called_setup) eslog::error("setup was already called\n");

    n_domains = feti.K.size();

    domain_data.resize(n_domains);

    ator_ws_gpu_persistent = std::make_unique<AllocatorArena_new>(AllocatorGPU_new::get_singleton());

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

        my.Kperm.set(my.n_dofs_domain, my.n_dofs_domain, my.K_new.nnz, my.K_new.order, AllocatorCPU_new::get_singleton());
        my.Kperm.prop = my.K_new.prop;
        my.Kperm.alloc();

        my.op_perm_K.set_matrix_src(&my.K_new);
        my.op_perm_K.set_matrix_dst(&my.Kperm);
        my.op_perm_K.set_perm_rows(&my.perm_surface_to_botright);
        my.op_perm_K.set_perm_cols(&my.perm_surface_to_botright);
        my.op_perm_K.perform_pattern();

        Allocator_new * ator_wz = (assemble_apply_where == 'G') ? (Allocator_new*)AllocatorHostPinned_new::get_singleton() : (Allocator_new*)AllocatorCPU_new::get_singleton();
        my.apply_w.set(my.n_dofs_surface, ator_wz);
        my.apply_z.set(my.n_dofs_surface, ator_wz);
        if(assemble_apply_where == 'G') {
            my.apply_w_gpu.set(my.n_dofs_surface, ator_ws_gpu_persistent.get());
            my.apply_z_gpu.set(my.n_dofs_surface, ator_ws_gpu_persistent.get());
        }
    }
    if(!cfg.inner_timers) stacktimer::enable();

    are_all_hermitian = std::all_of(domain_data.begin(), domain_data.end(), [](const per_domain_data & my) { return is_hermitian<T>(my.K_new.prop.symm); });
    Allocator_new * ator_sc = ((assemble_apply_where == 'G') ? (Allocator_new*)ator_ws_gpu_persistent.get() : (Allocator_new*)AllocatorCPU_new::get_singleton());
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
                sc_allocd.set(data_bigger.n_dofs_surface + (cfg.order_sc == 'R'), data_bigger.n_dofs_surface + (cfg.order_sc == 'C'), cfg.order_sc, ator_sc);
            }
            data_di.sc.set_view(data_di.n_dofs_surface, data_di.n_dofs_surface, sc_allocd.ld, sc_allocd.order, nullptr, sc_allocd.ator);
            data_di.op_sc_sub_from_allocd.set_matrix_src(&sc_allocd);
            data_di.op_sc_sub_from_allocd.set_matrix_dst(&data_di.sc);
            if(di == di_bigger) {
                data_di.op_sc_sub_from_allocd.set_bounds(0, data_di.sc.nrows, 0, data_di.sc.ncols);
                data_di.sc.prop.uplo = (sc_allocd.order == 'R') ? 'U' : 'L';
            }
            else {
                data_di.op_sc_sub_from_allocd.set_bounds((int)(cfg.order_sc == 'R'), data_di.sc.nrows + (int)(cfg.order_sc == 'R'), (int)(cfg.order_sc == 'C'), data_di.sc.ncols + (int)(cfg.order_sc == 'C'));
                data_di.sc.prop.uplo = (sc_allocd.order == 'R') ? 'L' : 'U';
            }
            data_di.sc.prop.symm = MatrixSymmetry_new::hermitian;
        }
    }
    else {
        sc_allocated.resize(n_domains);
        for(size_t di = 0; di < n_domains; di++) {
            per_domain_data & my = domain_data[di];
            sc_allocated[di].set(my.n_dofs_surface, my.n_dofs_surface, cfg.order_sc, ator_sc);
            my.sc.set_view(my.n_dofs_surface, my.n_dofs_surface, sc_allocated[di].ld, cfg.order_sc, nullptr, sc_allocated[di].ator);
            my.op_sc_sub_from_allocd.set_matrix_src(&sc_allocated[di]);
            my.op_sc_sub_from_allocd.set_matrix_dst(&my.sc);
            my.op_sc_sub_from_allocd.set_bounds(0, my.sc.nrows, 0, my.sc.ncols);
            my.sc.prop.symm = MatrixSymmetry_new::hermitian;
        }
    }

    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_set)
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_data & my = domain_data[di];

        stacktimer::info("subdomain %zu", di);

        if(assemble_apply_where == 'C') {
            my.op_sc = math::operations::schur_csx_dny<T,I>::make(cfg.schur_impl_cpu);
            my.op_sc->set_coefficients(utils::remove_complex_t<T>{1});
            my.op_sc->set_matrix(&my.Kperm, my.n_dofs_surface);
            my.op_sc->set_sc(&my.sc);
            my.op_sc->set_need_solve_A11(false);
        }
        if(assemble_apply_where == 'G') {
            my.op_sc_gpu = gpu::operations::schur_hcsx_ddny<T,I>::make(cfg.schur_impl_gpu);
            my.op_sc_gpu->set_handles(feti.queues[di % feti.queues.size()], feti.handles_sparse[di % feti.handles_sparse.size()], feti.handles_dense[di % feti.handles_dense.size()]);
            my.op_sc_gpu->set_coefficients(utils::remove_complex_t<T>{1});
            my.op_sc_gpu->set_matrix(&my.Kperm, my.n_dofs_surface);
            my.op_sc_gpu->set_sc(&my.sc);
            my.op_sc_gpu->set_need_solve_A11(false);
            my.op_sc_gpu->setup();
        }
    }
    if(!cfg.inner_timers) stacktimer::enable();

    wss_gpu_persistent = 0;
    wss_gpu_internal = 0;
    if(assemble_apply_where == 'G') {
        for(auto & sc : sc_allocated) wss_gpu_persistent += sc.get_memory_impact();
        for(size_t di = 0; di < n_domains; di++) {
            per_domain_data & my = domain_data[di];
            wss_gpu_persistent += my.op_sc_gpu->get_wss_persistent();
            wss_gpu_internal += my.op_sc_gpu->get_wss_internal();
        }
    }
    if(assemble_apply_where == 'G') {
        for(size_t di = 0; di < n_domains; di++) {
            per_domain_data & my = domain_data[di];
            wss_gpu_persistent += my.apply_w_gpu.get_memory_impact();
            wss_gpu_persistent += my.apply_z_gpu.get_memory_impact();
        }
    }

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();

    called_setup = true;
}



template<typename T, typename I>
void DirichletGeneralSchur<T,I>::set(const step::Step &/*step*/)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("DirichletGeneralSchur::set");

    if(!called_setup) eslog::error("setup was not called\n");

    ator_ws_gpu_persistent->set(ws_gpu_persistent, wss_gpu_persistent);

    for(auto & sc : sc_allocated) {
        sc.alloc();
    }

    AllocatorCBMB_new ator_gpu_tmp_cbmba(AllocatorGPU_new::get_singleton(), feti.gpu_tmp_mem, feti.gpu_tmp_size);

    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_set)
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_data & my = domain_data[di];

        stacktimer::info("subdomain %zu", di);

        my.apply_w.alloc();
        my.apply_z.alloc();

        if(assemble_apply_where == 'G') {
            my.apply_w_gpu.alloc();
            my.apply_z_gpu.alloc();
        }

        my.op_sc_sub_from_allocd.perform();

        if(assemble_apply_where == 'C') {
            my.op_sc->preprocess();
        }
        if(assemble_apply_where == 'G') {
            my.op_sc_gpu->set_ws_persistent(ator_ws_gpu_persistent->alloc(my.op_sc_gpu->get_wss_persistent()));
            void * ws_tmp = ator_gpu_tmp_cbmba.alloc(my.op_sc_gpu->get_wss_tmp_preprocess());
            my.op_sc_gpu->preprocess_submit(ws_tmp);
            gpu::mgm::submit_host_function(feti.queues[di % feti.queues.size()], [&,ws_tmp](){
                void * ws_tmp_ = ws_tmp;
                ator_gpu_tmp_cbmba.free(ws_tmp_);
            });
        }
    }
    if(!cfg.inner_timers) stacktimer::enable();

    if(assemble_apply_where == 'G') {
        gpu::mgm::device_wait();
    }

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void DirichletGeneralSchur<T,I>::update(const step::Step &/*step*/)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("DirichletGeneralSchur::update");

    AllocatorCBMB_new ator_gpu_tmp_cbmba(AllocatorGPU_new::get_singleton(), feti.gpu_tmp_mem, feti.gpu_tmp_size);

    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_set)
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_data & my = domain_data[di];

        stacktimer::info("subdomain %zu", di);

        my.op_perm_K.perform_values();

        if(assemble_apply_where == 'C') {
            my.op_sc->perform();
        }
        if(assemble_apply_where == 'G') {
            void * ws_tmp = ator_gpu_tmp_cbmba.alloc(my.op_sc_gpu->get_wss_tmp_perform());
            my.op_sc_gpu->perform_1_submit();
            my.op_sc_gpu->perform_2_submit(ws_tmp);
            gpu::mgm::submit_host_function(feti.queues[di % feti.queues.size()], [&,ws_tmp](){
                void * ws_tmp_ = ws_tmp;
                ator_gpu_tmp_cbmba.free(ws_tmp_);
            });
        }
    }
    if(!cfg.inner_timers) stacktimer::enable();

    if(assemble_apply_where == 'G') {
        gpu::mgm::device_wait();
    }

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void DirichletGeneralSchur<T,I>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    if(cfg.outer_timers) stacktimer::enable();
    stacktimer::push("DirichletGeneralSchur::apply");

    // TODO improve on gpu

    std::fill_n(y.vals, y.size, T{0});
    
    stacktimer::push("apply_all_subdomains");
    if(!cfg.inner_timers) stacktimer::disable();
    #pragma omp parallel for schedule(static,1) if(cfg.parallel_apply)
    for(size_t di = 0; di < n_domains; di++) {
        per_domain_data & my = domain_data[di];
        std::vector<int> & myD2C = feti.D2C[di];

        stacktimer::push("apply subdomain");

        stacktimer::push("C2D,Bt,restrict");
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
        stacktimer::pop();

        // apply_w -> apply SC -> apply_z
        stacktimer::push("mv");
        if(assemble_apply_where == 'C') {
            math::operations::mv_dnx<T>::do_all(&my.sc, &my.apply_w, &my.apply_z, T{1}, T{0});
        }
        if(assemble_apply_where == 'G') {
            gpu::mgm::queue & q = feti.queues[di % feti.queues.size()];
            gpu::dnblas::handle & hd = feti.handles_dense[di % feti.handles_dense.size()];
            gpu::mgm::copy_submit(q, my.apply_w, my.apply_w_gpu);
            if(is_hermitian<T>(my.sc.prop.symm)) {
                gpu::dnblas::hemv<T,I>(hd, my.n_dofs_surface, my.sc.vals, my.sc.ld, my.sc.order, 'N', my.sc.prop.uplo, my.apply_w_gpu.vals, my.apply_z_gpu.vals);
            }
            else {
                gpu::dnblas::gemv<T,I>(hd, my.sc.nrows, my.sc.ncols, my.sc.vals, my.sc.ld, my.sc.order, 'N', my.apply_w_gpu.vals, my.apply_z_gpu.vals);
            }
            gpu::mgm::copy_submit(q, my.apply_z_gpu, my.apply_z);
            gpu::mgm::queue_wait(q);
        }
        stacktimer::pop();

        // apply_z -> extend surface to domain, apply B, distribute to dual vector -> y
        stacktimer::push("extend,B,D2C");
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
        stacktimer::pop();

        stacktimer::pop();
    }
    if(!cfg.inner_timers) stacktimer::enable();
    stacktimer::pop();
    
    y.synchronize();

    stacktimer::pop();
    if(cfg.outer_timers) stacktimer::disable();
}



template<typename T, typename I>
void DirichletGeneralSchur<T,I>::setup_config()
{
    using ecf_config = DirichletGeneralSchurConfig;
    const ecf_config & ecf = feti.configuration.dirichlet_generalschur_config;
    
    switch(ecf.parallel_set) {
        case ecf_config::AUTOBOOL::AUTO:  cfg.parallel_set = true;  break;
        case ecf_config::AUTOBOOL::TRUE:  cfg.parallel_set = true;  break;
        case ecf_config::AUTOBOOL::FALSE: cfg.parallel_set = false; break;
    }

    switch(ecf.parallel_update) {
        case ecf_config::AUTOBOOL::AUTO:  cfg.parallel_update = true;  break;
        case ecf_config::AUTOBOOL::TRUE:  cfg.parallel_update = true;  break;
        case ecf_config::AUTOBOOL::FALSE: cfg.parallel_update = false; break;
    }

    switch(ecf.parallel_apply) {
        case ecf_config::AUTOBOOL::AUTO:  cfg.parallel_apply = true;  break;
        case ecf_config::AUTOBOOL::TRUE:  cfg.parallel_apply = true;  break;
        case ecf_config::AUTOBOOL::FALSE: cfg.parallel_apply = false; break;
    }

    switch(ecf.timers_outer) {
        case ecf_config::AUTOBOOL::AUTO:  cfg.outer_timers = false; break;
        case ecf_config::AUTOBOOL::TRUE:  cfg.outer_timers = true;  break;
        case ecf_config::AUTOBOOL::FALSE: cfg.outer_timers = false; break;
    }

    switch(ecf.timers_inner) {
        case ecf_config::AUTOBOOL::AUTO:  cfg.inner_timers = false; break;
        case ecf_config::AUTOBOOL::TRUE:  cfg.inner_timers = true;  break;
        case ecf_config::AUTOBOOL::FALSE: cfg.inner_timers = false; break;
    }

    switch(ecf.print_config) {
        case ecf_config::AUTOBOOL::AUTO:  cfg.print_config = false; break;
        case ecf_config::AUTOBOOL::TRUE:  cfg.print_config = true;  break;
        case ecf_config::AUTOBOOL::FALSE: cfg.print_config = false; break;
    }

    switch(ecf.order_sc) {
        case ecf_config::MATRIX_ORDER::AUTO:      cfg.order_sc = 'R'; break;
        case ecf_config::MATRIX_ORDER::ROW_MAJOR: cfg.order_sc = 'R'; break;
        case ecf_config::MATRIX_ORDER::COL_MAJOR: cfg.order_sc = 'C'; break;
    }

    switch(ecf.schur_impl_cpu) {
        case ecf_config::SCHUR_IMPL_CPU::AUTO:          cfg.schur_impl_cpu = schur_impl_cpu_t::autoselect;    break;
        case ecf_config::SCHUR_IMPL_CPU::MANUAL_SIMPLE: cfg.schur_impl_cpu = schur_impl_cpu_t::manual_simple; break;
        case ecf_config::SCHUR_IMPL_CPU::TRIANGULAR:    cfg.schur_impl_cpu = schur_impl_cpu_t::triangular;    break;
        case ecf_config::SCHUR_IMPL_CPU::MKLPARDISO:    cfg.schur_impl_cpu = schur_impl_cpu_t::mklpardiso;    break;
        case ecf_config::SCHUR_IMPL_CPU::SPARSE_SOLVER: cfg.schur_impl_cpu = schur_impl_cpu_t::sparse_solver; break;
        case ecf_config::SCHUR_IMPL_CPU::MUMPS:         cfg.schur_impl_cpu = schur_impl_cpu_t::mumps;         break;
        case ecf_config::SCHUR_IMPL_CPU::PASTIX:        cfg.schur_impl_cpu = schur_impl_cpu_t::pastix;        break;
    }

    switch(ecf.schur_impl_gpu) {
        case ecf_config::SCHUR_IMPL_GPU::AUTO:          cfg.schur_impl_gpu = schur_impl_gpu_t::autoselect;    break;
        case ecf_config::SCHUR_IMPL_GPU::MANUAL_SIMPLE: cfg.schur_impl_gpu = schur_impl_gpu_t::manual_simple; break;
        case ecf_config::SCHUR_IMPL_GPU::TRIANGULAR:    cfg.schur_impl_gpu = schur_impl_gpu_t::triangular;    break;
    }
}



#define INSTANTIATE_T_I(T,I) \
template struct DirichletGeneralSchur<T,I>;

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
