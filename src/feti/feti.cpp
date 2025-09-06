
#include "feti.h"
#include "dualoperator/dualoperator.h"
#include "iterativesolver/pcpg.h"
#include "projector/projector.h"
#include "preconditioner/preconditioner.h"
#include "math/math.h"
#include "esinfo/eslog.hpp"
#include "esinfo/systeminfo.h"
#include "wrappers/mpi/communication.h"

namespace espreso {

template <typename T>
FETI<T>::FETI(FETIConfiguration &configuration)
: configuration(configuration), decomposition(nullptr)
{
    use_gpu = gpu::mgm::is_linked() && gpu::mgm::is_available();
    if(use_gpu) {
        device = gpu::mgm::get_device_by_mpi(info::mpi::rank, info::mpi::size);

        gpu::mgm::init_gpu(device);
        gpu::mgm::set_device(device);

        size_t n_queues = omp_get_max_threads();

        queues.resize(n_queues);
        handles_dense.resize(n_queues);
        handles_sparse.resize(n_queues);

        gpu::mgm::queue_create(main_q);
        for(gpu::mgm::queue & q : queues) gpu::mgm::queue_create(q);
        for(size_t i = 0; i < n_queues; i++) gpu::dnblas::handle_create(handles_dense[i], queues[i]);
        for(size_t i = 0; i < n_queues; i++) gpu::spblas::handle_create(handles_sparse[i], queues[i]);

        gpu::dnblas::init_library(main_q);
        gpu::mgm::queue_wait(main_q);
    }
}

template <typename T>
FETI<T>::~FETI()
{
    delete iterativeSolver;
    delete projector;
    delete dualOperator;
    delete preconditioner;

    if(use_gpu) {
        gpu::mgm::memfree_device(gpu_mem_allocd);
    }
}

template <typename T>
bool FETI<T>::set(const step::Step &step)
{
    double start = eslog::time();

    if (K.size() == 1) {
        configuration.method = FETIConfiguration::METHOD::TOTAL_FETI;
    }

    for (size_t d = 0; d < K.size(); ++d) {
        if (K[d].nrows < x[d].size) { // it is possible when the solver is called with BEM
            math::set(x[d], 0.); // set inner DOFs to zero and resize 'f' to correct size
            f[d].size = K[d].nrows;
            x[d].size = K[d].nrows;
        }
    }

    int size = K.size();
    Communication::allReduce(&size, &sinfo.domains, 1, MPITools::getType(size).mpitype, MPI_SUM);

    eslog::checkpointln("FETI: SET INFO");

    dualOperator = DualOperator<T>::create(*this, step);
    projector = Projector<T>::create(*this, step);
    preconditioner = Preconditioner<T>::create(*this, step);
    iterativeSolver = IterativeSolver<T>::create(*this, step);

    dualOperator->setup();

    if(use_gpu) {
        size_t total_wss_internal = 0;
        total_wss_internal += dualOperator->get_wss_gpu_internal();
        
        size_t free_mem = gpu::mgm::get_device_memory_free();
        size_t mem_capacity = gpu::mgm::get_device_memory_capacity();
        size_t reserve = (mem_capacity * 5) / 100;
        size_t to_alloc = utils::round_down(free_mem - reserve - total_wss_internal, gpu::mgm::get_natural_pitch_align());
        gpu_mem_allocd = gpu::mgm::memalloc_device(to_alloc);
        ator_gpu_arena = std::make_unique<AllocatorArena_new>(AllocatorGPU_new::get_singleton());
        ator_gpu_arena->set(gpu_mem_allocd, to_alloc);

        dualOperator->set_ws_gpu_persistent(ator_gpu_arena->alloc(dualOperator->get_wss_gpu_persistent()));

        gpu_tmp_size = ator_gpu_arena->get_remaining_capacity();
        gpu_tmp_mem = ator_gpu_arena->alloc(gpu_tmp_size);
    }

    dualOperator->set(step);
    projector->set(step);
    preconditioner->set(step);
    iterativeSolver->set(step);

    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
    eslog::info(" = EXTERNAL LINEAR SOLVER %*s = \n", 66, DirectSparseSolver<T>::name());
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
    iterativeSolver->info();
    projector->info();
    dualOperator->info();
    preconditioner->info();

    eslog::info(" = FETI SOLVER SET                                                                %8.3f s = \n", eslog::time() - start);
    if (MPITools::node->rank == 0) {
        info::system::memory::solver = info::system::memoryAvail();
    }
    eslog::info(" = FETI SOLVER MEMORY FOOTPRINT [GB] %55.2f = \n", (info::system::memory::physics - info::system::memory::solver) / 1024. / 1024.);

    return true;
}

template <typename T>
bool FETI<T>::update(const step::Step &step)
{
    double start = eslog::time();
    if (updated.B) {
        Dual_Map::set(*this);
        Vector_Dual<T>::initBuffers();
    }

    if (configuration.check_input_matrices) {
        check();
    }

    projector->orthonormalizeKernels(step);
    dualOperator->update(step);
    projector->update(step);
    preconditioner->update(step);
    iterativeSolver->update(step);
    eslog::info("       = FETI SOLVER UPDATED                                                %8.3f s = \n", eslog::time() - start);
    return true;
}

template <typename T>
bool FETI<T>::solve(const step::Step &step)
{
    double start = eslog::time();
    IterativeSolverInfo info;
    iterativeSolver->solve(step, info);

    T kkt[2] = { T{0}, T{0} };
    #pragma omp parallel for
    for (size_t di = 0; di < K.size(); ++di) {
        SpBLAS<Matrix_CSR, T> spblas(K[di]);
        Vector_Dense<T> Ku; Ku.resize(K[di].nrows);
        Vector_Dense<T> fBtL; fBtL.resize(K[di].nrows);
        spblas.apply(Ku, T{1}, T{0}, x[di]);
        math::copy(fBtL, BtL[di]);
        math::add(fBtL, T{-1}, f[di]);
        #pragma omp atomic
        kkt[1] += math::dot(fBtL, fBtL);
        math::add(Ku, T{1}, fBtL);
        #pragma omp atomic
        kkt[0] += math::dot(Ku, Ku);
    }
    Communication::allReduce(kkt, nullptr, 2, MPITools::getType<T>().mpitype, MPI_SUM);

    eslog::info("       = ITERATIONS TOTAL                                                    %9d = \n", info.iterations);
    eslog::info("       = NORM(K * U - F + BT * LAMBDA) / NORM(F - BT * LAMBDA)               %.3e = \n", std::sqrt(kkt[0]) / std::sqrt(kkt[1]));
    eslog::info("       = FETI SOLVER TIME                                                   %8.3f s = \n", eslog::time() - start);

    std::string error;
    switch (info.error) {
    case IterativeSolverInfo::ERROR::OK:                 eslog::info("       = ----------------------------------------------------------------------------- = \n"); return true;
    case IterativeSolverInfo::ERROR::STAGNATION:             error = "       = FETI SOLVER ERROR                                   NONDECREASING CONVERGENCE = \n"; break;
    case IterativeSolverInfo::ERROR::MAX_ITERATIONS_REACHED: error = "       = FETI SOLVER ERROR                                  MAXIMUM ITERATIONS REACHED = \n"; break;
    case IterativeSolverInfo::ERROR::INVALID_DATA:           error = "       = FETI SOLVER ERROR                                          INVALID INPUT DATA = \n"; break;
    case IterativeSolverInfo::ERROR::CONVERGENCE_ERROR:      error = "       = FETI SOLVER ERROR              SOLVER DOES NOT CONVERGE TO THE REQUESTED NORM = \n"; break;
    }
    if (configuration.exit_on_nonconvergence) {
        eslog::globalerror(error.c_str());
    } else {
        eslog::warning(error.c_str());
    }
    eslog::info("       = ----------------------------------------------------------------------------- = \n");
    return false;
}

template <typename T>
void FETI<T>::check()
{
    T eps = std::numeric_limits<T>::epsilon();

    int invalid = false;
    std::vector<int> zero_s(K.size()), zero_s_reg(K.size());
    std::vector<T> maxK(K.size()), normK(K.size());
    std::vector<std::vector<T> > normKN(K.size());
    std::vector<Vector_Dense<T> > eig(K.size());
    std::vector<Vector_Dense<T> > eig_reg(K.size());
    std::vector<Vector_Dense<T> > s(K.size());
    std::vector<Vector_Dense<T> > s_reg(K.size());
    std::vector<Matrix_Dense<T> > U(K.size()), V(K.size());
    std::vector<Matrix_Dense<T> > U_reg(K.size()), V_reg(K.size());

    #pragma omp parallel for
    for (size_t di = 0; di < K.size(); ++di) {
        T tol = std::max(K[di].nrows, K[di].ncols) * eps;
        maxK[di] = *std::max_element(K[di].vals, K[di].vals + K[di].nnz);

        // K * N
        SpBLAS<Matrix_CSR, T> spblas(K[di]);
        for (int r = 0; r < K[di].nrows; ++r) {
            normK[di] += K[di].vals[K[di].rows[r] - Indexing::CSR] * K[di].vals[K[di].rows[r] - Indexing::CSR];
        }
        normK[di] = std::sqrt(normK[di]);
        Vector_Dense<T> R, KR; KR.resize(R1[di].ncols);
        for (int r = 0; r < R1[di].nrows; ++r) {
            R.size = R1[di].ncols; R.vals = R1[di].vals + R1[di].ncols * r;
            spblas.apply(KR, T{1}, T{0}, R);
            normKN[di].push_back(math::norm(KR));
        }

        // eigenvalues
        Matrix_Dense<T> m; m.resize(K[di].nrows, K[di].ncols); math::set(m, T{0});
        Matrix_Dense<T> m_reg; m_reg.resize(K[di].nrows, K[di].ncols); math::set(m_reg, T{0});
        eig[di].resize(K[di].nrows);
        eig_reg[di].resize(K[di].nrows);

        for (esint r = 0; r < K[di].nrows; ++r) {
            for (esint c = K[di].rows[r]; c < K[di].rows[r + 1]; ++c) {
                m.vals[r * K[di].ncols + K[di].cols[c - Indexing::CSR] - Indexing::CSR] = K[di].vals[c - Indexing::CSR];
                m_reg.vals[r * K[di].ncols + K[di].cols[c - Indexing::CSR] - Indexing::CSR] = K[di].vals[c - Indexing::CSR];
            }
        }
        for (esint r = 0; r < RegMat[di].nrows; ++r) {
            for (esint c = RegMat[di].rows[r]; c < RegMat[di].rows[r + 1]; ++c) {
                m_reg.vals[r * RegMat[di].ncols + RegMat[di].cols[c - Indexing::CSR] - Indexing::CSR] += RegMat[di].vals[c - Indexing::CSR];
            }
        }
        if (K[di].shape == Matrix_Shape::UPPER) {
            for (esint r = 0; r < m.nrows; ++r) {
                for (esint c = r; c < m.ncols; ++c) {
                    m.vals[c * m.ncols + r] = m.vals[r * m.ncols + c];
                    m_reg.vals[c * m.ncols + r] = m_reg.vals[r * m.ncols + c];
                }
            }
        }
        if (K[di].shape == Matrix_Shape::LOWER) {
            for (esint r = 0; r < m.nrows; ++r) {
                for (esint c = r; c < m.ncols; ++c) {
                    m.vals[r * m.ncols + c] = m.vals[c * m.ncols + r];
                    m_reg.vals[r * m.ncols + c] = m_reg.vals[c * m.ncols + r];
                }
            }
        }

        Matrix_Dense<T> n(m), n_reg(m_reg);
        math::lapack::get_eig_sym(m, eig[di]);
        math::lapack::get_eig_sym(m_reg, eig_reg[di]);
        math::lapack::get_svd(n, s[di], U[di], V[di]);
        math::lapack::get_svd(n_reg, s_reg[di], U_reg[di], V_reg[di]);

        T s_max = *std::max_element(s[di].vals, s[di].vals + s[di].size);
        T s_reg_max = *std::max_element(s_reg[di].vals, s_reg[di].vals + s_reg[di].size);
        for (int i = 0; i < s[di].size; ++i) {
            if (s[di].vals[i] <= tol * s_max) {
                zero_s[di] += 1;
            }
        }
        for (int i = 0; i < s_reg[di].size; ++i) {
            if (s_reg[di].vals[i] <= tol * s_reg_max) {
                zero_s_reg[di] += 1;
            }
        }
    }

    Communication::serialize([&] () {
        printf(" !! CHECK \n");
        printf(" !! ================\n");
        for (size_t di = 0; di < K.size(); ++di) {

            printf(" !! DOMAIN          : %d \n", (int)di);
            printf(" !! norm(K)         : %+e\n", normK[di]);
            printf(" !! norm(KN)        :");
            for (size_t v = 0; v < normKN[di].size(); ++v) {
                printf(" %+.2e", normKN[di][v]);
            }
            printf("\n");
            printf(" !! norm(KN)/norm(K):");
            for (size_t v = 0; v < normKN[di].size(); ++v) {
                printf(" %+.2e", normKN[di][v] / normK[di]);
            }
            printf("\n");
            printf(" !! eig(K)          :");
            for (int v = 0; v < 8 && v < eig[di].size; ++v) {
                printf(" %+.2e", eig[di].vals[v]);
            }
            if (8 < eig[di].size) {
                printf(" ... %+.2e", eig[di].vals[eig[di].size - 1]);
            }
            printf("\n");
            printf(" !! eig(K+RegMat)   :");
            for (int v = 0; v < 8 && v < eig_reg[di].size; ++v) {
                printf(" %+.2e", eig_reg[di].vals[v]);
            }
            if (8 < eig_reg[di].size) {
                printf(" ... %+.2e", eig_reg[di].vals[eig_reg[di].size - 1]);
            }
            printf("\n");
            printf(" !! SVD(K)        S :");
            printf(" %+.2e ...", s[di].vals[0]);
            for (int v = std::max(1, s[di].size - 8); v < s[di].size; ++v) {
                printf(" %+.2e", s[di].vals[v]);
            }
            printf("\n");
            printf(" !! SVD(K+RegMat) S :");
            printf(" %+.2e ...", s_reg[di].vals[0]);
            for (int v = std::max(1, s_reg[di].size - 8); v < s_reg[di].size; ++v) {
                printf(" %+.2e", s_reg[di].vals[v]);
            }
            printf("\n");
            printf(" !! S[0] / S[-1]    : %e\n", s_reg[di].vals[0] / s_reg[di].vals[s_reg[di].size - 1]);
            printf(" !! R.cols          : %d\n", R1[di].nrows);
            printf(" !! zero in S       : %d\n", zero_s[di]);
            printf(" !! zero in S_REG   : %d\n", zero_s_reg[di]);
            if (zero_s[di] != R1[di].nrows) {
                invalid = 1;
                printf(" !! INVALID NUMBER OF KERNELS\n");
            }
            if (zero_s_reg[di] && MoorePenroseInv[di].nrows == 0) {
                invalid = 1;
                printf(" !! INVALID REGULARIZATION\n");
            }
            printf(" !! ----------------\n");
        }
    });
    Communication::allReduce(&invalid, nullptr, 1, MPITools::getType(invalid).mpitype, MPI_SUM);
    if (invalid) {
        eslog::globalerror("invalid input matrices for FETI solver.\n");
    }
}

template struct FETI<double>;

}
