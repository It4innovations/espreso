
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

}

template <typename T>
FETI<T>::~FETI()
{
    delete iterativeSolver;
    delete projector;
    delete dualOperator;
    delete preconditioner;
}

template <typename T>
bool FETI<T>::set(const step::Step &step)
{
    double start = eslog::time();

    if (K.size() == 1) {
        configuration.method = FETIConfiguration::METHOD::TOTAL_FETI;
    }

    BtL.resize(K.size());
    for (size_t d = 0; d < K.size(); ++d) {
        if (K[d].nrows < x[d].size) { // it is possible when the solver is called with BEM
            math::set(x[d], 0.); // set inner DOFs to zero and resize 'f' to correct size
            f[d].size = K[d].nrows;
            x[d].size = K[d].nrows;
        }
        BtL[d].resize(K[d].nrows);
    }

    int size = K.size();
    Communication::allReduce(&size, &sinfo.domains, 1, MPITools::getType(size).mpitype, MPI_SUM);

    eslog::checkpointln("FETI: SET INFO");

    dualOperator = DualOperator<T>::create(*this, step);
    projector = Projector<T>::create(*this, step);
    preconditioner = Preconditioner<T>::create(*this, step);
    iterativeSolver = IterativeSolver<T>::create(*this, step);

    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
    eslog::info(" = EXTERNAL LINEAR SOLVER %*s = \n", 66, DirectSparseSolver<T>::name());
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");

    double dualop_start = eslog::time();
    dualOperator->set(step);
    double dualop_stop = eslog::time();
    printf("TMP DUAL OPERATOR SET TIME:    %12.6f ms\n", (dualop_stop - dualop_start) * 1000.0);
    projector->set(step);
    preconditioner->set(step);
    iterativeSolver->set(step);

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

    projector->orthonormalizeKernels(step);

    #warning this is temporary just for testing, remove later
    for(int rep = 0; rep < 2; rep++)
    {
        double start = eslog::time();
        dualOperator->update(step);
        double stop = eslog::time();
        printf("TMP DUAL OPERATOR UPDATE TIME: %12.6f ms\n", (stop - start) * 1000.0);
        if(rep == 0)
        {
            // try to clear cpu cache
            size_t size = (size_t{1} << 27) - rand() % 1000;
            double * data = new double[size]; // ~1 GiB
            #pragma omp parallel for
            for(size_t i = 0; i < size; i++) data[i] = std::sin((double)(i+1));
            double result = 0;
            #pragma omp parallel for reduction(+:result)
            for(size_t i = 0; i < size; i++) result += data[i] * std::cos((double)(i+2));
            if(std::abs(result) < 0.001) printf("hello %f\n", result);
            delete[] data;

            // try to clear gpu cache
            dualOperator->clear_gpu_cache();
        }
    }

    // dualOperator->update(step);
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

    std::string error;
    switch (info.error) {
    case IterativeSolverInfo::ERROR::OK:{
        eslog::info("       = ITERATIONS TOTAL                                                    %9d = \n", info.iterations);
        eslog::info("       = NORM(K * U - F + BT * LAMBDA) / NORM(F - BT * LAMBDA)               %.3e = \n", std::sqrt(kkt[0]) / std::sqrt(kkt[1]));
        eslog::info("       = FETI SOLVER TIME                                                   %8.3f s = \n", eslog::time() - start);
        eslog::info("       = ----------------------------------------------------------------------------- = \n");
        return true;
    }
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
    return false;
}

template struct FETI<double>;

}
