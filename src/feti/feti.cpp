
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

    int offset[2] = { 0, 0 }, ooffset[2];
    int size[4] = { 0, 0, 0, 0 };
    size[0] = K.size();

    bool fetiProjection =
            configuration.method == FETIConfiguration::METHOD::TOTAL_FETI ||
            configuration.projector == FETIConfiguration::PROJECTOR::ORTHOGONAL_FULL ||
            configuration.projector == FETIConfiguration::PROJECTOR::ORTHOGONAL_FULL_WITH_FACTORS;

    if (fetiProjection) {
        for (size_t d = 0; d < K.size(); ++d) {
            sinfo.R1size = offset[0] = size[2] += R1[d].nrows;
            sinfo.R2size = offset[1] = size[3] += R2[d].nrows;
        }
    } else {
        for (size_t d = 0; d < K.size(); ++d) {
            sinfo.R1size = offset[0] = size[2] = std::max(size[2], R1[d].nrows);
            sinfo.R2size = offset[1] = size[3] = std::max(size[3], R2[d].nrows);
        }
    }

    for (size_t d = 0; d < K.size(); ++d) {
        if (K[d].nrows < x[d].size) { // it is possible when the solver is called with BEM
            math::set(x[d], 0.); // set inner DOFs to zero and resize 'f' to correct size
            f[d].size = K[d].nrows;
            x[d].size = K[d].nrows;
        }
    }
    sinfo.eq_offset = lambdas.eq_size;
    sinfo.nc_offset = lambdas.nc_size;
    sinfo.eq_total = Communication::exscan(sinfo.eq_offset);
    sinfo.nc_total = Communication::exscan(sinfo.nc_offset);
    sinfo.dual_total = sinfo.eq_total + sinfo.nc_total;

    Communication::exscan(offset, ooffset, 2, MPITools::getType<int>().mpitype, MPI_SUM);
    Communication::allReduce(size, NULL, 5, MPITools::getType<int>().mpitype, MPI_SUM);
    sinfo.domains = size[0];
    sinfo.R1totalSize = size[2];
    sinfo.R2totalSize = size[3];
    sinfo.R1offset = info::mpi::rank ? ooffset[0] : 0;
    sinfo.R2offset = info::mpi::rank ? ooffset[1] : 0;

    Dual_Map::set(*this);
    Vector_Dual<T>::initBuffers();
    Vector_Kernel<T>::set(sinfo.R1offset, sinfo.R1size, sinfo.R1totalSize);

    eslog::checkpointln("FETI: SET INFO");

    iterativeSolver = IterativeSolver<T>::set(*this, step);
    projector = Projector<T>::set(*this, step);
    dualOperator = DualOperator<T>::set(*this, step);
    preconditioner = Preconditioner<T>::set(*this, step);

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
    projector->update(step);
    dualOperator->update(step);
    preconditioner->update(step);
    eslog::info("       = FETI SOLVER UPDATED                                                %8.3f s = \n", eslog::time() - start);
    return true;
}

template <typename T>
bool FETI<T>::solve(const step::Step &step)
{
    double start = eslog::time();
    IterativeSolverInfo info;
    iterativeSolver->solve(step, info);

    std::string error;
    switch (info.error) {
    case IterativeSolverInfo::ERROR::OK:{
        eslog::info("       = ITERATIONS TOTAL                                                    %9d = \n", info.iterations);
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
