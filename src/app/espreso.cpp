
#include "analysis/analysis.h"
#include "esinfo/eslog.hpp"
#include "esinfo/stepinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/systeminfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "wrappers/mpi/communication.h"

#include "basis/logging/logger.h"
#include "basis/logging/progresslogger.h"
#include "basis/logging/timelogger.h"
#include "basis/logging/papicounters.h"
#include "basis/logging/profiler.h"
#include "basis/utilities/stacktimer.h"

#include "config/reader/reader.h"
#include "config/configuration.h"
#include "mesh/mesh.h"
#include "output/output.h"


using namespace espreso;

int main(int argc, char **argv)
{
    profiler::syncstart("espreso");

    profiler::syncstart("initialization");
    info::system::setSignals();
    info::env::set();
    profiler::synccheckpoint("set_signals_and_env");
    info::mpi::init(&argc, &argv);
    profiler::synccheckpoint("mpi_init");
    MPITools::init();
    profiler::synccheckpoint("mpi_init_tools");

    eslog::init(new Logger<TimeLogger, ProgressTerminalLogger, ProgressFileLogger, PAPICounters>);
    profiler::synccheckpoint("init_loggers");
    eslog::startln("ESPRESO: STARTED", "ESPRESO");

    stacktimer::init();

    ECF::init(&argc, &argv, "espreso");
    profiler::synccheckpoint("init_configuration");
    eslog::checkpointln("ESPRESO: CONFIGURATION READ");
    eslog::startln("CONFIGURATION STARTED", "CONFIGURATION");

    bool divided = info::mpi::divide(info::ecf->input.decomposition.mesh_duplication);
    MPITools::setSubset(info::ecf->input.third_party_scalability_limit);
    eslog::initFiles();
    profiler::synccheckpoint("divide_mpi");
    eslog::printRunInfo(&argc, &argv);
    profiler::synccheckpoint("init_run_info");
    if (!divided) {
        eslog::globalerror("Cannot set MESH DUPLICATION: the number of MPI processes is not divisible by %d\n", info::ecf->input.decomposition.mesh_duplication);
    }
    eslog::checkpointln("CONFIGURATION: RUN INFO INITIALIZED");

    Mesh::init();
    profiler::synccheckpoint("init_mesh");
    eslog::endln("CONFIGURATION: MESH INITIALIZED");
    eslog::checkpointln("ESPRESO: RUN INITIALIZED");
    profiler::syncend("initialization");

    if (info::mpi::irank == 0) {
        profiler::syncstart("mesh_preprocessing");
        info::mesh->load();
        eslog::checkpointln("ESPRESO: MESH LOADED");
        info::mesh->preprocess();
        eslog::checkpointln("ESPRESO: MESH PREPROCESSED");
    }
    if (info::mpi::isize > 1) {
        info::mesh->duplicate();
        eslog::checkpoint("ESPRESO: MESH DUPLICATED");
        eslog::param("COPY", info::mpi::irank);
        eslog::ln();
    }
    info::mesh->printMeshStatistics();
    info::mesh->printDecompositionStatistics();
    profiler::syncend("mesh_preprocessing");

    profiler::syncstart("mesh_output");
    info::mesh->output->updateMesh();
    if (info::ecf->output.mode == OutputConfiguration::MODE::SYNC) {
        eslog::checkpointln("ESPRESO: MESH STORED");
    }
    profiler::syncend("mesh_output");

    if (Mesh::convertDatabase()) {
        if (info::ecf->output.store_decomposition > 1) {
            info::mesh->output->updateMonitors(step::Step());
            info::mesh->output->updateSolution(step::Step(), step::Time());
        }
        eslog::endln("ESPRESO: DATABASE CONVERTED");
    } else {
        profiler::syncstart("physical_solver");
        Analysis looper;
        looper.run();
        profiler::syncend("physical_solver");
        eslog::endln("ESPRESO: SIMULATION FINISHED");
    }

    Mesh::finish();
    eslog::finish();
    profiler::syncend("espreso");
    profiler::print(); // need to be printed before MPI_Finalize

    ECF::finish();
    MPITools::finish();
    info::mpi::finish();
    stacktimer::finish();
    return 0;
}
