
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/systeminfo.h"
#include "wrappers/mpi/communication.h"

#include "basis/logging/logger.h"
#include "basis/logging/progresslogger.h"
#include "basis/logging/timelogger.h"
#include "basis/logging/profiler.h"

#include "config/reader/reader.h"
#include "config/configuration.h"
#include "output/output.h"
#include "wrappers/precice/w.precice.h"

using namespace espreso;

int main(int argc, char **argv)
{
    info::system::setSignals();
    info::env::set();
    info::mpi::init(&argc, &argv);
    MPITools::init();

    eslog::init(new Logger<TimeLogger, ProgressTerminalLogger, ProgressFileLogger>);
    eslog::startln("COUPLER: STARTED", "COUPLER");

    ECF::init(&argc, &argv, "dummycoupler");
    MPITools::setSubset(info::ecf->input.third_party_scalability_limit);
    eslog::initFiles();
    info::ecf->output.results_store_frequency = OutputConfiguration::STORE_FREQUENCY::NEVER;
    info::ecf->output.mode = OutputConfiguration::MODE::SYNC;
    eslog::checkpointln("COUPLER: CONFIGURATION READ");

    eslog::printRunInfo(&argc, &argv);
    Mesh::init();
    eslog::checkpointln("COUPLER: RUN INITIALIZED");

    info::mesh->load();
    info::mesh->withSurface = true;
    eslog::checkpointln("COUPLER: MESH LOADED");

    info::mesh->preprocess();
    eslog::checkpointln("COUPLER: MESH PREPROCESSED");
    info::mesh->printMeshStatistics();

    Precice::dummy();
    eslog::endln("COUPLER: PRECICE STARTED");

    Mesh::finish();
    eslog::finish();
    MPITools::finish();
    ECF::finish();

    profiler::syncend("mesio");
    profiler::print(); // need to be printed before MPI_Finalize

    info::mpi::finish();
}
