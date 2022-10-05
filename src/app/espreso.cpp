
#include "esinfo/eslog.hpp"
#include "esinfo/stepinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/systeminfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "wrappers/mpi/communication.h"
#include "wrappers/cuda/w.cuda.h"

#include "basis/logging/logger.h"
#include "basis/logging/progresslogger.h"
#include "basis/logging/timelogger.h"
#include "basis/logging/oldtimelogger.h"
#include "basis/logging/profiler.h"

#include "config/reader/reader.h"
#include "config/configuration.h"
#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "output/output.h"
#include "physics/physicalsolver.h"

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

	eslog::init(new Logger<TimeLogger, ProgressTerminalLogger, ProgressFileLogger, OldTimeLogger>);
	profiler::synccheckpoint("init_loggers");
	eslog::startln("ESPRESO: STARTED", "ESPRESO");

	ECF::init(&argc, &argv, "espreso");
	profiler::synccheckpoint("init_configuration");
	eslog::checkpointln("ESPRESO: CONFIGURATION READ");
	eslog::startln("CONFIGURATION STARTED", "CONFIGURATION");


	bool divided = info::mpi::divide(info::ecf->input.decomposition.mesh_duplication);
	MPITools::setSubset(info::ecf->input.third_party_scalability_limit);
	eslog::initFiles();
	profiler::synccheckpoint("divide_mpi");
	eslog::printRunInfo(&argc, &argv);
	cuda::fillDeviceInfo(); // it also prints the info
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
		Mesh::load();
	}

	if (!info::ecf->input.convert_database) {
		profiler::syncstart("physical_solver");
		PhysicalSolver::run();
		profiler::syncend("physical_solver");
	}

	Mesh::finish();
	eslog::endln("ESPRESO: ESPRESO FINISHED");
	eslog::finish();
	profiler::syncend("espreso");
	profiler::print(); // need to be printed before MPI_Finalize

	ECF::finish();
	MPITools::finish();
	info::mpi::finish();
	return 0;
}
