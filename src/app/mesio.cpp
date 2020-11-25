
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/systeminfo.h"
#include "basis/utilities/communication.h"

#include "basis/logging/logger.h"
#include "basis/logging/progresslogger.h"
#include "basis/logging/timelogger.h"
#include "basis/logging/profiler.h"

#include "config/reader/reader.h"
#include "config/configuration.h"
#include "output/resultstore.h"

using namespace espreso;

int main(int argc, char **argv)
{
	profiler::syncstart("mesio");

	info::system::setSignals();
	info::env::set();
	info::mpi::init(&argc, &argv);
	MPITools::init();

	eslog::init(new Logger<TimeLogger, ProgressTerminalLogger, ProgressFileLogger>);
	eslog::startln("MESIO: STARTED", "MESIO");

	ECF::init(&argc, &argv, "mesio");
	MPITools::setSubset(info::ecf->input.third_party_scalability_limit);
	eslog::initFiles();
	info::ecf->output.mode = OutputConfiguration::MODE::SYNC;
	eslog::checkpointln("MESIO: CONFIGURATION READ");

	eslog::printRunInfo(&argc, &argv);
	Mesh::init();
	ResultStore::createAsynchronizedStore();
	eslog::checkpointln("MESIO: RUN INITIALIZED");

	Mesh::load();
	eslog::checkpointln("MESIO: MESH LOADED");

	info::mesh->preprocess();
	eslog::checkpointln("MESIO: MESH PREPROCESSED");
	info::mesh->printMeshStatistics();
	info::mesh->printDecompositionStatistics();

	info::mesh->store->updateMesh();
	eslog::endln("MESIO: MESH STORED");

	eslog::finish();
	ResultStore::destroyAsynchronizedStore();
	Mesh::finish();
	MPITools::finish();
	ECF::finish();

	profiler::syncend("mesio");
	profiler::print(); // need to be printed before MPI_Finalize

	info::mpi::finish();

	return 0;
}
