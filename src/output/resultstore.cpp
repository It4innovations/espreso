
#include "resultstore.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.hpp"
#include "config/ecf/output.h"
#include "basis/utilities/sysutils.h"
#include "basis/utilities/communication.h"

#include "executor/asyncexecutor.h"
#include "executor/directexecutor.h"

#include "monitors/monitoring.h"
#include "visualization/vtklegacy.h"
#include "visualization/ensightgold.h"
#include "visualization/xdmf.h"
#include "visualization/stl.h"
#include "visualization/netgen.h"
#include "visualization/insitu.h"
#include <string>

#ifndef HAVE_ASYNC
struct espreso::Dispatcher {
	void init() {}
	void dispatch() {}
	bool isExecutor() { return false; }
	void setThreadMode() {}
};
#else
#include "async/Dispatcher.h"
struct espreso::Dispatcher: public async::Dispatcher {
void setThreadMode() { async::Config::setMode(async::THREAD); }
};
#endif /* HAVE_ASYNC */


using namespace espreso;

Dispatcher* ResultStore::_dispatcher = NULL;

ResultStoreBase::ResultStoreBase(const Mesh &mesh)
: _mesh(mesh), _path(info::ecf->outpath + "/"), _directory("PREPOSTDATA/"), _name(info::ecf->name), _measure(Mesh::convertDatabase())
{

}

void ResultStoreBase::createOutputDirectory()
{
	utils::createDirectory({ _path, _directory });
}

void ResultStore::createAsynchronizedStore()
{
	if (info::mesh->store != NULL) {
		delete info::mesh->store;
	}

	info::mesh->store = new ResultStore();

	info::mesh->store->_direct = new DirectExecutor(*info::mesh);
	ResultStoreExecutor *executor = info::mesh->store->_direct;
	switch (info::ecf->output.mode) {
	case OutputConfiguration::MODE::SYNC:
		break;
	case OutputConfiguration::MODE::THREAD:
		_dispatcher = new Dispatcher();
		info::mesh->store->_async = new AsyncStore(*info::mesh);
		_dispatcher->setThreadMode();
		executor = info::mesh->store->_async;
		break;
//		case OutputConfiguration::MODE::MPI:
//			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: not implemented OUTPUT::MODE==MPI.";
//			_dispatcher = new async::Dispatcher();
//			_asyncStore->_async = new AsyncStore(mesh, configuration);
//			if (info::mpi::MPIsize == 1) {
//				ESINFO(GLOBAL_ERROR) << "Invalid number of MPI processes. OUTPUT::MODE==MPI required at least two MPI processes.";
//			}
//			if (configuration.output_node_group_size == 0) {
//				ESINFO(GLOBAL_ERROR) << "OUTPUT::OUTPUT_NODE_GROUP_SIZE cannot be 0.";
//			}
//			async::Config::setMode(async::MPI);
//			async::Config::setGroupSize(configuration.output_node_group_size);
//			_dispatcher->setGroupSize(configuration.output_node_group_size);
//			break;
	}

	// TODO: optimize
	if (info::ecf->output.results_store_frequency != OutputConfiguration::STORE_FREQUENCY::NEVER) {
		switch (info::ecf->output.format) {
		case OutputConfiguration::FORMAT::VTK_LEGACY:
			executor->addResultStore(new VTKLegacy(executor->mesh(), !info::ecf->input.convert_database)); break;
		case OutputConfiguration::FORMAT::ENSIGHT:
			executor->addResultStore(new EnSightGold(executor->mesh(), !info::ecf->input.convert_database)); break;
		case OutputConfiguration::FORMAT::XDMF:
			executor->addResultStore(new XDMF(executor->mesh())); break;
		case OutputConfiguration::FORMAT::STL_SURFACE:
			executor->addResultStore(new STL(*info::mesh)); break;
		case OutputConfiguration::FORMAT::NETGEN:
			executor->addResultStore(new Netgen(*info::mesh)); break;
		default:
			eslog::internalFailure("implement the selected output format.\n");
		}

	}
	if (info::ecf->output.monitors_store_frequency != OutputConfiguration::STORE_FREQUENCY::NEVER && info::ecf->output.monitoring.size()) {
		info::mesh->store->_direct->addResultStore(new Monitoring(*info::mesh));
	}
	if (info::ecf->output.catalyst) {
		info::mesh->store->_direct->addResultStore(new InSitu(*info::mesh));
	}

	if (!info::mesh->store->_direct->hasStore()) {
		delete info::mesh->store->_direct;
		info::mesh->store->_direct = NULL;
	}

	if (info::mesh->store->_async != NULL && !info::mesh->store->_async->hasStore()) {
		delete info::mesh->store->_dispatcher;
		delete info::mesh->store->_async;
		info::mesh->store->_dispatcher = NULL;
		info::mesh->store->_async = NULL;
	}

	if (_dispatcher != NULL) {
		_dispatcher->init();
		_dispatcher->dispatch();

//		if (configuration.mode == OutputConfiguration::MODE::MPI) {
//			int computeSize;
//			MPI_Comm_size(_dispatcher->commWorld(), &computeSize);
//			_asyncStore->computeProcesses = computeSize;
//			_asyncStore->storeProcesses = info::mpi::MPIsize - computeSize;
//		}
//
//		info::mpi::comm = _dispatcher->commWorld();
//		MPI_Comm_rank(info::mpi::comm, &info::mpi::rank);
//		MPI_Comm_size(info::mpi::comm, &info::mpi::size);
	}
}

void ResultStore::destroyAsynchronizedStore()
{
	if (info::mesh->store) {
		delete info::mesh->store;
	}
	if (_dispatcher) {
		delete _dispatcher;
	}

	info::mesh->store = NULL;
	_dispatcher = NULL;
}

bool ResultStore::isStoreNode()
{
	return _dispatcher != NULL && _dispatcher->isExecutor();
}

bool ResultStore::isComputeNode()
{
	return !isStoreNode();
}

void ResultStore::suppress()
{
	_allowed = false;
}

void ResultStore::permit()
{
	_allowed = true;
}

bool ResultStore::storeStep()
{
	bool store = false;
	if (_async) store |= _async->storeStep();
	if (_direct) store |= _direct->storeStep();
	return store && _allowed;
}

void ResultStore::updateMesh()
{
	if (_async && _async->hasStore()) _async->updateMesh();
	if (_direct && _direct->hasStore()) _direct->updateMesh();
}

void ResultStore::updateMonitors()
{
	if (_async && _async->hasStore()) _async->updateMonitors();
	if (_direct && _direct->hasStore()) _direct->updateMonitors();
}

void ResultStore::updateSolution()
{
	if (_async && storeStep()) _async->updateSolution();
	if (_direct && storeStep()) _direct->updateSolution();
}

ResultStore::ResultStore(): _allowed(true), _async(NULL), _direct(NULL)
{

}

ResultStore::~ResultStore()
{
	if (_async) delete _async;
	if (_direct) delete _direct;
}


