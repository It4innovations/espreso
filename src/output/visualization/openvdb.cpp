
#include "openvdb.h"
#include "wrappers/openvdb/w.openvdb.h"

#include "basis/containers/volumepacker.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "basis/containers/allocators.h"
#include "basis/containers/serializededata.h"
#include "wrappers/mpi/communication.h"

#include <vector>
#include <iomanip>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <chrono>

using namespace espreso;

struct SharedVolume {
	static int counter;
	static std::mutex mutex;
	static std::condition_variable cv;

	SharedVolume(const std::string &prefix, int step, int root)
	: root(root), step(step),
	  packer(info::mesh->elements->volumeGrid)
	{
		std::stringstream name;
		name << prefix + "." << std::setw(4) << std::setfill('0') << std::to_string(step) << ".vdb";
		filename = name.str();
	}

	int root, step;
	std::string filename;
	VolumePacker packer;
};

int SharedVolume::counter = 0;
std::mutex SharedVolume::mutex;
std::condition_variable SharedVolume::cv;

OpenVDB::OpenVDB()
: _filename(_path + _directory + _name), _step(0)
{
	nranks.reserve(MPITools::node->within.size);
	nranks.push_back(0);
	int step = MPITools::node->within.size;
	while (step) {
		if (step % 2 == 0) {
			step /= 2;
			size_t size = nranks.size();
			for (size_t i = 0; i < size; ++i) {
				nranks.push_back(nranks[i] + step);
			}
		} else {
			size_t size = nranks.size();
			for (int s = 1; s < step; ++s) {
				for (size_t i = 0; i < size; ++i) {
					if (nranks[i] + s < MPITools::node->within.size) {
						nranks.push_back(nranks[i] + s);
					}
				}
			}
			step = 0;
		}
	}
}

OpenVDB::~OpenVDB()
{
	while (!_postponed.empty()) {
		if (_postponed.front().test()) {
			_postponed.pop();
		}
	}
	std::unique_lock<std::mutex> lk(SharedVolume::mutex);
	SharedVolume::cv.wait(lk, [] { return SharedVolume::counter == 0; });
}

void OpenVDB::updateMesh()
{
	// nothing is done here since we gather mesh to a different process each time step
	std::lock_guard<std::mutex> lk(info::mesh->voxelization.mutex);
	--info::mesh->voxelization.counter;
	info::mesh->voxelization.cv.notify_one();
}

void OpenVDB::updateMonitors()
{
	for (size_t di = 0; di < info::mesh->elements->data.size(); ++di) {
		if (info::mesh->elements->data[di]->name.size()) {
			activeVariables.push_back(info::mesh->elements->data[di]);
		}
	}
}

void store(SharedVolume *volume)
{
	if (volume->root == info::mpi::rank) {
		OpenVDBWrapper wrapper(info::mesh->elements->volumeOrigin, info::mesh->elements->volumeSize, info::mesh->elements->volumeGrid);
		std::vector<OpenVDBWrapper::FloatData*> data(volume->packer.variableSize);
		for (size_t v = 0; v < volume->packer.variableSize; ++v) {
			data[v] = wrapper.addFloat(volume->packer.names[v]);
		}

		volume->packer.unpack([&data] (const _Point<short> &voxel, const size_t &vindex, const float &value) { data[vindex]->insert(voxel, value); });
		for (int r = 1; r < info::mpi::size; ++r) {
			MPI_Recv(volume->packer.packedData, volume->packer.packSize, MPI_BYTE, MPI_ANY_SOURCE, volume->step, MPITools::asynchronous->communicator, MPI_STATUS_IGNORE);
			volume->packer.unpack([&data] (const _Point<short> &voxel, const size_t &vindex, const float &value) { data[vindex]->insert(voxel, value); });
		}
		wrapper.store(volume->filename.c_str());
	}

	delete volume; // we are in detached thread
	{
		std::lock_guard<std::mutex> lk(SharedVolume::mutex);
		--SharedVolume::counter;
	}
	SharedVolume::cv.notify_one();
}

OpenVDB::OpenVDBData::OpenVDBData(SharedVolume *volume)
: volume(volume), root(volume->root), req{}
{
	if (root == info::mpi::rank) {
		{
			std::lock_guard<std::mutex> lk(SharedVolume::mutex);
			++SharedVolume::counter;
		}
		std::thread t(store, volume);
		t.detach();
	} else {
		MPI_Isend(volume->packer.packedData, volume->packer.packSize, MPI_BYTE, root, volume->step, MPITools::asynchronous->communicator, &req);
	}
}

bool OpenVDB::OpenVDBData::test()
{
	if (root == info::mpi::rank) {
		return true;
	} else {
		if (Communication::test(&req)) {
			delete volume;
			return true;
		}
	}
	return false;
}

void OpenVDB::lock()
{
	std::lock_guard<std::mutex> lk(info::mesh->voxelization.mutex);
	++info::mesh->voxelization.counter;
}

void OpenVDB::updateSolution()
{
	if (info::mesh->elements->volumeIndices == nullptr) {
		eslog::warning("SET volumeIndices to store OpenVDB\n");
		return;
	}

	if (_measure) { eslog::startln("OPENVDB: STORING STARTED", "OPENVDB"); }

	if (info::ecf->output.volume_sleep) {
		std::this_thread::sleep_for(std::chrono::milliseconds(info::ecf->output.volume_sleep));
	}

	if (_measure) { eslog::checkpointln("OPENVDB: SLEPT"); }

	while (!_postponed.empty() && _postponed.front().test()) {
		_postponed.pop();
	}

	int node = _step % MPITools::node->across.size;
	int rank = nranks[(_step / MPITools::node->across.size) % MPITools::node->within.size];

	SharedVolume *volume = new SharedVolume(_filename, _step++, node * MPITools::node->within.size + rank);
	volume->packer.analyze(info::mesh->elements->volumeIndices, activeVariables);
	Communication::reduce(&volume->packer.packSize, nullptr, 1, MPITools::getType(volume->packer.packSize).mpitype, MPI_MAX, volume->root, MPITools::asynchronous);
	volume->packer.allocate();
	volume->packer.pack(info::mesh->elements->volumeIndices, activeVariables);

	{ // we can unlock mesh (values are copied by asynchronous output automatically)
		std::lock_guard<std::mutex> lk(info::mesh->voxelization.mutex);
		--info::mesh->voxelization.counter;
		info::mesh->voxelization.cv.notify_one();
	}

	if (_measure) { eslog::checkpointln("OPENVDB: DATA SERIALIZED"); }

	_postponed.push(volume);

	if (_measure) { eslog::endln("OPENVDB: THREAD DETACHED"); }
}
