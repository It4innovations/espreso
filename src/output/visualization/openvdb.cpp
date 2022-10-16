
#include "openvdb.h"
#include "wrappers/openvdb/w.openvdb.h"

#include "basis/containers/volumepacker.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"
#include "basis/containers/allocators.h"
#include "basis/containers/serializededata.h"
#include "wrappers/mpi/communication.h"

#include <vector>
#include <iomanip>
#include <thread>
#include <mutex>
#include <condition_variable>

using namespace espreso;

struct SharedVolume {
	static int counter;
	static std::mutex mutex;
	static std::condition_variable cv;

	SharedVolume(const std::string &prefix, int step, int root)
	: step(step), root(root),
	  packer(info::mesh->elements->volumeGrid), nvoxels{}, ndata{}, nvalues{}, pack{}, values{}
	{
		std::stringstream name;
		name << prefix + "." << std::setw(4) << std::setfill('0') << std::to_string(step) << ".vdb";
		filename = name.str();
	}

	std::string filename;
	int step, root;

	VolumePacker packer;
	size_t nvoxels, ndata, nvalues;
	char *pack;
	float* values;
};

int SharedVolume::counter = 0;
std::mutex SharedVolume::mutex;
std::condition_variable SharedVolume::cv;

OpenVDB::OpenVDB()
: _filename(_path + _directory + _name), _step(0)
{

}

OpenVDB::~OpenVDB()
{
	while (!_postponed.empty()) {
		if (_postponed.front().call()) {
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

void store(SharedVolume *volume)
{
	if (volume->root == info::mpi::rank) {
		OpenVDBWrapper wrapper(info::mesh->elements->volumeOrigin, info::mesh->elements->volumeSize, info::mesh->elements->volumeGrid);
		int ndata = 0;
		for (size_t di = 0; di < info::mesh->elements->data.size(); di++) { // go through all element values
			if (info::mesh->elements->data[di]->name.size()) {
				OpenVDBWrapper::FloatData *data = wrapper.addFloat(info::mesh->elements->data[di]->name);
				char *pack = volume->pack;
				float* values = volume->values + ndata * volume->nvalues;
				for (int r = 0; r < info::mpi::size; ++r) {
					data->insert(volume->packer, pack, values);
					pack += volume->nvoxels;
					values += volume->ndata * volume->nvalues;
				}
				++ndata;
			}
		}

		wrapper.store(volume->filename.c_str());
	}

	delete[] volume->pack;
	delete[] volume->values;
	delete volume; // we are in detached thread
	{
		std::lock_guard<std::mutex> lk(SharedVolume::mutex);
		--SharedVolume::counter;
	}
	SharedVolume::cv.notify_one();
}

bool OpenVDB::OpenVDBData::call()
{
	if (Communication::testAll(2, req)) {
		{
			std::lock_guard<std::mutex> lk(SharedVolume::mutex);
			++SharedVolume::counter;
		}
		std::thread t(store, volume);
		t.detach();
		return true;
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

	while (!_postponed.empty() && _postponed.front().call()) {
		_postponed.pop();
	}

	int node = _step % MPITools::node->across.size;
	int rank = (_step / MPITools::node->across.size) % MPITools::node->within.size;

	SharedVolume *volume = new SharedVolume(_filename, _step++, node * MPITools::node->within.size + rank);

	size_t datasize[2] = { 0, 0 };
	datasize[0] = volume->packer.analyze(info::mesh->elements->volumeIndices);
	for (auto e = info::mesh->elements->volumeIndices->cbegin(); e != info::mesh->elements->volumeIndices->cend(); ++e) {
		if (e->size()) ++datasize[1];
	}
	Communication::allReduce(datasize, nullptr, 2, MPITools::getType<size_t>().mpitype, MPI_MAX, MPITools::asynchronous);
	volume->nvoxels = datasize[0];
	volume->nvalues = datasize[1];

	for (size_t di = 0; di < info::mesh->elements->data.size(); di++) { // go through all element values
		if (info::mesh->elements->data[di]->name.size()) {
			++volume->ndata;
		}
	}

	int mult = volume->root == info::mpi::rank ? info::mpi::size : 1;
	volume->pack = new char[volume->nvoxels * mult];
	volume->values = new float[volume->ndata * volume->nvalues * mult];

	int offset = volume->root == info::mpi::rank ? info::mpi::rank : 0;
	volume->packer.pack(info::mesh->elements->volumeIndices, volume->pack + offset * volume->nvoxels);

	{ // we can unlock mesh (values are copied by asynchronous output automatically)
		std::lock_guard<std::mutex> lk(info::mesh->voxelization.mutex);
		--info::mesh->voxelization.counter;
		info::mesh->voxelization.cv.notify_one();
	}

	offset *= volume->ndata * volume->nvalues;
	for (size_t di = 0; di < info::mesh->elements->data.size(); di++) { // go through all element values
		if (info::mesh->elements->data[di]->name.size()) {
			size_t i = 0, e = 0;
			for (auto indices = info::mesh->elements->volumeIndices->cbegin(); indices != info::mesh->elements->volumeIndices->cend(); ++indices, ++e) {
				if (indices->size()) {
					if (info::mesh->elements->data[di]->dimension == 1) {
						volume->values[offset + i] = info::mesh->elements->data[di]->data[e];
					} else {
						float value = 0;
						for (int d = 0; d < info::mesh->elements->data[di]->dimension; ++d) {
							double &v = info::mesh->elements->data[di]->data[e * info::mesh->elements->data[di]->dimension + d];
							value += v * v;
						}
						volume->values[offset + i] = std::sqrt(value);
					}
					++i;
				}
			}
			offset += volume->nvalues;
		}
	}

	if (_measure) { eslog::checkpointln("OPENVDB: DATA SERIALIZED"); }

	_postponed.emplace();
	_postponed.back().volume = volume;
	Communication::igather(volume->pack, nullptr, volume->nvoxels, MPI_BYTE, volume->root, _postponed.back().req[0], MPITools::asynchronous);
	Communication::igather(volume->values, nullptr, volume->ndata * volume->nvalues, MPI_FLOAT, volume->root, _postponed.back().req[1], MPITools::asynchronous);

	if (_measure) { eslog::checkpointln("OPENVDB: DATA GATHERED"); }

	if (_postponed.front().call()) {
		_postponed.pop();
	}

	if (_measure) { eslog::endln("OPENVDB: THREAD DETACHED"); }
}
