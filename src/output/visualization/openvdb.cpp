
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
	: root(root),
	  packer(info::mesh->elements->volumeGrid), nvoxels{}, ndata{}, nvalues{}, split{1}, pack{}, values{}
	{
		std::stringstream name;
		name << prefix + "." << std::setw(4) << std::setfill('0') << std::to_string(step) << ".vdb";
		filename = name.str();
	}

	std::string filename;
	int root;

	VolumePacker packer;
	size_t nvoxels, ndata, nvalues, split;
	char *pack;
	float* values;
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
		for (size_t s = 0; s < volume->split; ++s) {
			int ndata = 0;
			for (size_t di = 0; di < info::mesh->elements->data.size(); di++) { // go through all element values
				if (info::mesh->elements->data[di]->name.size()) {
					OpenVDBWrapper::FloatData *data = wrapper.addFloat(info::mesh->elements->data[di]->name);
					char *pack = volume->pack + s * volume->nvoxels * info::mpi::size;
					float* values = volume->values + ndata * volume->nvalues + s * volume->ndata * volume->nvalues * info::mpi::size;
					for (int r = 0; r < info::mpi::size; ++r) {
						data->insert(volume->packer, pack, values);
						pack += volume->nvoxels;
						values += volume->ndata * volume->nvalues;
					}
					++ndata;
				}
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
	if (Communication::testAll(2 * volume->split, req.data())) {
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

	if (info::ecf->output.volume_sleep) {
		std::this_thread::sleep_for(std::chrono::milliseconds(info::ecf->output.volume_sleep));
	}

	if (_measure) { eslog::startln("OPENVDB: STORING STARTED", "OPENVDB"); }

	while (!_postponed.empty() && _postponed.front().call()) {
		_postponed.pop();
	}

	int node = _step % MPITools::node->across.size;
	int rank = nranks[(_step / MPITools::node->across.size) % MPITools::node->within.size];

	SharedVolume *volume = new SharedVolume(_filename, _step++, node * MPITools::node->within.size + rank);

	size_t datasize[2] = { 0, 0 };
	datasize[0] = volume->packer.analyze(info::mesh->elements->volumeIndices, 0, info::mesh->elements->volumeIndices->structures());
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

	size_t mult = volume->root == info::mpi::rank ? info::mpi::size : 1;
	size_t maxsize = std::max(datasize[0] * info::mpi::size, datasize[1] * info::mpi::size * sizeof(float));
	volume->split = maxsize / (1UL << 30) + 1;
	volume->nvoxels = volume->nvoxels / volume->split + 1;
	volume->nvalues = volume->nvalues / volume->split + 1;

	std::vector<size_t> splitters = { 0UL, info::mesh->elements->volumeIndices->structures() };
	if (volume->split > 1) {
		splitters.pop_back();
		size_t counter = 0, i = 0;
		for (auto e = info::mesh->elements->volumeIndices->cbegin(); e != info::mesh->elements->volumeIndices->cend(); ++e, ++i) {
			if (e->size()) {
				++counter;
			}
			if (counter && counter % volume->nvalues == 0) {
				splitters.push_back(i);
				counter = 0;
			}
		}
		splitters.resize(volume->split + 1, info::mesh->elements->volumeIndices->structures());

		for (size_t s = 0; s < volume->split; ++s) {
			volume->nvoxels = std::max(volume->nvoxels, volume->packer.analyze(info::mesh->elements->volumeIndices, splitters[s], splitters[s + 1]));
		}
		Communication::allReduce(&volume->nvoxels, nullptr, 1, MPITools::getType<size_t>().mpitype, MPI_MAX, MPITools::asynchronous);
	}

	volume->pack = new char[volume->split * volume->nvoxels * mult];
	volume->values = new float[volume->split * volume->ndata * volume->nvalues * mult];

	int offset = volume->root == info::mpi::rank ? info::mpi::rank : 0;
	for (size_t s = 0; s < volume->split; ++s) {
		volume->packer.pack(info::mesh->elements->volumeIndices, splitters[s], splitters[s + 1], volume->pack + offset * volume->nvoxels + s * mult * volume->nvoxels);
	}

	offset *= volume->ndata * volume->nvalues;
	for (size_t s = 0; s < volume->split; ++s) {
		for (size_t di = 0; di < info::mesh->elements->data.size(); di++) { // go through all element values
			if (info::mesh->elements->data[di]->name.size()) {
				size_t i = 0, e = splitters[s];
				for (auto indices = info::mesh->elements->volumeIndices->cbegin() + splitters[s]; indices != info::mesh->elements->volumeIndices->cbegin() + splitters[s + 1]; ++indices, ++e) {
					if (indices->size()) {
						if (info::mesh->elements->data[di]->dimension == 1) {
							volume->values[offset + i] = info::mesh->elements->data[di]->store[e];
						} else {
							float value = 0;
							for (int d = 0; d < info::mesh->elements->data[di]->dimension; ++d) {
								double &v = info::mesh->elements->data[di]->store[e * info::mesh->elements->data[di]->dimension + d];
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
		offset -= volume->ndata * volume->nvalues;
		offset += volume->ndata * volume->nvalues * mult;
	}

	{ // we can unlock mesh (values are copied by asynchronous output automatically)
		std::lock_guard<std::mutex> lk(info::mesh->voxelization.mutex);
		--info::mesh->voxelization.counter;
		info::mesh->voxelization.cv.notify_one();
	}

	if (_measure) { eslog::checkpointln("OPENVDB: DATA SERIALIZED"); }

	_postponed.emplace();
	_postponed.back().volume = volume;
	_postponed.back().req.resize(2 * volume->split);
	for (size_t s = 0; s < volume->split; ++s) {
		Communication::igather(volume->pack + s * volume->nvoxels * mult, nullptr, volume->nvoxels, MPI_BYTE, volume->root, _postponed.back().req[2 * s], MPITools::asynchronous);
		Communication::igather(volume->values + s * volume->ndata * volume->nvalues * mult, nullptr, volume->ndata * volume->nvalues, MPI_FLOAT, volume->root, _postponed.back().req[2 * s + 1], MPITools::asynchronous);
	}

	eslog::info(" == VOXELS GATHERED      ROOT %6d, NODE %4d, VOXELS %10.2f MB, VALUES %10.2f MB == \n", volume->root, node, info::mpi::size * volume->split * volume->nvoxels / 1024. / 1024., info::mpi::size * volume->split * volume->ndata * volume->nvalues * 4 / 1024. / 1024.);
	if (_measure) { eslog::checkpointln("OPENVDB: DATA GATHERED"); }

	if (_postponed.front().call()) {
		_postponed.pop();
	}

	if (_measure) { eslog::endln("OPENVDB: THREAD DETACHED"); }
}
