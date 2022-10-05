
#include "openvdb.h"
#include "wrappers/openvdb/w.openvdb.h"

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
	static esint size;
	static std::vector<esint> sizes;

	SharedVolume(const std::string &prefix, int step, int root)
	: step(step), root(root),
	  nvoxels{}, nvalues{}, displacement{}, voxels{}, values{}
	{
		std::stringstream name;
		name << prefix + "." << std::setw(4) << std::setfill('0') << std::to_string(step) << ".vdb";
		filename = name.str();
	}

	std::string filename;
	int step, root;

	esint nvoxels, nvalues;
	esint* displacement;
	_Point<short>* voxels;
	float* values;
};

int SharedVolume::counter = 0;
std::mutex SharedVolume::mutex;
std::condition_variable SharedVolume::cv;
esint SharedVolume::size;
std::vector<esint> SharedVolume::sizes;

OpenVDB::OpenVDB()
: _filename(_path + _directory + _name), _step(0)
{

}

OpenVDB::~OpenVDB()
{
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
	OpenVDBWrapper wrapper(info::mesh->elements->volumeOrigin, info::mesh->elements->volumeSize, info::mesh->elements->volumeGrid);
	for (size_t di = 0; di < info::mesh->elements->data.size(); di++) { // go through all element values
		if (info::mesh->elements->data[di]->name.size()) {
			OpenVDBWrapper::FloatData *data = wrapper.addFloat(info::mesh->elements->data[di]->name);
			esint* displacement = volume->displacement;
			_Point<short>* voxels = volume->voxels;
			float* values = volume->values;
			for (int r = 0; r < info::mpi::size; ++r) {
				data->insert(volume->sizes[r], displacement, voxels, values);
				displacement += volume->size + 1;
				voxels += volume->nvoxels;
				values += volume->size;
			}
		}
	}

	wrapper.store(volume->filename.c_str());

	delete[] volume->displacement;
	delete[] volume->voxels;
	delete[] volume->values;
	delete volume; // we are in detached thread
	{
		std::lock_guard<std::mutex> lk(SharedVolume::mutex);
		--SharedVolume::counter;
	}
	SharedVolume::cv.notify_one();
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

	if (_step == 0) {
		SharedVolume::sizes.resize(info::mpi::size); // store max in the last element
		Communication::allGather(&info::mesh->elements->distribution.process.size, SharedVolume::sizes.data(), 1, MPITools::getType<esint>().mpitype, MPITools::asynchronous);
		SharedVolume::size = *std::max_element(SharedVolume::sizes.begin(), SharedVolume::sizes.end());
	}

	int node = _step % MPITools::node->across.size;
	int rank = (_step / MPITools::node->across.size) % MPITools::node->within.size;

	SharedVolume *volume = new SharedVolume(_filename, _step++, node * MPITools::node->within.size + rank);

	volume->nvoxels = info::mesh->elements->volumeIndices->datatarray().size();
	Communication::allReduce(&volume->nvoxels, nullptr, 1, MPITools::getType<size_t>().mpitype, MPI_MAX, MPITools::asynchronous);

	ivector<esint> displacement;
	ivector<_Point<short> > voxels;
	ivector<float> values;

	displacement.reserve(SharedVolume::size + 1);
	displacement.insert(displacement.end(), info::mesh->elements->volumeIndices->boundarytarray().cbegin(), info::mesh->elements->volumeIndices->boundarytarray().cend());
	displacement.resize(SharedVolume::size + 1);
	voxels.reserve(volume->nvoxels);
	voxels.insert(voxels.end(), info::mesh->elements->volumeIndices->datatarray().cbegin(), info::mesh->elements->volumeIndices->datatarray().cend());
	voxels.resize(volume->nvoxels);
	volume->nvalues = 0;
	for (size_t di = 0; di < info::mesh->elements->data.size(); di++) { // go through all element values
		if (info::mesh->elements->data[di]->name.size()) {
			++volume->nvalues;
		}
	}
	values.reserve(SharedVolume::size * volume->nvalues);
	volume->nvalues = 0;
	for (size_t di = 0; di < info::mesh->elements->data.size(); di++) { // go through all element values
		if (info::mesh->elements->data[di]->name.size()) {
			if (info::mesh->elements->data[di]->dimension == 1) {
				std::copy(info::mesh->elements->data[di]->data.cbegin(), info::mesh->elements->data[di]->data.cend(), values.end());
			} else {
				for (esint e = 0; e < info::mesh->elements->distribution.process.size; ++e) {
					float value = 0;
					for (int d = 0; d < info::mesh->elements->data[di]->dimension; ++d) {
						double &v = info::mesh->elements->data[di]->data[e * info::mesh->elements->data[di]->dimension + d];
						value += v * v;
					}
					values.push_back(std::sqrt(value));
				}
			}
			++volume->nvalues;
			values.resize(SharedVolume::size * volume->nvalues);
		}
	}
	{
		std::lock_guard<std::mutex> lk(info::mesh->voxelization.mutex);
		--info::mesh->voxelization.counter;
		info::mesh->voxelization.cv.notify_one();
	}

	if (_measure) { eslog::checkpointln("OPENVDB: DATA SERIALIZED"); }

	if (volume->root == info::mpi::rank) {
		volume->displacement = new esint[displacement.size() * info::mpi::size];
		volume->voxels = new _Point<short>[voxels.size() * info::mpi::size];
		volume->values = new float[values.size() * info::mpi::size];
	}
	Communication::gather(displacement.data(), volume->displacement, displacement.size(), MPITools::getType<esint>().mpitype, volume->root, MPITools::asynchronous);
	Communication::gather(voxels.data(), volume->voxels, voxels.size() * sizeof(_Point<short>), MPI_BYTE, volume->root, MPITools::asynchronous);
	Communication::gather(values.data(), volume->values, values.size(), MPI_FLOAT, volume->root, MPITools::asynchronous);

	if (_measure) { eslog::checkpointln("OPENVDB: DATA GATHERED"); }

	if (volume->root == info::mpi::rank) {
		{
			std::lock_guard<std::mutex> lk(SharedVolume::mutex);
			++SharedVolume::counter;
		}
		std::thread t(store, volume);
		t.detach();
	} else {
		delete volume;
	}

	if (_measure) { eslog::endln("OPENVDB: THREAD DETACHED"); }
}
