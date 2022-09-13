
#include "openvdb.h"
#include "wrappers/openvdb/w.openvdb.h"

#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"
#include "basis/containers/allocators.h"
#include "basis/containers/serializededata.h"
#include "wrappers/mpi/communication.h"

#include <float.h>
#include <vector>

using namespace espreso;

OpenVDB::OpenVDB()
{
	_filename = _path + _directory + _name;
	_step = 0;
}

OpenVDB::~OpenVDB()
{

}

void OpenVDB::updateMesh()
{

}

void OpenVDB::updateSolution()
{
	if (info::mesh->elements->volumeIndices == nullptr) {
		printf("SET volumeIndices to store OpenVDB\n");
		return;
	}

	if (_measure) { eslog::startln("OPENVDB: STORING STARTED", "OPENVDB"); }

	// gather volume data
	size_t size[2] = { info::mesh->elements->volumeIndices->boundarytarray().size() + 1, info::mesh->elements->volumeIndices->datatarray().size() };
	Communication::allReduce(size, nullptr, 2, MPITools::getType<size_t>().mpitype, MPI_MAX);
	if (_measure) { eslog::checkpointln("OPENVDB: VOLUME SIZE REDUCED"); }

	ivector<esint> localDist, globalDist;
	ivector<_Point<short> > localData, globalData;
	ivector<float> localValue(size[0]), globalValues;

	localDist.reserve(size[0]);
	localDist.push_back(info::mesh->elements->volumeIndices->boundarytarray().size() - 1);
	localDist.insert(localDist.end(), info::mesh->elements->volumeIndices->boundarytarray().cbegin(), info::mesh->elements->volumeIndices->boundarytarray().cend());
	localDist.resize(size[0]);

	localData.reserve(size[1]);
	localData.insert(localData.end(), info::mesh->elements->volumeIndices->datatarray().cbegin(), info::mesh->elements->volumeIndices->datatarray().cend());
	localData.resize(size[1]);
	if (info::mpi::rank == 0) {
		globalDist.resize(info::mpi::size * size[0]);
		globalData.resize(info::mpi::size * size[1]);
		globalValues.resize(info::mpi::size * size[0]);
	}

	MPI_Gather(localDist.data(), size[0], MPITools::getType<esint>().mpitype, globalDist.data(), size[0], MPITools::getType<esint>().mpitype, 0, info::mpi::comm);
	MPI_Gather(localData.data(), 3 * size[1], MPITools::getType<short>().mpitype, globalData.data(), 3 * size[1], MPITools::getType<short>().mpitype, 0, info::mpi::comm);

	if (_measure) { eslog::checkpointln("OPENVDB: VOLUME GATHERED"); }

	OpenVDBWrapper wrapper;
	for (size_t di = 0; di < info::mesh->elements->data.size(); di++) { // go through all element values
		if (info::mesh->elements->data[di]->name.size()) {
			if (info::mesh->elements->data[di]->dimension == 1) {
				std::copy(info::mesh->elements->data[di]->data.cbegin(), info::mesh->elements->data[di]->data.cend(), localValue.begin());
			} else {
				for (esint e = 0; e < info::mesh->elements->distribution.process.size; ++e) {
					float value = 0;
					for (int d = 0; d < info::mesh->elements->data[di]->dimension; ++d) {
						double &v = info::mesh->elements->data[di]->data[e * info::mesh->elements->data[di]->dimension + d];
						value += v * v;
					}
					localValue[e] = std::sqrt(value);
				}
			}
		}
		MPI_Gather(localValue.data(), size[0], MPI_FLOAT, globalValues.data(), size[0], MPI_FLOAT, 0, info::mpi::comm);

		if(info::mpi::rank == 0) {
			wrapper.add_grid(size[0], size[1], globalDist.data(), globalData.data(), globalValues.data(), info::mesh->elements->data[di]->name);
		}
	}
	if (_measure) { eslog::checkpointln("OPENVDB: VOLUME INSERTED"); }

	if(info::mpi::rank == 0) {
		wrapper.store_grids((_filename + std::to_string(_step) + ".vdb").c_str());
		_step++;
	}
	if (_measure) { eslog::endln("OPENVDB: DATA STORED"); }
}
