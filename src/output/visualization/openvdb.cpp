
#include "openvdb.h"
#include "wrappers/openvdb/w.openvdb.h"

#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"
#include "basis/containers/serializededata.h"

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
	auto volume_data = info::mesh->elements->volumeIndices;
	auto dist_array = volume_data->boundarytarray();
	auto data_array = volume_data->datatarray();

	esint dist_array_size = dist_array.size();
	esint max_dist_array_size;
	MPI_Allreduce(&dist_array_size, &max_dist_array_size, 1, MPI_UNSIGNED, MPI_MAX, info::mpi::comm);
	esint data_array_size = data_array.size();
	esint max_data_array_size;
	MPI_Allreduce(&data_array_size, &max_data_array_size, 1, MPI_UNSIGNED, MPI_MAX, info::mpi::comm);

	if (_measure) { eslog::checkpointln("OPENVDB: VOLUME SIZE REDUCED"); }

	int filling_inx = -1;
	std::vector<int> dist_vector(dist_array.begin(), dist_array.end());
	dist_vector.resize(max_dist_array_size, filling_inx);
	std::vector<int> all_dist;
	all_dist.resize(info::mpi::size * max_dist_array_size);
	MPI_Gather(dist_vector.data(), max_dist_array_size, MPI_INT, all_dist.data(), max_dist_array_size, MPI_INT, 0, info::mpi::comm);
	_Point<int> filling_point = _Point<int>(-1, -1 , -1);
	std::vector<_Point<int>> data_vector(data_array.begin(), data_array.end());
	data_vector.resize(max_data_array_size, filling_point);
	std::vector<_Point<int>> all_data;
	all_data.resize(info::mpi::size * max_data_array_size);
	MPI_Gather(data_vector.data(), max_data_array_size*3, MPI_INT, all_data.data(), max_data_array_size*3, MPI_INT, 0, info::mpi::comm);

	if (_measure) { eslog::checkpointln("OPENVDB: VOLUME GATHERED"); }

	OpenVDBWrapper wrapper;
	for (size_t di = 0; di < info::mesh->elements->data.size(); di++) { // go through all element values
		// gather value data
		auto elements_data = info::mesh->elements->data[di]->data;
		esint max_data_size = max_dist_array_size * info::mesh->elements->data[di]->dimension;		

		double filling_value = DBL_MAX;
		elements_data.resize(max_data_size, filling_value);
		std::vector<double> all_elements_data;
		all_elements_data.resize(info::mpi::size * max_data_size);
		MPI_Gather(elements_data.data(), max_data_size, MPI_DOUBLE, all_elements_data.data(), max_data_size, MPI_DOUBLE, 0, info::mpi::comm);

		if(info::mpi::rank == 0) {
			wrapper.add_grid(all_dist, all_data, all_elements_data, info::mesh->elements->data[di]->name, info::mesh->elements->data[di]->dimension);
		}
	}
	if (_measure) { eslog::checkpointln("OPENVDB: VOLUME INSERTED"); }

	if(info::mpi::rank == 0) {
		wrapper.store_grids((_filename + std::to_string(_step) + ".vdb").c_str());
		_step++;
	}
	if (_measure) { eslog::endln("OPENVDB: DATA STORED"); }
}
