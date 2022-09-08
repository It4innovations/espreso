

#include "ensightgold.volume.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/ecfinfo.h"
#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <functional>

using namespace espreso;

EnSightGoldVolume::EnSightGoldVolume()
: _step(-1)
{
	_geometry = _directory + _name + ".geo";
	_fixedDataPath = _directory;
}

EnSightGoldVolume::~EnSightGoldVolume()
{

}

void EnSightGoldVolume::updateMesh()
{
	if (isRoot()) {
		casefile();

		_writer.format();
		_writer.description("EnSight Gold geometry format");
		_writer.description("----------------------------");
		_writer.description("node id off");
		_writer.description("element id off");
		_writer.description("part");
		_writer.int32ln(1);
		_writer.description("3D uniform grid");
		_writer.description("block uniform");

		_writer.int32(info::mesh->elements->volumeGrid.x);
		_writer.int32(info::mesh->elements->volumeGrid.y);
		_writer.int32(info::mesh->elements->volumeGrid.z);
		_writer.ln();

		_writer.float32ln(info::mesh->elements->volumeOrigin.x);
		_writer.float32ln(info::mesh->elements->volumeOrigin.y);
		_writer.float32ln(info::mesh->elements->volumeOrigin.z);

		_writer.float32ln(info::mesh->elements->volumeSize.x / info::mesh->elements->volumeGrid.x);
		_writer.float32ln(info::mesh->elements->volumeSize.y / info::mesh->elements->volumeGrid.y);
		_writer.float32ln(info::mesh->elements->volumeSize.z / info::mesh->elements->volumeGrid.z);
	}

	_writer.groupData();
	_writer.commitFile(_path + _geometry);
	_writer.reorder();
	_writer.write();
}

void EnSightGoldVolume::updateSolution()
{
	if (_measure) { eslog::startln("ENSIGHT VOLUME: STORING STARTED", "ENSIGHT VOLUME"); }
	switch (step::outstep.type) {
	case step::TYPE::TIME:      _times.push_back(step::outtime.current); break;
	case step::TYPE::FREQUENCY: _times.push_back(step::outfrequency.current); break;
	case step::TYPE::FTT: break;
	}

	int datasize = 0;
	for (size_t di = 0; di < info::mesh->elements->data.size(); di++) {
		datasize += info::mesh->elements->data[di]->dimension;
	}
	if (info::mesh->elements->data.size() != _variables.size()) { // add support for unnamed data
		_variables.clear();
		for (size_t i = 0; i < info::mesh->elements->data.size(); i++) {
			NamedData *data = info::mesh->elements->data[i];
			if (data->dataType == NamedData::DataType::VECTOR) {
				_variables.push_back("vector per node: 1 " + data->name + "   " + _directory + data->name + ".****" );
			}
			if (data->dimension == 1) {
				_variables.push_back("scalar per node: 1 " + data->name + "   " + _directory + data->name + ".****" );
			}
		}
	}

	const auto &grid = info::mesh->elements->volumeGrid;
	size_t gridsize = grid.x * grid.y * grid.z;
	size_t chunk = gridsize / info::mpi::size + ((gridsize % info::mpi::size) ? 1 : 0);
	size_t coffset = chunk * info::mpi::rank;
	size_t mychunk = std::min(gridsize - coffset, chunk);

	ivector<int> sBuffer, rBuffer;
	std::vector<esint> offset(info::mpi::size + 1, 2);
	for (auto voxels = info::mesh->elements->volumeIndices->begin(); voxels != info::mesh->elements->volumeIndices->end(); ++voxels) {
		for (auto voxel = voxels->begin(); voxel != voxels->end(); ++voxel) {
			offset[*voxel / chunk] += 1 + datasize;
		}
	}
	utils::sizesToOffsets(offset);
	sBuffer.resize(offset.back());

	for (int r = 0; r < info::mpi::size; ++r) {
		sBuffer[offset[r]] = offset[r + 1] - offset[r];
		++offset[r];
		sBuffer[offset[r]++] = r;
	}

	esint e = 0;
	for (auto voxels = info::mesh->elements->volumeIndices->begin(); voxels != info::mesh->elements->volumeIndices->end(); ++voxels, ++e) {
		for (auto voxel = voxels->begin(); voxel != voxels->end(); ++voxel) {
			int target = *voxel / chunk;
			sBuffer[offset[target]++] = *voxel;
			for (size_t di = 0; di < info::mesh->elements->data.size(); di++) {
				NamedData *data = info::mesh->elements->data[di];
				for (int d = 0; d < data->dimension; ++d) {
					float value = data->data[data->dimension * e + d];
					memcpy(sBuffer.data() + offset[target], &value, sizeof(float));
					++offset[target];
				}
			}
		}
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer, 0, info::mpi::size, MPITools::asynchronous)) {
		eslog::error("cannot exchange voxels\n");
	}
	utils::clearVector(sBuffer);

	if (_measure) { eslog::checkpointln("ENSIGHT VOLUME: DATA CHUNKED"); }

	std::vector<float*> output(info::mesh->elements->data.size());
	for (size_t di = 0; di < info::mesh->elements->data.size(); di++) {
		NamedData *data = info::mesh->elements->data[di];
		if (isRoot()) {
			_writer.description(data->name);
			_writer.description("part");
			_writer.int32ln(1);
			_writer.description("coordinates");
		}

		// TODO: multi-dimensional data
		output[di] = (float*)_writer.enlarge(mychunk * info::mesh->elements->data[di]->dimension * sizeof(float));
		_writer.groupData();

		std::stringstream file;
		file << _path + _directory + info::mesh->elements->data[di]->name + "." << std::setw(4) << std::setfill('0') << _times.size() + step::outduplicate.offset;
		_writer.commitFile(file.str());
	}

	for (esint r = 0, offset = 0, next = 0; r < info::mpi::size; ++r) {
		next += rBuffer[offset++];
		offset++; // source
		while (offset < next) {
			esint voxel = rBuffer[offset++] - coffset;
			for (size_t di = 0; di < info::mesh->elements->data.size(); di++) {
				NamedData *data = info::mesh->elements->data[di];
				for (int d = 0; d < data->dimension; ++d) {
					output[di][d * mychunk + voxel] = reinterpret_cast<float&>(rBuffer[offset]);
					++offset;
				}
			}
		}
	}

	if (_measure) { eslog::checkpointln("ENSIGHT VOLUME: DATA INSERTED"); }

	_writer.directWrite();
	if (isRoot()) {
		casefile();
	}

	if (_measure) { eslog::endln("ENSIGHT VOLUME: DATA STORED"); }
}

std::string EnSightGoldVolume::dataname(const NamedData *data, int d)
{
	if (data->dataType == NamedData::DataType::VECTOR) {
		return data->name;
	}
	if (data->dimension == 1) {
		return data->name;
	}
	return data->name + data->suffix(d);
}

void EnSightGoldVolume::casefile()
{
	std::ofstream os(_path + _name + ".case");

	os << "\n# output from ESPRESO library (espreso.it4.cz)";
	os << "\n#";
	os << "\n";
	os << "\nFORMAT";
	os << "\ntype: ensight gold";
	os << "\n";
	os << "\nGEOMETRY";
	os << "\nmodel: " << _geometry;
	os << "\n";
	os << "\nVARIABLE";
	os << "\n";
	if (_times.size()) {
		for (size_t i = 0; i < _variables.size(); ++i) {
			os << "\n" << _variables[i];
		}
		os << "\n";
		os << "\nTIME";
		os << "\ntime set:               1";
		os << "\nnumber of steps:        " << _times.size();
		os << "\nfilename start numbers: 1";
		os << "\nfilename increment:     1";
		os << "\ntime values:";
		for (size_t i = 0; i < _times.size(); ++i) {
			if (i % 10 == 0) {
				os << "\n" << _times[i];
			} else {
				os << " " << _times[i];
			}
		}
	}
}




