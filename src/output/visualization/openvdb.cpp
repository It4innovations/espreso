
#include "openvdb.h"
#include "wrappers/openvdb/w.openvdb.h"

#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

using namespace espreso;

OpenVDB::OpenVDB()
{
	_filename = _path + _directory + _name + ".vdb";
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

	OpenVDBWrapper wrapper;
	for (size_t di = 0; di < info::mesh->elements->data.size(); di++) { // go through all element values
		wrapper.add_grid(info::mesh->elements->volumeIndices, info::mesh->elements->data[di]);
	}
	wrapper.store_grids(_filename.c_str());
}