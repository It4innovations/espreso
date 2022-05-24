
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
	wrapper.setTransformation();
	for (size_t di = 0; di < info::mesh->elements->data.size(); di++) { // go through all element values
		wrapper.store(_filename.c_str(), info::mesh->elements->volumeIndices, info::mesh->elements->data[di]);
	}
}