
#ifndef SRC_CONFIG_ECF_INPUT_INPUT_H_
#define SRC_CONFIG_ECF_INPUT_INPUT_H_

#include "config/description.h"
#include "contactinterface.h"
#include "decomposition.h"
#include "meshgeneration.h"
#include "transformation.h"

#include <string>
#include <map>

namespace espreso {

struct ClippingBox: public ECFDescription {
	double min[3], max[3];
	bool apply;

	ClippingBox();
};

struct InputConfiguration: public ECFDescription {

	enum class FORMAT {
		ANSYS_CDB,
		OPENFOAM,
		ABAQUS,
		XDMF,
		ENSIGHT,
		VTK_LEGACY,
		NETGET,
		GMSH,
		NGLIB
	};

	enum class LOADER {
		MPI,
		MPI_COLLECTIVE,
		POSIX
	};

	std::string path;
	FORMAT format;

	ClippingBox clipping_box;

	bool omit_midpoints, insert_midpoints;
	bool omit_face_sets;
	bool keep_material_sets;
	bool convert_database;
	double duplication_tolerance;

	LOADER loader;
	size_t stripe_size;
	int third_party_scalability_limit;

	std::map<std::string, InputTransformationConfiguration> transformations;
	DecompositionConfiguration decomposition;
	MeshGenerationConfiguration generation;
	std::map<std::string, ContactInterfaceConfiguration> contact_interfaces;

	InputConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_INPUT_H_ */
