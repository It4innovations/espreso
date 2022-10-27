
#ifndef SRC_CONFIG_ECF_INPUT_INPUT_H_
#define SRC_CONFIG_ECF_INPUT_INPUT_H_

#include "config/description.h"
#include "contactinterface.h"
#include "decomposition.h"
#include "meshgeneration.h"
#include "noderegion.h"
#include "transformation.h"

#include <string>
#include <map>

namespace espreso {

struct ClippingBox: public ECFDescription {
	double min[3], max[3];
	bool apply;

	ClippingBox();
};

struct EnSightInputConfiguration: public ECFDescription {
	int timeset;
	std::string time, variables;

	EnSightInputConfiguration();
};

struct OpenFOAMInputConfiguration: public ECFDescription {
	bool direct_load;

	std::string time, variables;

	OpenFOAMInputConfiguration();
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
		NEPER,
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

	EnSightInputConfiguration ensight;
	OpenFOAMInputConfiguration openfoam;

	ClippingBox clipping_box;

	bool omit_midpoints, insert_midpoints;
	bool omit_face_sets;
	bool keep_material_sets;
	bool convert_database;
	bool compute_bodies;
	double duplication_tolerance;

	bool insert_orientation;

	LOADER loader;
	size_t stripe_size;
	int third_party_scalability_limit;

	std::map<std::string, InputTransformationConfiguration> transformations;
	std::map<std::string, InputNodeRegionConfiguration> node_regions;
	DecompositionConfiguration decomposition;
	MeshGenerationConfiguration generation;
	std::map<std::string, ContactInterfaceConfiguration> contact_interfaces;

	InputConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_INPUT_H_ */
