
#ifndef SRC_CONFIG_ECF_INPUT_INPUT_H_
#define SRC_CONFIG_ECF_INPUT_INPUT_H_

#include "config/description.h"
#include "decomposition.h"
#include "transformation.h"

#include <string>
#include <map>

namespace espreso {

struct InputConfiguration: public ECFDescription {

	enum class FORMAT {
		ANSYS_CDB,
		OPENFOAM,
		ABAQUS,
		XDMF,
		ENSIGHT,
		VTK_LEGACY,
		NETGET
	};

	enum class LOADER {
		MPI,
		MPI_COLLECTIVE,
		POSIX
	};

	std::string path;
	FORMAT format;

	bool omit_midpoints, insert_midpoints;
	bool keep_material_sets;
	bool convert_database;
	double duplication_tolerance;

	LOADER loader;
	size_t stripe_size;
	int third_party_scalability_limit;

	std::map<std::string, InputTransformationConfiguration> transformations;
	DecompositionConfiguration decomposition;

	InputConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_INPUT_H_ */
