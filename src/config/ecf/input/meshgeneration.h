
#ifndef SRC_CONFIG_ECF_INPUT_MESHGENERATION_H_
#define SRC_CONFIG_ECF_INPUT_MESHGENERATION_H_

#include "config/description.h"

#include <string>

namespace espreso {

struct GMSHConfiguration: public ECFDescription {

	struct CharacteristicLength: public ECFDescription {
		double extern_from_boundary, from_points, from_curvature, min, max;

		CharacteristicLength();
	};

	CharacteristicLength characteristic_length;
	int algorithm3D, subdivisionAlgorithm, optimize;
	double stl_angle, stl_precision;

	GMSHConfiguration();
};

struct MeshGenerationConfiguration: public ECFDescription {

	GMSHConfiguration gmsh_options;

	MeshGenerationConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_MESHGENERATION_H_ */
