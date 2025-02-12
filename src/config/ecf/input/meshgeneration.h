
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

struct NGLibConfiguration: public ECFDescription {

    int uselocalh;
    double maxh, minh;

    double fineness, grading;
    double elementsperedge, elementspercurve;

    int closeedgeenable;
    double closeedgefact;

    int minedgelenenable;
    double minedgelen;

    int second_order, quad_dominated;
    int optsurfmeshenable, optvolmeshenable, optsteps_3d, optsteps_2d;
    int invert_tets, invert_trigs;
    int check_overlap, check_overlapping_boundary;

    NGLibConfiguration();
};


struct MeshGenerationConfiguration: public ECFDescription {

    GMSHConfiguration gmsh_options;
    NGLibConfiguration nglib_options;

    MeshGenerationConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_MESHGENERATION_H_ */
