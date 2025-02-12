
#ifndef SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_TRANSIENTSECONDORDERIMPLICIT_H_
#define SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_TRANSIENTSECONDORDERIMPLICIT_H_

#include "damping.h"
#include "transientfirstorderimplicit.h"
#include "config/holders/expression.h"

namespace espreso {

struct TransientSecondOrderImplicitSolverConfiguration: public ECFDescription {

    enum class METHOD {
        NEWMARK
    };

    enum class MASS_MATRIX_TYPE {
        CONSISTENT,
        DIAGONAL,
        HRZDIAGONAL
    };

    double numerical_damping;

    METHOD method;
    AutoTimeSteppingConfiguration auto_time_stepping;
    double alpha, delta, alphaF, alphaM, time_step;

    DampingConfiguration damping;

    MASS_MATRIX_TYPE mass_matrix_type;

    TransientSecondOrderImplicitSolverConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_TRANSIENTSECONDORDERIMPLICIT_H_ */
