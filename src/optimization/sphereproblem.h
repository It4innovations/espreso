#ifndef _SRC_OPTIMIZATION_SPHEREPROBLEM_H_
#define _SRC_OPTIMIZATION_SPHEREPROBLEM_H_

#include "../config/configuration.h"

namespace espreso
{

struct SphereProblem : public ECFObject
{
    int x, y;
    // enum class SphereEnum
    // {
    //     ZERO, ONE, TWO, THREE, FOUR, FIVE
    // };

    // SphereEnum x, y;

    SphereProblem();

    double evaluate();
};

}

#endif /* _SRC_OPTIMIZATION_SPHEREPROBLEM_H_ */