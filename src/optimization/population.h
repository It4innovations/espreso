#ifndef SRC_OPTIMIZATION_POPULATION_H_
#define SRC_OPTIMIZATION_POPULATION_H_

#include <vector>

#include "problem.h"

namespace espreso
{

class PopulationGenerator
{
public:
    static std::vector<std::vector<int> > generate(OptimizationProblem& problem);
};

}

#endif /* SRC_OPTIMIZATION_POPULATION_H_ */