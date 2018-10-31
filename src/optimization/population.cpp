#include "population.h"

#include <cstdlib>

using namespace espreso;

std::vector<std::vector<int> > PopulationGenerator::generate(OptimizationProblem& problem)
{
    std::vector<std::vector<int> > population;
    population.resize(problem.population());
    for (int i = 0; i < problem.population(); i++)
    {
        std::vector<int> specimen;
        specimen.resize(problem.dimension());
        int p = 0;
        for (auto it = problem.paramBegin(); 
        it != problem.paramEnd();
        ++it
        )
        {
            specimen[p++] = it->minValue() + rand() % ( it->maxValue() - it->minValue() );
        }
        population[i] = specimen;
    }

    return population;
}